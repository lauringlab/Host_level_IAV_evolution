message(
  paste0(
    "##############################################################################\n",
    "############################### Distance work ################################\n",
    "##############################################################################\n"))# ============================ Packages =======================================
require(magrittr)
require(tidyverse)
require(HIVEr)
require(multidplyr)
require(testthat)

if(!("package:tidyverse" %in% search())){
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(readr)
  require(rlang)
}
require(doMC)
doMC::registerDoMC(cores=8)
options(warn=1)
# =========================== Reading in the data ===============================

meta<-read_csv("./data/reference/all_meta.sequence_success.csv")
qual<-read_csv("./data/processed/secondary/qual.snv.csv",
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all 
stopifnot(min(qual$freq.var)>0.02)

# ============================== meta one =======================================
# Here we get one SPECID for each ENROLLID. We choose the highest quality SPECID
# avaliable. see documentation for details.

meta_one<-only_one(meta)
# Now we will make all possible pair comparisions within a season and strain. 
# Note - because we are using the meta_one data frame we are not accounting for
# intrahost dynamics (in season 2014-2015). We have set up meta_one to give us 
# the highest quality iSNV data.

# ======================= Getting all SPECID comparisions =====================

all_possible_pairs<- meta_one %>%
  filter(snv_qualified==T) %>%
  group_by(season,pcr_result) %>%
  do({expand.grid(.$SPECID,.$SPECID)}) %>%
  ungroup()# there is only 1 SPECID/ ERNOLLID in meta_one

# Expand grid leaves the last columns names as "V1" and "V2"
names(all_possible_pairs)<-c("season","pcr_result","SPECID1","SPECID2")

# Expands gives us all possible comparisions and that's too much. We take care of that
# below

# This ensures that
# 1) we are not comparing a sample against itself
# 2) both samples were well sequenced and had a chance of calling variants.
# Recall meta_one contains the highest quality isolate for each person
all_possible_pairs<-subset(all_possible_pairs,SPECID1!=SPECID2 & 
                             SPECID1 %in% meta_one$SPECID[meta_one$snv_qualified==T] & 
                             SPECID2 %in% meta_one$SPECID[meta_one$snv_qualified==T] ) 


## One of the draw backs of expand grid is we get all possibilities so there are 2 pairs for each
# combination of specid (bidirectional). We only want 1.

# This requires that the second sample was taken from the person who was sick second. 
# This matches how we set pairs in tranmission. the EnrollID is the tie-braker.

# First we need to add ENROLLID and onset  columns
all_possible_pairs<-mutate(all_possible_pairs,
                           # Get ENROLLID from meta_one
                           ENROLLID1=meta_one$ENROLLID[match(x = SPECID1,meta_one$SPECID)],
                           ENROLLID2=meta_one$ENROLLID[match(x = SPECID2,meta_one$SPECID)],
                           # Get *difference* in time of onset from meta_1
                           # And the difference in collection times
                           time_onset=meta_one$onset[match(x = SPECID2,meta_one$SPECID)] - 
                             meta_one$onset[match(SPECID1,meta_one$SPECID)],
                           time_collect=meta_one$collect[match(x = SPECID2,meta_one$SPECID)] - 
                             meta_one$collect[match(SPECID1,meta_one$SPECID)]
)

write.csv(all_possible_pairs,file="./data/processed/secondary/every_possible_pair.csv")

# only those with SPECID2 sick later than or equal to SPECID1
all_possible_pairs<-subset(all_possible_pairs,time_onset>=0) 

# Now for those with a difference of 0. 
zeros_p<-subset(all_possible_pairs,time_onset==0)

# This is the tie break in the transmission pair search.
zeros_pass<- zeros_p %>% 
  filter(as.character(ENROLLID1)<as.character(ENROLLID2)) 

# now get rid of the 0's from the all_possible_pairs data frame
fine_p<-subset(all_possible_pairs,time_onset>0) 

# add back the zeros that work
possible_pairs<-rbind(fine_p,zeros_pass) 

# ======================== Distance calculations ===================================
# I know how to use plyr in parallel so we'll use it here
possible_pairs.dist<-plyr::adply(possible_pairs,1,
                                 summarize,L1_norm=dist_tp(c(SPECID1,SPECID2),
                                                           snv = qual),
               .parallel = TRUE) 

# if you wanted to do dplyr here it is.
# possible_pairs.dist <- possible_pairs %>% rowwise() %>%
#  mutate(L1_norm= dist_tp(c(SPECID1,SPECID2),snv = qual)) %>% ungroup()

# ====================== Working with transmission pairs ========================

# Get the houses for each pair. 
possible_pairs.dist<-mutate(possible_pairs.dist,
                            HOUSE_ID1=meta_one$HOUSE_ID[match(x = SPECID1,meta_one$SPECID)],
                            HOUSE_ID2=meta_one$HOUSE_ID[match(x = SPECID2,meta_one$SPECID)],
                            Household = HOUSE_ID1==HOUSE_ID2)

# Now Id the valid transmission pairs in the possible_pairs.

# These are the valid transmission pairs that match the epidemiological definition.
all_pairs.tp <- getting_tp(meta_one)

all_pairs.tp$pair_id<-1:nrow(all_pairs.tp)

all_pairs.tp<-ungroup(all_pairs.tp)
valid_pairs<-subset(all_pairs.tp,valid==T)

# Here we add the valid designation for pairs (each row is a pair) where
# the recipient (ENROLLID2) is in 
possible_pairs.dist<- left_join(possible_pairs.dist,
                                select(all_pairs.tp,season,pcr_result,
                                       ENROLLID1,ENROLLID2,valid))
# Now we correct for the pairs that had no chance of a valid transmission pair
possible_pairs.dist$valid[is.na(possible_pairs.dist$valid)]<-FALSE

valid_sequenced_pairs<-all_pairs.tp %>% filter(snv_qualified_pair==T) %>% select(valid)
stopifnot(sum(possible_pairs.dist$valid)==sum(valid_sequenced_pairs$valid))



# ====================== Distance cutoffs =====================================
# Now we need to add a quality cut off based on the distance distributions
# No household pairs included
cutoffs<-tibble(L1_norm=quantile(possible_pairs.dist$L1_norm[possible_pairs.dist$Household==F],
                                probs = seq(0,1,0.05))) # get the percentiles for the community pairs
names(cutoffs)<-"L1_norm"
cutoffs$threshold<-seq(0,1,0.05)
# How many valid pairs exist below each 5% threshold
cutoffs<-cutoffs %>% rowwise() %>%
  mutate(valid_pairs=
           nrow(possible_pairs.dist[(possible_pairs.dist$valid==T & 
                                       possible_pairs.dist$L1_norm<L1_norm),]))

message(paste0("The cutoff will be ", cutoffs$L1_norm[cutoffs$threshold==0.05]))
possible_pairs.dist<-mutate(possible_pairs.dist,
                            quality_distance=L1_norm<cutoffs$L1_norm[cutoffs$threshold==0.05])

# Add distance cutoff data to all_pairs.tp

all_pairs.tp<-left_join(all_pairs.tp,
                    select(possible_pairs.dist,ENROLLID1,ENROLLID2,quality_distance))

# We only want quality distances were we had a chance to measure the distance
all_pairs.tp$quality_distance[all_pairs.tp$snv_qualified_pair==F]<-NA
# ============================== Save the output ===============================
write.csv(possible_pairs.dist,file = "./data/processed/secondary/possible.pairs.dist.csv")
write.csv(all_pairs.tp,file = "./data/processed/secondary/transmission_pairs.csv")
write.csv(meta_one,file = "./data/processed/secondary/meta_one.sequence.success.csv")
