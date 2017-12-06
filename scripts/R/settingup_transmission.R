message(
  paste0(
    "##############################################################################\n",
    "########################## Setting up transmission ###########################\n",
    "##############################################################################\n"))


# ============================ Packages =======================================
require(knitr)
require(ggplot2)
require(magrittr)
require(tidyverse)
if(!("package:tidyverse" %in% search())){
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(readr)
  require(rlang)
}
require(HIVEr)
require(extrafont)
require(wesanderson)
require(doMC)
doMC::registerDoMC(cores=8)
options(warn=1)

set.seed(42) # Set seed so randomization is reproducible
# ============================ Functions =========================================
write_to_summary<-function(line_pattern,value){
  file = readLines("../../results/results.table.tsv")
  line_pattern_regex = paste0("^",line_pattern)
  line = grep(line_pattern_regex,file)
  file[line] = paste0(line_pattern,"\t",value)
  writeLines(file,"../../results/results.table.tsv")
}
# =========================== Reading in the data ===============================
message("Reading in the data")

meta<-read_csv("../../data/reference/all_meta.sequence_success.csv")
qual<-read_csv("../../data/processed/secondary/qual.snv.csv",
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               ))

no_cut_qual <- read_csv("../../data/processed/secondary/no_freq_cut.qual.snv.csv",
                        col_types = list(
                          ENROLLID= col_character(),
                          SPECID = col_character(),
                          LAURING_ID = col_character(),
                          Id = col_character()
                        ))

  
trans_pairs<-read_csv("../../data/processed/secondary/transmission_pairs.csv")



# =========================== Initial processing ===============================
message("Initial processing")
useful_trans_pairs<-subset(trans_pairs,
                           snv_qualified_pair==T & valid==T & quality_distance==T,
                           select=-c(gc_ul1,gc_ul2))
write_to_summary("Sequence validated pairs:",nrow(useful_trans_pairs))

## Now lets add the flipped pairs when applicatble

flipped<-subset(useful_trans_pairs,double==T)

flipped<-plyr::rename(flipped, c("ENROLLID1"="ENROLLID2","ENROLLID2"="ENROLLID1",
                                 "onset1"="onset2","onset2"="onset1",
                                 "sequenced1"="sequenced2","sequenced2"="sequenced1",
                                 "snv_qualified1"="snv_qualified2","snv_qualified2"="snv_qualified1"))
flipped$pair_id=flipped$pair_id+0.5
useful_trans_pairs<-rbind(useful_trans_pairs,flipped)

write_to_summary("Duel directionality:",nrow(flipped))

# ==================================== Selecting a SPECID =======================
# Getting the SPECID for each pair.
# 
# Here we use the following criteria - note this may result in different SPECID 
# than used in the L1-norm measurements.
# 
# 
# 1) Sample closest to transmission that are on the "right side" of
# tranmission (before transmission for donor if possible)
# 
# 2) Titer is the tie breaker when applicable.
# 
# 3) as in the get close docs. We take the sample with iSNV if the 
# ideal one doesn't have any
# Here we are going to count the minority iSNV in each sample.

message("Selecting SPECID")

min_qual<-subset(qual,freq.var<0.5)
min.count<-plyr::ddply(min_qual,~SPECID,summarize,iSNV=length(unique(mutation)))
meta<-left_join(meta,min.count)

meta$iSNV[is.na(meta$iSNV)]<-0

# Select the SPECID
useful_trans_pairs<-useful_trans_pairs %>% rowwise()%>%
  mutate(SPECID1=get_close(meta,transmission,ENROLLID1,"donor"),
         SPECID2=get_close(meta,transmission,ENROLLID2,"recipient"))

# Add the dates for collection from the meta data file.
useful_trans_pairs<-mutate(useful_trans_pairs,
                           collect1=meta$collect[match(SPECID1,meta$SPECID)],
                           collect2=meta$collect[match(SPECID2,meta$SPECID)]) 

# verify no mixed infections
putative_mixed<-c("HS1530", "M54062" ,"MH8125", "MH8137" ,"MH8156" ,"MH8390")
pmixed_tp <- putative_mixed[which(putative_mixed %in% 
                       c(useful_trans_pairs$SPECID2,useful_trans_pairs$SPECID1))]
useful_trans_pairs<-subset(useful_trans_pairs,
                           !(SPECID1 %in% putative_mixed) & 
                             !(SPECID2 %in% putative_mixed))
write_to_summary("PT mixed infection:",paste(pmixed_tp,collapse = ", "))



# ========================  Compare the frequency of mutations =================
message("Comparing frequencies")
trans_freq<-plyr::adply(useful_trans_pairs,1,function(x){
  get_freqs(c(x$SPECID1,x$SPECID2),qual)},
  .parallel = T)
write.csv(trans_freq,file = "../../data/processed/secondary/trans_freq.csv")

# Reduce to sites that are polymorphic in the donor.
trans_freq.comp<-polish_freq(trans_freq,freq1,0.02)
trans_freq.comp$found=trans_freq.comp$freq2>0.02 # was it found in the second sample
write.csv(x = trans_freq.comp,
          file = "../../data/processed/secondary/transmission_pairs_freq.poly.donor.csv")

# ==========  Compare the frequency of mutations no frequency cutoff ==========
message("Comparing frequencies without a frequency cut off")
no_cut_trans_freq<-plyr::adply(useful_trans_pairs,1,function(x){
  get_freqs(c(x$SPECID1,x$SPECID2),no_cut_qual)},
  .parallel = T)
write.csv(no_cut_trans_freq,file = "../../data/processed/secondary/no_cut_trans_freq.csv")

# Reduce to sites that are polymorphic in the donor.
no_cut_trans_freq.comp<-polish_freq(no_cut_trans_freq,freq1,0)
no_cut_trans_freq.comp$found=no_cut_trans_freq.comp$freq2>0 # was it found in the second sample
write.csv(x = no_cut_trans_freq.comp,
         file =  "../../data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv")

# ================================= Community pairs ===========================
message("Community pair comparisions")

possible_pairs<-read_csv("../../data/processed/secondary/possible.pairs.dist.csv") # All possible pairs (1 SPECID/person)
community_pairs<-subset(possible_pairs,Household==F) # not household pairs
community_pairs$pair_id=1:nrow(community_pairs)
community_pairs.freq<-plyr::adply(community_pairs,1,function(x) {
  get_freqs(c(x$SPECID1,x$SPECID2),qual)
  },
  .parallel = T)

write.csv(community_pairs.freq,file = "../../data/processed/secondary/community_pairs.freq.csv")

# only polymorphic sites in sample 1 
community_pairs.freq.comp<-polish_freq(community_pairs.freq,freq1,0.02) 
# was it found in the second sample
community_pairs.freq.comp$found=community_pairs.freq.comp$freq2>0.02 

write.csv(community_pairs.freq.comp,
         file =  "../../data/processed/secondary/community_pairs_freq.poly.donor.csv")






