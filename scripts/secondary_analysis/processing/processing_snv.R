# Better comments and partitions

message(
  paste0(
"##############################################################################\n",
"############################## Processing snv ###############################\n",
"##############################################################################\n"))
# ============================ Packages =======================================
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
options(warn=1)
# =========================== Reading in the data ===============================

# Read in all the files we'll need 
variants_csv<-c("./data/processed/HK_1/all.variants.csv",
                "./data/processed/HK_2/all.variants.csv",
                "./data/processed/HK_6/all.variants.csv",
                "./data/processed/HK_7/all.variants.csv",
                "./data/processed/HK_8/all.variants.csv",
                "./data/processed/cali09/all.variants.csv",
                "./data/processed/cali09_2/all.variants.csv",
                "./data/processed/victoria/all.variants.csv",
                "./data/processed/victoria_2/all.variants.csv",
                "./data/processed/perth/all.variants.csv",
                "./data/processed/perth_2/all.variants.csv")
# This will overwrite the run column to be consistant with the coverage file.

variants<-read_rbind(variants_csv,n=5,cols = list(Id=col_character()))

meta<-read_csv("./data/reference/all_meta.csv")
cov_sample <- read_csv("./data/processed/secondary/average_coverages.csv",
                       col_types = list(Id=col_character(),
                      LAURING_ID=col_character()))


# ========================= Updating the meta data file ==========================

# The season were not in a uniform format
meta$pcr_result[meta$pcr_result=="H3N2"]<-"A/H3N2"
meta$pcr_result[meta$pcr_result=="A/H3N3 & B/YAM"] <-"A/H3N2" # this one didn't qualify for iSNV identification
meta$season[meta$season=="10-11"]<-"2010-2011"
meta$season[meta$season=="11-12"]<-"2011-2012"
meta$season[meta$season=="12-13"]<-"2012-2013"
meta$season[meta$season=="13-14"]<-"2013-2014"   

# some samples had poor coverage and so did not qualify for variants calling
poor_coverage_Id<-subset(cov_sample,coverage<1000,select=c(Id,LAURING_ID,run))
meta<-mutate(meta,
             snv_qualified=(gc_ul>1e3 & 
                              sequenced==T & 
                              !(LAURING_ID %in% poor_coverage_Id$LAURING_ID)))

# using LAURING_ID is ok as if a sample was sequenced twice then it needs to have
# good coverage in both samples by definition.





# ========================= Cleaning up the variant data =====================

# Id refers to the sequenced sample - 
# LAURING_ID is analogous to the SPECID - the actual nose and throat sample we
# processed in the lab. There can be mulitple Id for the same LAURING_ID
# 
# Here we separate the Id into the LAURING_ID and duplicate where needed.
# this leads to SPECID==LAURING_ID for season 2014-2014 and dup of A or B.
# For the other seasons this leaves Id=LAURING_ID and a dup of NA. Other seasons
# Id do not contain duplicate labeling. That information is held in the sequencing
# run. i.e. duplicates were sequenced on separate runs. (this is true for 2014-2015
# as well)


# some Ids include a decimal - this removes that
variants <- variants %>%
  do(correct_id(.,Id))

variants <- variants %>% 
  tidyr::separate(Id,into=c("LAURING_ID","dup"),
                  sep="_",remove=F,fill="right")

# Now we add the meta data to the variant file
variants<-left_join(variants,meta,by="LAURING_ID")

# If the LAUIRNG_ID is not in the meta file the sample will 
# have a lot of NA columns. gc_ul is just one of those. This should only remove
# plasmid controls. We check that with the messages below.

extra<-subset(variants,is.na(gc_ul))
variants<-subset(variants,!(is.na(gc_ul)))

message(paste0("We are removing sample ",
               unique(extra$Id), " becasue they weren't found in the meta data.\n"))

# This makes sure every sample that was sequenced twice was sequenced on separate runs 
# - This is (not any more) assumed in the qual function below. only HK runs have duplicates. 
# The others were separated by time so I know they were sequenced on different run.

wrong <- variants %>% group_by(LAURING_ID) %>%
  summarize(runs = length(unique(run)),
            dup_labels=length(unique(dup))) %>%
  filter(runs<dup_labels)
# This is the case where the sample was sequenced on fewer runs than duplicates 
# i.e. both duplicates on 1 run
stopifnot(nrow(wrong)==0)

# Id refers to the sequenced sample - 
# LAURING_ID refers to the SPECID - the actual nose and through sample

# we will remove samples that lack sufficient coverage.
variants<-left_join(variants,cov_sample,by = c("Id", "run", "LAURING_ID", "dup"))
variants <-rename(variants,sample_coverage = coverage )
stopifnot(any(!is.na(variants$sample_coverage)))


poor_coverage.variants<-variants %>% filter(sample_coverage<1000)


variants <- variants %>% filter(sample_coverage>=1000)


######## limit the work ########

# Let's remove the sites that only have one allele present in a season.

diverse_variants<-diverse_sites(variants,1,season,pcr_result,pos,chr)

# Ensure this doesn't eliminate samples

stopifnot(all(unique(paste0(variants$Id,variants$run))==unique(paste0(diverse_variants$Id,diverse_variants$run))))
# join alleles to variants
# remove those with only 1 allele

#--------------------- Variant calling ---------------------------------------#



############### THE LONG LOOP ##############

qual <- diverse_variants %>% group_by(LAURING_ID) %>%
   do(quality(.))




qual$class_factor=NA
qual$class_factor[grep("Noncoding",qual$Class)]<-"Noncoding"
qual$class_factor[grep('Syn', qual$Class)]<-"Synonymous"
qual$class_factor[grep('Nonsyn',qual$Class)]<-"Nonsynonymous" 
# if it is nonsynonymous in any OR we will catch it

qual$class_factor<-as.factor(qual$class_factor)
qual<-subset(qual,class_factor!="Noncoding") # Eliminate the noncoding

qual<-diverse_sites(qual,1,season,pcr_result,pos,chr) 
# This is valid as we have not removed any sites so fixed differences within a season will be
# maintained. We only remove sites with one allele shared in all samples. In
# no way will these sites enter into our calculations. By definition they are 
# present at 100% in all samples (otherwise there would be other alleles).
# (See monomorphic below - this command does remove some sites with freq.var<0.90 but in
# our analysis that is equivalent to setting the frequency to 1 since it is shared by all samples
# and there are no minority alleles.)

no_freq_cut<-qual

################## Frequency cut off ############################################
qual<- qual %>% filter(freq.var>0.02)
# This line is replaced by those below setting monorphic sites to frequencies of 1.
#qual$freq.var[qual$freq.var>0.98]<-1 # a corolary of the rule above
# This removes the consensus bases where there are no polymorphism 
# (at least at a 2% threshold)

# Here we set the frequency of monomorphic sites to 1. These are sites where there is 
# only one allele. I have investigated these sites by hand. There are cases where 
# their frequencies are <1. I have looked at all sites with such frequencies between 50-90%.
# (There are none between 1-50%.) In every case there were false positive 
# minor variants that failed either due to artifacts near the end of the read which were either
# not called by deepSNV (likely due to stran bias) or because of our benchmarked quality controls.
# Some were also indels - we don't look at indels in this work and never have. Thus we are right
# to only report 1 allele; however, the frequencies are based on the raw data which encludes erroneous
# reads and false positives. This line corrects these frequencies. This work is saved in 
# results/notebook/2017_11_13.Rmd through 2017_11_16.Rmd. It relies on the original quality 
# data frame which did not have the diverse sites removed and employed a frequency cut off of 2%. The 
# diverse sites command also fixes this issue as many of these sites are monomorphic with a season.
# That work is verified in 2017_11_16.Rmd. 

qual<-monomorphic(qual,SPECID,season,chr,pos,pcr_result)
no_freq_cut<-monomorphic(no_freq_cut,SPECID,season,chr,pos,pcr_result)


# we can also removed the fixed reference alleles as these are accounted for in the equal compare 
# fucntion and there can be no minor variants here.
qual<-subset(qual,!(ref==var & freq.var==1)) 

no_freq_cut<-subset(no_freq_cut,!(ref==var & freq.var==1)) 

# now we can cut off the loci that now only have 1 allele per season per strain.
#qual<-diverse_sites(qual,1,season,pcr_result,pos,chr) 
## This is not valid because we have removed reference sites and so we will remove
## fixed differences. I leave this as a warning. this must be done above. prior to
## removing sites.

###################### Saving the data #########################################
write.csv(x=no_freq_cut,file="./data/processed/secondary/no_freq_cut.qual.snv.csv")
write.csv(x=qual,file="./data/processed/secondary/qual.snv.csv")
write.csv(x = meta,file = "./data/reference/all_meta.sequence_success.csv")
