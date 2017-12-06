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
options(warn =1)
# =========================== Reading in the data ===============================
message(
  paste0(
    "##############################################################################\n",
    "############################## Intrahost setup ###############################\n",
    "##############################################################################\n"))
meta<-read_csv("../../data/reference/all_meta.sequence_success.csv")
qual<-read_csv("../../data/processed/secondary/qual.snv.csv",
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all 

# Snv.intra is the data frame that contains meta data on intrahost pairs 
# that qualified for snv identificiation
snv.intra<-subset(meta,snv_qualified==T) %>%
  get_double() # Get those samples that come from longitundal pairs

snv.intra <- dplyr::select(snv.intra,c(HOUSE_ID,ENROLLID,season,pcr_result,
                                       SPECID,collect,gc_ul,onset,home_collected))

intra<-short_pairs(snv.intra,c("SPECID","collect","gc_ul","home_collected"),
                   HOUSE_ID,ENROLLID,season,pcr_result,onset)

intra_freq<-plyr::adply(intra,1,
                        function(x) get_freqs(c(x$SPECID1,x$SPECID2),qual)) 
# Make the comparision

intra_freq<-as_tibble(intra_freq)
intra_freq<-dplyr::mutate(intra_freq,within_host_time=collect2-collect1)
just_one<-intra_freq %>% dplyr::group_by(SPECID1,SPECID2) %>%
  dplyr::summarize(within_host_time=unique(within_host_time)) # Just get one within_host_time/ pair

# require(ggplot2)
# intra_time.p<-ggplot(just_one,aes(x=within_host_time))+geom_histogram(binwidth =1, color="white")+xlab("Days between samples")+ ylab("Sample Pairs")+stat_bin(aes(y=..count.., label=..count..), binwidth = 1,geom="text", vjust=-.5)
# intra_time.p



# We are only looking at polymorphic sites between 2% and 98%.


intra_freq.comp<-polish_freq(intra_freq,freq1,0.02)
# Now we add the class of the isnv in the first sample ie the donor sample
intra_freq.comp<-plyr::adply(intra_freq.comp,1,function(x){
  donor_class = subset(qual,SPECID==x$SPECID1[1] & 
                         mutation == x$mutation[1],select=c(class_factor))
  #print(x)
  x$donor_class = as.character(donor_class$class_factor)
  return(x)
})

write.csv(x = intra_freq,file = "../../data/processed/secondary/Intrahost_all.csv")
write.csv(x = intra,file =  "../../data/processed/secondary/Intrahost_pairs.csv")
write.csv(x = intra_freq.comp,file = "../../data/processed/secondary/Intrahost_initially_present.csv")

