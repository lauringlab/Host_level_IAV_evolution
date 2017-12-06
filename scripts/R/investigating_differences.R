require(magrittr)
require(tidyverse)
require(HIVEr)
require(multidplyr)
require(testthat)
require(ggplot2)


if(!("package:tidyverse" %in% search())){
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(readr)
  require(rlang)
}
qual<-read_csv("../../data/processed/secondary/qual.snv.csv",
col_types = list(
ENROLLID= col_character(),
SPECID = col_character(),
LAURING_ID = col_character(),
Id = col_character()
)) # read in quality variant calls from all


old_qual<-read_csv("../../../qual.snv.csv",
col_types = list(
ENROLLID= col_character(),
SPECID = col_character(),
LAURING_ID = col_character(),
Id = col_character()
)) # read in quality variant calls from all
dim(old_qual)
dim(qual)
min_qual<-subset(qual,freq.var<0.5)
min_old_qual<-subset(old_qual,freq.var<0.5)
dim(min_old_qual)
dim(min_qual)

anti_join(min_old_qual,min_qual,by = c("SPECID","chr","pos","ref","var","ENROLLID"))->lost
dim(lost)

lost
lost %>% select(SPECID,chr,pos,ref,var,freq.var,Id)
# Lost should just be MH8156 variants 
# 
# variants %>% filter(SPECID=="MH8156",chr=="PB1",pos==1392|pos==1394) %>% 
#   select(SPECID,chr,pos,ref,var,freq.var,Id,run,cov.tst.fw) %>% 
#   arrange(pos,var,run)
#   
#   
# Are lost_major variants at positions that are not variable in the dataset?
lost_all <- anti_join(old_qual,qual,by = c("SPECID","chr","pos","ref","var","ENROLLID"))
lost_all
summary(lost_all$freq.var)


lost_all %>% ggplot(aes(x=freq.var))+
  geom_histogram(binwidth = 0.02,color='white')+scale_y_log10()

# Are lost_major variants at positions that are not variable in the dataset?

alleles <- old_qual %>% group_by(season,pcr_result,chr,pos)%>%
  summarize(alleles=length(unique(var)))

lost_all<-left_join(lost_all,alleles)

lost_all %>% ggplot(aes(x=alleles))+
  geom_histogram(binwidth = 0.5,color='white')+scale_y_log10()

summary(lost_all$alleles)

lost_all %>% filter(alleles==2) %>% 
  select(SPECID,ENROLLID,freq.var,chr,pos,ref,var)
# There is one allele that had variability.
# 

old_qual %>%filter(SPECID=="HS1448",chr=="HA",pos==1598)%>%
  select(SPECID,ENROLLID,freq.var,chr,pos,ref,var)

old_qual %>%filter(chr=="HA",pos==1598,season=="2014-2015")%>%
  select(SPECID,ENROLLID,freq.var,chr,pos,ref,var)

qual %>%filter(chr=="HA",pos==1598,season=="2014-2015")%>%
  select(SPECID,ENROLLID,freq.var,chr,pos,ref,var)

### Distances 
### 

# It is clear why the MH8156 is different. The others must be different because of 
# sites with one variant present (in the data set and individual) but a freq.var of less than 98%. 
# These would be removed in the new analysis and set to 100% in equal_compare but 
# could have remained in the old.

l1<-read_csv("../../data/processed/secondary/possible.pairs.dist.csv")
old_l1<-read_csv("../../results/possible.pairs.dist.csv")
dim(old_l1)
dim(l1)
anti_join(l1,old_l1,by=c("SPECID1","SPECID2","L1_norm")) %>% 
  select("SPECID1","SPECID2","L1_norm","HOUSE_ID1","HOUSE_ID2","ENROLLID1","ENROLLID2")
anti_join(old_l1,l1,by=c("SPECID1","SPECID2","L1_norm"))%>% 
  select("SPECID1","SPECID2","L1_norm","HOUSE_ID1","HOUSE_ID2","ENROLLID1","ENROLLID2")

# Intrahost - needs to be run 
intra<-read_csv("../../data/processed/secondary/Intrahost_initially_present.csv")
intra_old<-read_csv("../../results/Intrahost_initially_present.csv")

#There is only 1 difference here it is in both directions
anti_join(intra_old,intra,by=c("SPECID1","SPECID2","chr","pos","ref","var",
                               "freq1","freq2","gc_ul1","gc_ul2","donor_class","mutation"))

different<-anti_join(intra,intra_old,by=c("SPECID1","SPECID2","chr","pos","ref","var",
                               "freq1","freq2","gc_ul1","gc_ul2","donor_class","mutation"))

intra %>% filter(chr=="NR",pos==168) %>% select(chr,pos,SPECID1,SPECID1,ref,var,freq1,freq2,within_host_time)
intra_old %>% filter(chr=="NR",pos==168) %>% select(chr,pos,SPECID1,SPECID1,ref,var,freq1,freq2,within_host_time)

# Transmission
trans <- read_csv("../../data/processed/secondary/trans_freq.csv")
trans_old <- read_csv("../../results/trans_freq.csv")

trans<- trans %>% select(-X1,-X1_1)
trans_old <- trans_old %>% select(-X1)

dim(trans)
dim(trans_old)

anti_join(trans,trans_old,by=c("ENROLLID1","ENROLLID2",
                               "SPECID1","SPECID2",
                               "chr","pos","ref","var",
                               "freq1","freq2","HOUSE_ID",
                               "onset1","onset2","collect1","collect2")) 

### Trans freq comp

trans.comp <- read_csv("../../data/processed/secondary/transmission_pairs_freq.poly.donor.csv")
trans_old.comp <- read_csv("../../results/transmission_pairs_freq.poly.donor.csv")

trans.comp<- trans.comp %>% select(-X1,-X1_1)
trans_old.comp <- trans_old.comp %>% select(-X1)

dim(trans.comp)
dim(trans_old.comp)

anti_join(trans.comp,trans_old.comp,by=c("ENROLLID1","ENROLLID2",
                               "SPECID1","SPECID2",
                               "chr","pos","ref","var",
                               "freq1","freq2","HOUSE_ID",
                               "onset1","onset2","collect1","collect2")) 

################### alleles and all that jazz #######################

# How many loci have 1 allele but questionable frequencies.
# 

add_alleles<-function(df,...){
  group_var <- rlang::quos(...)
  df_alleles <- df %>% dplyr::group_by(!!!group_var) %>% 
    dplyr::summarize(alleles = length(unique(var)))
  df <- dplyr::left_join(df, df_alleles)
}

qual<-add_alleles(qual,SPECID,chr,pos,pcr_result)

one_allele<-subset(qual,alleles==1 & freq.var<0.98)

ggplot(one_allele,aes(x=freq.var))+geom_histogram(binwidth = 0.02,color="white")+
  scale_y_log10()


old_qual<-add_alleles(old_qual,SPECID,chr,pos,pcr_result)

old_one_allele<-subset(old_qual,alleles==1 & freq.var<0.9)

ggplot(old_one_allele,aes(x=freq.var))+geom_histogram(binwidth = 0.02,color="white")+
  scale_y_log10()



old_one_allele %>% filter(SPECID %in% trans_old$SPECID1| SPECID %in% trans_old$SPECID2) %>%
  ggplot(aes(x=freq.var))+geom_histogram(binwidth = 0.02,color="white")+ scale_y_log10()

old_one_allele %>% filter(SPECID %in% trans_old$SPECID1| SPECID %in% trans_old$SPECID2) %>%
  select(SPECID) %>%unique()%>% nrow()

# See if these are present or not in the deepSNV output
subset(old_one_allele,freq.var<0.6) %>% head() %>% select(Id,chr,pos,freq.var,ref,var)


vic <-read_csv("../../data/processed/victoria/deepSNV/1230.removed.csv")
vic %>% filter(chr=="PB2", pos==1921)
vic_2<-read_csv("../../data/processed/victoria_2/deepSNV/1230.removed.csv")
vic_2 %>% filter(chr=="PB2", pos==1921)

# I need the bamfile in igv for this Yes deletions play a role
# 


one_allele_distributions<- old_one_allele %>% group_by(season,pcr_result,chr,pos,var,ref) %>% 
  summarize(samples=length(unique(SPECID)))

one_allele_distributions$chr<-factor(one_allele_distributions$chr,
                           levels = rev(c("PB2","PB1","PA","HA","NP","NR","M","NS"))) 
# Set the segments as factors with PB2 on top
chrs<-read.csv("../../data/reference/segs.csv",
               stringsAsFactors = T) 
# get the start and stop of each OR for each segment (2014-2015 used as reference)

chrs$chr<-factor(chrs$chr,levels=levels(one_allele_distributions$chr)) # set factors on the is meta data
genome_loc.p<-ggplot(one_allele_distributions,aes(x=pos,y=chr))+
  geom_point(aes(color=samples),shape=108,size=5)+
  geom_segment(data=chrs,aes(x = start, y = chr, xend = stop, yend = chr))+
  ylab("")+
  xlab("")+
  theme(axis.ticks =element_blank(),
        axis.line.x = element_blank(),
        axis.line.y=element_blank())+
  scale_x_continuous(breaks=c())
genome_loc.p+facet_wrap(~season)

one_allele_distributions %>% ggplot(aes(x=samples))+
  geom_histogram(color="white",binwidth = 1)+scale_y_log10()
