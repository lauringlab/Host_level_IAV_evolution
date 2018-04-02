require(tidyverse)
require(magrittr)
require(cowplot)
require(extrafont)
require(HIVEr)
cbPalette<-wesanderson::wes_palette("Zissou") # Set up figures

# --------------------------------- Functions ---------------------------------
write_to_summary<-function(line_pattern,value){
  file = readLines("./results/results.table.tsv")
  line_pattern_regex = paste0("^",line_pattern)
  line = grep(line_pattern_regex,file)
  file[line] = paste0(line_pattern,"\t",value)
  writeLines(file,"./results/results.table.tsv")
}


# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
# A
var<-read_csv("./data/processed/secondary/diverse_sites_isnv.csv", # The  iSNV
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all
# 
# get rid of variants that are only noncoding
var<- var%>% filter(coding_pos !="['Noncoding']")
chrs<-read.csv("./data/reference/segs.csv",stringsAsFactors = T) 

qual<-read_csv("./data/processed/secondary/qual.snv.csv", # The  iSNV
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all

intra<-read_csv("./data/processed/secondary/Intrahost_all.csv")

# ----------------------------- data processing ------------------------------
#  All duplicates were sequenced on separate runs so we will use run as the 
#  duplicate label. Id may also be an option.
# ----------------------------------------------------------------------------



#var <- filter(var,season=="2010-2011",
#                    pcr_result=="A/H3N2")

sequenced<-var %>%
  group_by(SPECID) %>%
  summarize(replicates = length(unique(run)))

var<-left_join(var,sequenced) %>% filter(replicates>1)

# This is to ensure the SPECID now refers to the sample and run
var <- var %>%
  rename("SPECID_original" = "SPECID") %>%
  mutate(SPECID=paste(run,SPECID_original,sep="_"))
# set up

replicates<-var %>% select(SPECID,SPECID_original,home_collected,collect,
                  gc_ul,season,pcr_result)%>% distinct()

dup_col<-c("SPECID")

pairs<- short_pairs(replicates,dup_col,SPECID_original)


freqs<-pairs %>%
  rowwise() %>%
  do(get_freqs(c(.$SPECID1,.$SPECID2),var))
freqs<-left_join(freqs,pairs)

write.csv(freqs,"./data/processed/secondary/duplicate_sequences.csv")


# ---------------------------------- figure ----------------------------------
#  Read in the csv from above 
#  ToDo
#  Read in qual - only use mutations in qual
#  use frequency reported in qual to report relative difference
# ----------------------------------------------------------------------------



freqs<-read_csv("./data/processed/secondary/duplicate_sequences.csv")
ggplot(freqs,aes(x=freq1,y=freq2))+
  geom_point()+geom_abline(slope=1,intercept = 0,lty=2)+
  #scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  #scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  geom_smooth(data=filter(freqs,freq1>0,freq2>0),method=lm,se=T)+
  xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 2")

# What's with the odd ones on the axis at >50% in one replicate?
odd<-freqs %>% filter(freq1>0.5, freq2==0)

odd$chr<-factor(odd$chr,levels = rev(c("PB2","PB1","PA","HA","NP","NR","M","NS")))

chrs$chr<-factor(chrs$chr,levels=levels(odd$chr)) # set factors on the is meta data
# Plot on genome
genome_loc.p<-ggplot(odd,aes(x=pos,y=chr))+
  geom_point(shape=108,size=5)+
  geom_segment(data=chrs,aes(x = start, y = chr, xend = stop, yend = chr))+
  ylab("")+
  xlab("")+
  scale_color_manual(name="",values=cbPalette[c(1,4)])+
  theme(axis.ticks =element_blank(),
        axis.line.x = element_blank(), axis.line.y=element_blank())+
  scale_x_continuous(breaks=c())+
  theme(legend.position = "none")
genome_loc.p

# They are just on the ends of segements (except for a few in NS. They can be
# removed by only looking in the OR).
# A few remain. They have very low coverage in the test sample.
# freqs %>% filter(freq1>0.8, freq2==0)
# Freq 1 =0
# var %>% filter(chr=="NS",freq.var>0.95, SPECID=="perth_2_MH1346",
#  pos %in% c(90,168,263,699)) %>% select(cov.tst.fw,cov.tst.bw,sample_coverage)
#  4 of these are the inverse of those variants that were only found int replicate1 
#  in MH1346
#  The others match cases where there are biases in one sample. These are removed
#  from the others.
#  
#  


quality_mut<-paste0(qual$mutation,qual$SPECID)

freqs<-mutate(freqs, mut_specid = paste0(mutation,SPECID_original),
              used= mut_specid %in% quality_mut,
              difference = abs(freq1-freq2)/freq1)


ggplot(freqs,aes(x=freq1,y=freq2))+
  geom_point(aes(color=used))+geom_abline(slope=1,intercept = 0,lty=2)+
  scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  geom_smooth(data=filter(freqs,freq1>0,freq2>0),method=lm,se=T)+
  xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 2")+
  scale_color_manual(values=c("gray","black"),name="",labels=c("Removed","Used"))


# Just to make sure
intra<-filter(intra,freq1<0.5) 
# Remove duplicate infections
intra<- intra %>% filter(!(SPECID2 %in% c("HS1530","MH8137","MH8390")) &
                           !(SPECID1 %in% c("HS1530","MH8137","MH8390")))
intra %>% mutate(DPS1 = collect1-onset,DPS2 = collect2-onset) ->intra
intra<-intra %>% mutate(Endpoint="Persistent") %>%
  mutate(Endpoint = if_else(freq1==0,"Arisen",Endpoint)) %>%
  mutate(Endpoint = if_else(freq2==0,"Lost",Endpoint))
# no arisen
intra <- intra %>% filter(Endpoint !="Arisen") %>% 
  mutate(difference = abs(freq1-freq2)/freq1)

difference <- rbind(
  tibble(difference = 
           freqs$difference[freqs$used==T],
         group="Measure"),
  tibble(difference = intra$difference,group="Intrahost"))
require(ggridges)
ggplot(difference,aes(x=difference,y=group))+geom_density_ridges2(scale=6)+
  xlab("Relative difference")+
  ylab("")+scale_y_discrete(labels=c("Intrahost dynamics","Measurment Error"))


ggplot(difference,aes(x=difference,fill=group))+
  geom_histogram(aes(y=..ncount..), position="dodge")+
  xlab("Relative Frequency difference")+
  ylab("Normalized count")+
  scale_x_log10()+
  scale_fill_manual(values = cbPalette[c(1,2)],
                    labels = c("Intrahost dynamics","Measurment Error"),
                    name="")+
  theme(legend.position = c(0.2,0.5))
