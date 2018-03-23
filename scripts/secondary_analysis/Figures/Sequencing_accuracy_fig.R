require(tidyverse)
require(magrittr)
require(cowplot)
require(extrafont)
loadfonts()

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
# -------------------------------- Comparison --------------------------------
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
  
write.csv(freqs,"./data/processed/secondary/duplicate_sequences.csv")

ggplot(freqs,aes(x=freq1,y=freq2))+
  geom_point()+geom_abline(slope=1,intercept = 0,lty=2)+
  scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  geom_smooth(data=filter(freqs,freq1>0,freq2>0),method=lm)+
  xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 1")

# used 
# # this is not correct as writted. 
# We should calculate the coverage as we do in the analysis.
# 
used <- freqs %>%
  filter(max(freq1,freq2)>0.02) 
ggplot(used,aes(x=freq1,y=freq2))+
  geom_point()+geom_abline(slope=1,intercept = 0,lty=2)+
  #scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  #scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1))+
  geom_smooth(data=filter(used,freq1>0,freq2>0),method=lm)+
  xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 1")

freqs<- freqs %>%
  mutate(difference =  abs(freq1-freq2))

freqs %>% filter(freq1>0,freq2>0,min(freq1,freq2)<0.5)->x

x<-x  %>% mutate(percent_diff = difference/mean(freq1,freq2))

ggplot(x,aes(x=difference))+ geom_histogram(color="white")+scale_x_log10()

ggplot(filter(x,max(freq1,freq2)>0.02),aes(x=difference))+ 
  geom_histogram(color="white")+scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1))


ggplot(filter(x,max(freq1,freq2)>0.02),aes(x=percent_diff))+ 
  geom_histogram(color="white",bins=100)+scale_x_log10()+xlab("relative difference")


