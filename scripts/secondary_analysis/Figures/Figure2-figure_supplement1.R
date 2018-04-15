require(tidyverse)
require(magrittr)
require(extrafont)
require(HIVEr)
require(cowplot)

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
              # this data file is made in processing_snv.R it contains allele calls at postions
              # that had multiple alleles present in the data set (not necessisarily SNV it could
              # be fixed differences between samples). It is written to file right before the call
              # to the qual function that require duplicate calls for samples in lower titer 
              # samples.
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

qual<-read_csv("./data/processed/secondary/qual.snv.csv", 
               # These are the iSNV that were used in the study.
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all
# Intrahost data 
intra<-read_csv("./data/processed/secondary/Intrahost_all.csv")
# quality SNV that passed qualtiy thresholds but do not necessarily pass the frequency cut off
# It was not applied to these iSNV.
no_freq_cut<-read_csv("./data/processed/secondary/no_freq_cut.qual.snv.csv",
                      col_types = list(
                        ENROLLID= col_character(),
                        SPECID = col_character(),
                        LAURING_ID = col_character(),
                        Id = col_character()
                      ))

# ----------------------------- data processing ------------------------------
#  All duplicates were sequenced on separate runs so we will use run as the 
#  duplicate label. Id may also be an option. Again var refers to iSNV that
#  pass quality thresholds (excluding frequency) but did not have to be found 
#  in duplicate runs.
# ----------------------------------------------------------------------------

# # Just get samples that were sequenced in duplicate
# sequenced<-var %>%
#   group_by(SPECID) %>%
#   summarize(replicates = length(unique(run)))
# 
# var<-left_join(var,sequenced) %>% filter(replicates>1)
# 
# # This is to ensure the SPECID now refers to the sample and run
# var <- var %>%
#   rename("SPECID_original" = "SPECID") %>% # save the original
#   mutate(SPECID=paste(run,SPECID_original,sep="_"))
# # set up
# # get meta data on samples seqeunced in duplicate
# replicates<-var %>% select(SPECID,SPECID_original,home_collected,collect,
#                   gc_ul,season,pcr_result)%>% distinct()
# 
# dup_col<-c("SPECID")
# 
# # Get the SPECID_run column in SPECID1 and SPeCID2 columns for use by get_freqs below
# This will compare the frequencies of mutations found in both sequencing runs.
# pairs<- short_pairs(replicates,dup_col,SPECID_original) 
# 
# freqs<-pairs %>%
#   rowwise() %>%
#   do(get_freqs(c(.$SPECID1,.$SPECID2),var))
#   # add the meta data back
# freqs<-left_join(freqs,pairs)
# 
# # write this to limit run time in the future
# write.csv(freqs,"./data/processed/secondary/duplicate_sequences.csv")
# 

# ----------------------- Figure1A  -----------------------------
# Frequency measurements in biological replicates. ie 2 samples
# from the same person on the same day.
# ---------------------------------------------------------------
intra<- intra %>% filter(!(SPECID2 %in% c("HS1530","MH8137","MH8390")) &
                           !(SPECID1 %in% c("HS1530","MH8137","MH8390")))
intra %>% mutate(DPS1 = collect1-onset,DPS2 = collect2-onset) ->intra

intra<-filter(intra,freq1<0.5) 

same_day<-subset(intra,within_host_time==0 & freq1<0.5)
Sup_6<-ggplot(same_day,aes(x=freq1,y=freq2))+geom_point()+
  xlab("Frequency in home isolate") + ylab("Frequency in clinic isolate") + 
  geom_abline(slope=1,intercept = 0,lty=2)+
  scale_x_continuous(limits = c(0,0.5))+scale_y_continuous(limits = c(0,0.5))
Sup_6

lm_fit<-lm(freq2~freq1,same_day)
summary(lm_fit)->sum_fit
sum_fit

write_to_summary("R2 samples same day:",sum_fit$r.squared)
save_plot("./results/Figures/Figure2-figure_supplement1A.pdf", Sup_6,
          base_aspect_ratio = 1)
embed_fonts("./results/Figures/Figure2-figure_supplement1A.pdf")



# -------------------- intra processing --------------------------
# - taken from Figure2.R 1decf3dc6df24c6e2b4d58cd3fdea24652fd46d1
# Reidenfiy the frequency of lost and arisen mutation using the dataframe
# that does not include the frequency cut off.
# ----------------------------------------------------------------
# Remove mixed infections
intra<- intra %>% filter(!(SPECID2 %in% c("HS1530","MH8137","MH8390")) &
                           !(SPECID1 %in% c("HS1530","MH8137","MH8390")))
intra %>% mutate(DPS1 = collect1-onset,DPS2 = collect2-onset) ->intra

intra<-filter(intra,freq1<0.5) 

arisen<-intra %>% filter(freq1==0) %>% 
  rowwise() %>% 
  mutate(freq_close_look = 
           ifelse(length(which(no_freq_cut$SPECID==SPECID1 & no_freq_cut$mutation==mutation))==1,
                  no_freq_cut$freq.var[which(no_freq_cut$SPECID==SPECID1 & no_freq_cut$mutation==mutation)],
                  0),
         freq1=freq_close_look) %>%
  select(-freq_close_look)

lost<-intra %>% filter(freq2==0) %>% 
  rowwise() %>% 
  mutate(freq_close_look = 
           ifelse(length(which(no_freq_cut$SPECID==SPECID2 & no_freq_cut$mutation==mutation))==1,
                  no_freq_cut$freq.var[which(no_freq_cut$SPECID==SPECID2 & no_freq_cut$mutation==mutation)],
                  0),
         freq2=freq_close_look) %>%
  select(-freq_close_look)

intra_processed<-rbind(arisen,lost,filter(intra,freq1>0,freq2>0))

intra_processed<-intra_processed %>% mutate(Endpoint="Persistent") %>%
  mutate(Endpoint = if_else(freq1==0,"Arisen",Endpoint)) %>%
  mutate(Endpoint = if_else(freq2==0,"Lost",Endpoint))

intra_processed$Endpoint<-factor(intra_processed$Endpoint,
                                 levels = c("Persistent","Arisen","Lost"),ordered = T)
# this was written as the data file for Figure 2e in 1decf3dc6df24c6e2b4d58cd3fdea24652fd46d1

# ---------------------------------- figure ----------------------------------
#  Read in the csv from above 
#  ToDo
#  Read in qual - only use mutations in qual
#  use frequency reported in qual to report relative difference
# ----------------------------------------------------------------------------



freqs<-read_csv("./data/processed/secondary/duplicate_sequences.csv")

# now we will select only iSNV used in the analysis
# To do this I will make a SPECID_mutation column linking each iSNV with the isolate
# in which it was found.
# 
quality_mut<-paste0(qual$mutation,qual$SPECID)

# Freqs are the mutations found without frequency cutoff. In freq1 freq2 format
freqs<-mutate(freqs, mut_specid = paste0(mutation,SPECID_original), # here we use SPECID_original to match that used in qual above
              used= mut_specid %in% quality_mut)

# Those not included are mutations relative to the plasmid control that are fixed.
# qual<- qual %>% mutate(mut_specid = paste0(qual$mutation,qual$SPECID))
# not_included<-qual %>% filter(gc_ul<1e5,!(mut_specid %in% freqs$mut_specid))
# not_included %>% select(freq.var) %>% unique()
# # This plot is just for exploration
# freqs %>% filter(ref!=var) %>% 
#   ggplot(aes(x=freq1,y=freq2))+
#     geom_point(aes(color=used))+geom_abline(slope=1,intercept = 0,lty=2)+
#     scale_y_log10()+
#     scale_x_log10()+
#     #geom_smooth(data=filter(freqs,freq1>0,freq2>0),method=lm,se=T)+
#     xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 2")+
#     scale_color_manual(values=c("gray","black"),name="",labels=c("Removed","Used"))

# # This one uses a continous color gradiant
# dot_plot<-freqs%>% filter(used==TRUE)%>% # filter so only those used are plotted
# ggplot(aes(x=freq1,y=freq2))+
#   geom_point(aes(color=log10(gc_ul)))+
#   geom_abline(slope=1,intercept = 0,lty=2)+
#   scale_y_continuous(limits=c(0,0.5))+
#   scale_x_continuous(limits=c(0,0.5))+
#   xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 2")+
#   scale_color_continuous(name="Log(copies/ul)")
#   
# save_plot("./results/Figures/sequencing_dot_plot.pdf", dot_plot,
#           base_aspect_ratio = 1.3)
# embed_fonts("./results/Figures/sequencing_dot_plot.pdf")
# 
# write.csv(select(filter(freqs,used==T),gc_ul,freq1,freq2),
#           "./results/Figures/data/sequencing_dot_plot.csv")

dot_plot.discrete<-freqs%>% filter(used==TRUE)%>% # filter so only those used are plotted
  ggplot(aes(x=freq1,y=freq2))+
  geom_point(aes(color=as.factor(floor(log10(gc_ul)))))+
  geom_abline(slope=1,intercept = 0,lty=2)+
  scale_y_continuous(limits=c(0,0.5))+
  scale_x_continuous(limits=c(0,0.5))+
  xlab("Frequency in replicate 1")+ ylab("Frequency in replicate 2")+
  scale_color_manual(name="Log(copies/ul)",values=cbPalette[c(4,1)])

dot_plot.discrete
save_plot("./results/Figures/Figure2-figure_supplement1B.pdf", dot_plot.discrete,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2-figure_supplement1B.pdf")

intra<-intra_processed
intra <- intra %>% 
  mutate(rel_difference = if_else(freq1!=0,abs(freq1-freq2)/freq1,
                              abs(freq1-freq2)/freq2),
         difference = abs(freq1-freq2))
freqs<- freqs %>%
  mutate(rel_difference = abs(freq1-freq2)/freq1,
         difference = abs(freq1-freq2))

difference <- rbind(
  tibble(difference = 
           freqs$difference[which(freqs$used==T & freqs$freq1<0.5)],
         group="Measure"),
  tibble(difference = intra$difference,group="Intrahost"))
# require(ggridges)
# ggplot(difference,aes(x=difference,y=group))+geom_density_ridges2(scale=6)+
#   xlab("Relative difference")+
#   ylab("")+scale_y_discrete(labels=c("Intrahost dynamics","Measurment Error"))


difference_histogram<-ggplot(difference,aes(x=difference,fill=group))+
  geom_histogram(aes(y= ..ncount..), position="dodge")+
  xlab("Frequency difference")+
  ylab("Normalized count")+
  scale_x_log10(breaks = c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1))+
  scale_fill_manual(values = cbPalette[c(1,5)],
                    labels = c("Intrahost dynamics","Measurment Error"),
                    name="")+
  theme(legend.position = c(0.1,0.5))

save_plot("./results/Figures/Figure2-figure_supplement1C.pdf", difference_histogram,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2-figure_supplement1C.pdf")

# write.csv(difference,
#           "./results/Figures/data/sequencing_hist.csv")




