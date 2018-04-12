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
plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}

# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
# A
meta<-read_csv("./data/reference/all_meta.sequence_success.csv")
qual<-read_csv("./data/processed/secondary/qual.snv.csv", # The quality iSNV
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all 
# get the start and stop of each OR for each segment (2014-2015 used as reference)
chrs<-read.csv("./data/reference/segs.csv",stringsAsFactors = T) 
# --------------------------------- Figure 1A ---------------------------------
#   The distribution of titers over time
# -----------------------------------------------------------------------------

titers<-subset(meta,!(is.na(DPI))) #)& DPI>=0 ) 
# We only want the titers where we measured a titer 
# ( I beleive there are 2 samples here with NA). 
# Also 2 samples were collected before symptoms - according to the meta data. 
figure1A_data<-select(titers,SPECID,
                      days.post.symptom.onset=DPI,
                      genomes.per.ul = gc_ul)

titer.p<-ggplot(figure1A_data,aes(x=as.factor(days.post.symptom.onset),
                                  y=genomes.per.ul))+
  geom_boxplot(notch = T)+
  ylab(expression(paste(Genomes,"/" ,mu,L)))+
  scale_y_log10()+xlab("Days post symptom onset")
save_plot("./results/Figures/Figure1A.pdf", titer.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1A.pdf")

write.csv(figure1A_data,"./results/Figures/data/Figure_1-source_data_1.csv")

# --------------------------------- Figure 1B ---------------------------------
#   The number of iSNV in each isolate stratified by DPI
# -----------------------------------------------------------------------------
min.qual.o<-subset(qual,freq.var<0.5) # only looking at minor alleles original - includes all samples
min.qual<-subset(min.qual.o,!(SPECID %in% c("HS1530","MH8137","MH8390"))) # remove mixed infections
# How many rare mutations in the sample (SPECID)
min.count.sample.o<-min.qual.o %>% group_by(SPECID) %>%
  dplyr::summarize(iSNV=length(unique(mutation)),HA_iSNV=length(which(chr=="HA"))) 


snv_qual_meta.o<-subset(meta,snv_qualified==T)
snv_qual_meta.o<-merge(snv_qual_meta.o,min.count.sample.o,by="SPECID",all.x=T)

snv_qual_meta.o$iSNV[is.na(snv_qual_meta.o$iSNV)]<-0 # these are the ones with no diversity
snv_qual_meta.o$HA_iSNV[is.na(snv_qual_meta.o$HA_iSNV)]<-0

figure1B_data<-snv_qual_meta.o %>% select(days.post.symptom.onset=DPI,iSNV)
isnv_by_day.p<-ggplot(snv_qual_meta.o,aes(x=as.factor(DPI),y=iSNV))+
  #geom_quasirandom()+
  #stat_summary(fun.data="plot.median", geom="errorbar", colour="red", width=0.55, size=0.5)
  geom_boxplot()+xlab("Day post symptom onset")
write.csv(figure1B_data,"./results/Figures/data/Figure_1-source_data_2.csv")

save_plot("./results/Figures/Figure1B.pdf", isnv_by_day.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1B.pdf")

# --------------------------------- Figure 1C ---------------------------------
#   The number of iSNV in each isolate stratified by vaccination
# -----------------------------------------------------------------------------

isnv_by_vaccination<-ggplot(snv_qual_meta.o,aes(y=iSNV,x=as.factor(vaccination_status)))+
#geom_quasirandom()+
  geom_dotplot(stackdir = "center",binaxis = 'y',binwidth = 1,dotsize = 0.5)+
  stat_summary(fun.data="plot.median", geom="errorbar", colour="red", width=0.95, size=0.5)+
  scale_x_discrete(labels=c("Not Vaccinated","Vaccinated"))+xlab("")
figure1C_data<-snv_qual_meta.o %>% select(vaccination_status,iSNV)

save_plot("./results/Figures/Figure1C.pdf", isnv_by_vaccination,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1C.pdf")

write.csv(figure1C_data,"./results/Figures/data/Figure_1-source_data_3.csv")


# --------------------------------- Figure 1D ---------------------------------
#   iSNV position across the genome - highlighting mutations found in  more than
#   one individual
# -----------------------------------------------------------------------------
# Set the segments as factors with PB2 on top
min.qual$chr<-factor(min.qual$chr,levels = rev(c("PB2","PB1","PA","HA","NP","NR","M","NS")))

chrs$chr<-factor(chrs$chr,levels=levels(min.qual$chr)) # set factors on the is meta data
# Plot on genome
genome_loc.p<-ggplot(min.qual,aes(x=pos,y=chr))+
  geom_point(aes(color=class_factor),shape=108,size=5)+
  geom_segment(data=chrs,aes(x = start, y = chr, xend = stop, yend = chr))+
  ylab("")+
  xlab("")+
  scale_color_manual(name="",values=cbPalette[c(1,4)])+
  theme(axis.ticks =element_blank(),
        axis.line.x = element_blank(), axis.line.y=element_blank())+
  scale_x_continuous(breaks=c())+
  theme(legend.position = "none")

# add multiple iSNV
# In how many people (ENROLLID) is the muation found. It would have to be minor in both 
# people to count in this plot. The chr pos var is used so that called and infered variants 
# of the same var are counted together.

min.count<- min.qual %>% group_by(chr,pos,var,pcr_result,season) %>%
  dplyr::summarize(counts=length(unique(ENROLLID))) 
min.qual<-mutate(min.qual,scratch = paste(chr,pos,chr,var,season,pcr_result,sep="."))
min.count <-mutate(min.count,scratch = paste(chr,pos,chr,var,season,pcr_result,sep="."))
multiple<-subset(min.qual,scratch %in% min.count$scratch[min.count$counts>1])

min.qual<-subset(min.qual,select=-c(scratch))
min.count<-subset(min.count,select=-c(scratch))

genome_loc.p.dots<-genome_loc.p + 
  geom_point(data = multiple, aes(x=pos,y=as.numeric(chr)+0.3,color=class_factor),size=1,shape=6)

save_plot("./results/Figures/Figure1D.pdf", genome_loc.p.dots,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1D.pdf")
figure1D_data<-min.qual %>% select(chr,pos,var,freq.var,class=class_factor,ENROLLID,season,strain = pcr_result)
write.csv(figure1D_data,"./results/Figures/data/Figure_1-source_data_4.csv")
