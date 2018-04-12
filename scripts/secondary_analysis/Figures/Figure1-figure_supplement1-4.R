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
cov.files<-c("./data/processed/HK_1/all.coverage.csv","./data/processed/HK_2/all.coverage.csv","./data/processed/HK_6/all.coverage.csv","./data/processed/HK_7/all.coverage.csv","./data/processed/HK_8/all.coverage.csv","./data/processed/cali09/all.coverage.csv","./data/processed/cali09_2/all.coverage.csv","./data/processed/victoria/all.coverage.csv","./data/processed/victoria_2/all.coverage.csv","./data/processed/perth/all.coverage.csv","./data/processed/perth_2/all.coverage.csv")
cov<-read_rbind(cov.files,
                3,
                cols = cols(
                  chr = col_character(),
                  chr.pos = col_integer(),
                  coverage = col_integer(),
                  concat.pos = col_integer(),
                  Id = col_character()
                ))

meta<-read_csv("./data/reference/all_meta.sequence_success.csv")
qual<-read_csv("./data/processed/secondary/qual.snv.csv", # The quality iSNV
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all 

# --------------------------------- Functions --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------
# coverage
slide<-function(cov.df,setup.df){
  
  coverage = rep(NA,nrow(setup.df))
  for( i in 1:nrow(setup.df)){
    s = setup.df$starts[i]
    e = setup.df$ends[i]
    subset(cov.df,concat.pos>=s  & concat.pos<e,select=c(coverage)) -> position
    mean(position$coverage)->coverage[i]
  }
  out<-data.frame(mean = coverage,concat.pos = setup.df$concat.pos,chr=setup.df$chr)
  out$Id = unique(cov.df$Id)
  out$run = unique(cov.df$run)
  return(out)
}
cov_plot<-function(cov.df,title){
  
  ## Get the steps and positions for each chr
  
  cov.df %>% group_by(chr) %>% 
    summarize(first  = min(concat.pos),last = max(concat.pos)) %>% 
    plyr::adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=100))) %>%
    mutate(ends = ifelse(starts+200<last,starts+200,last)) %>%
    rowwise() %>%
    mutate(concat.pos = mean(c(starts,ends)))  -> setup
  plyr::ddply(cov.df,~Id+run,slide,setup) -> cov.slid.df
  
  # cov.slid.df<-ddply(cov.df,~Id+chr,function(x) slide(x,setup)
  x.labels<-plyr::ddply(cov.slid.df,~chr,plyr::summarize,concat.pos=concat.pos[which(abs(concat.pos-mean(concat.pos))==(min(abs(concat.pos-mean(concat.pos)))))]) # No idea if this will get what I want but heres to hoping!
  
  #of course sometimes there are 2 good choices I'll take the first one
  x.labels<-plyr::ddply(x.labels,~chr,function(x) return(x[1,]))
  
  #x.labels$chr[x.labels$chr %in% c("NR","N_A")]<-"NA"
  
  cov.plot<-ggplot(cov.slid.df, #subset(cov.slid.df,!(Sample%in%c("90","91","93"))),
                   mapping=aes(x=as.factor(concat.pos),
                               y=mean))+geom_boxplot(fill="white")
  
  cov.plot<-cov.plot+ggtitle(title)+ylab("Read depth")+scale_x_discrete(labels = x.labels$chr,breaks=x.labels$concat.pos)+xlab("Concatenated Genome Position")
  cov.plot<-cov.plot+theme(axis.title.y = element_text(vjust=1.2))
  cov.plot<-cov.plot+theme(legend.position="none")  
  return(cov.plot)
}

# --------------------------------- Figure 1 --------------------------------
#  Coverage
# ---------------------------------------------------------------------------


cov_plot(cov,title = "Coverage") ->coverage.plot


save_plot("./results/Figures/Figure1-figure_supplement1.pdf", coverage.plot,base_width = 15.0,
          base_aspect_ratio = 1.1)
embed_fonts("./results/Figures/Figure1-figure_supplement1.pdf")

# --------------------------------- Figure 2-3 --------------------------------
#  These are made by results/Investigating\ mixed\ infections.ipynb and saved as
#  tree files in results/coding_alignments/HXNX_coding.tree.
#  Various annoted versions also exist there.
#  sorry for the spaces in the file name!
# -----------------------------------------------------------------------------

# --------------------------------- Figure 4 --------------------------------
# The relationship between titer and iSNV richness
# ---------------------------------------------------------------------------

min.qual.o<-subset(qual,freq.var<0.5) # only looking at minor alleles original - includes all samples
min.count.sample.o<-min.qual.o %>% group_by(SPECID) %>%
  dplyr::summarize(iSNV=length(unique(mutation)),HA_iSNV=length(which(chr=="HA"))) # How many rare mutations in the sample (SPECID)
snv_qual_meta.o<-subset(meta,snv_qualified==T)
snv_qual_meta.o<-merge(snv_qual_meta.o,min.count.sample.o,by="SPECID",all.x=T)

snv_qual_meta.o$iSNV[is.na(snv_qual_meta.o$iSNV)]<-0 # these are the ones with no diversity
snv_qual_meta.o$HA_iSNV[is.na(snv_qual_meta.o$HA_iSNV)]<-0
isnv_titer<-ggplot(snv_qual_meta.o,aes(x=gc_ul,y=iSNV))+
  geom_point()+scale_x_log10()+
  xlab(expression(paste(Genomes,"/" ,mu,L)))


save_plot("./results/Figures/Figure1-figure_supplement4.pdf", isnv_titer,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure1-figure_supplement4.pdf")
 