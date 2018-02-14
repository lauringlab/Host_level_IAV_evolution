require(tidyverse)
require(extrafont)
loadfonts()
require(cowplot)
require(HIVEr)
require(ggbeeswarm)
require(directlabels)
require(lubridate)
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
qual<-read_csv("./data/processed/secondary/qual.snv.csv", # The quality iSNV
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all 
stopifnot(min(qual$freq.var)>0.02)
# B
antigenic<-read_csv("./data/processed/secondary/antigenic_isnv.csv")
# C
# This is for just antigenic mutations.
nextflu<-read.table("./data/processed/secondary/global_freq_antigenic.tsv"
                    ,header = T)

minor_HA.ns<-read.csv("./data/processed/secondary/minor_nonsynom.csv",
                stringsAsFactors = F,
                colClasses = c('onset'='Date','collect'='Date'))
#D
meta<-read_csv("./data/reference/all_meta.sequence_success.csv")
intra<-read_csv("./data/processed/secondary/Intrahost_all.csv")
# --------------------------------- Data processing ---------------------------
#   Any data processing beyond the usual data configuration goes here.
#   In this case this is the analysis used to get a list of HA1 minority 
#   variants that are found in our data set. This list was entered into nextflu
#   to obtain the global frequencies of the variants. Nextflu was accesed at 
#   http://nextflu.org/h3n2/ha/12y/ and data was download on February 2,2018
#    10:57
# -----------------------------------------------------------------------------

minor_HA.ns<-subset(minor_HA.ns,!is.na(Antigenic))
minor_HA.ns %>% subset(pcr_result=="A/H3N2") %>%.$H3_name ->H3
H3<-H3[H3!="-"]

paste(H3,sep=",",collapse = ", ")

minor_HA.ns %>% subset(pcr_result=="A/H1N1") %>%.$H1_name ->H1


# --------------------------------- Figure 2A ---------------------------------
#   The frequency distribution of NS and S data
# -----------------------------------------------------------------------------
# only looking at minor alleles original - includes all samples

min.qual.o<-subset(qual,freq.var<0.5) 
# remove mixed infections
min.qual<-subset(min.qual.o,!(SPECID %in% c("HS1530","MH8137","MH8390")))
# class factor is is Nonsym if in any open reading frame the mutation is Nonsym.
classes<-min.qual %>% group_by(class_factor) %>%
  dplyr::summarize(mutations=length(mutation)) 
#print(knitr::kable(classes))
write_to_summary("Ns/S:",
                 classes$mutations[classes$class_factor=="Nonsynonymous"]/
                   classes$mutations[classes$class_factor=="Synonymous"])

Figure2A_data<-min.qual %>% select(SPECID,chr,pos,freq.var,class_factor)

freq_hist.p<-ggplot(Figure2A_data,aes(x=freq.var,fill=class_factor))+
  geom_histogram(color="white",binwidth=.05,position=position_dodge(),
                 boundary = 0.02)+
  xlab("Frequency")+ylab("iSNV")+
  scale_fill_manual(name="",values=cbPalette[c(1,4)])+
  theme(legend.position = c(0.5, 0.5))

freq_hist.p

write.csv(Figure2A_data,"./results/Figures/data/Figure2A.csv") 

# Getting data from the plot
freq_hist.d<-ggplot_build(freq_hist.p)
freq_hist.d$data[[1]]->bins
bins$bin<-rep(1:10,each=2)

bin_ratio<- bins %>% group_by(bin) %>%
  dplyr::summarize(ratio = count[2]/count[1]) # NS/S
write_to_summary("Max Ns/S:",max(bin_ratio$ratio))

save_plot("./results/Figures/Figure2A.pdf", freq_hist.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2A.pdf")

# --------------------------------- Figure 2B ---------------------------------
#   The frequency distribution of NS mutations in antigenic epitopes vs.
#   The frequency of NS mutations in non antigenci epitopes on HA
# -----------------------------------------------------------------------------

antigenic_counts<-antigenic %>% group_by(mutation) %>% summarize(found = length(mutation))
print(knitr::kable(filter(antigenic,mutation %in% 
               antigenic_counts$mutation[antigenic_counts$found>1])))

min.qual<-left_join(min.qual,antigenic,
                    by = c("HOUSE_ID","ENROLLID","SPECID","mutation","pcr_result",
                           "vaccination_status","DPI","Ref_AA","Var_AA"))
min.qual$Antigenic[is.na(min.qual$Antigenic)]<-"None"
Figure2B_data<-filter(min.qual,class_factor=="Nonsynonymous",chr=="HA") %>% 
  select(HOUSE_ID,ENROLLID,SPECID,season,pcr_result,DPI,H3_pos,Ref_AA,Var_AA,Antigenic,freq.var=freq.var.x)

write.csv(Figure2B_data,"./results/Figures/data/Figure2B_data.csv") 


Figure2B<-ggplot(Figure2B_data, aes(y=freq.var,x=Antigenic=="None"))+
  geom_quasirandom(varwidth = TRUE)+
  stat_summary(fun.data="plot.median", geom="errorbar", colour="red", width=0.55, size=0.5)+
  ylab("iSNV frequency")+xlab(label = "")+
  scale_x_discrete(labels = c("Antigenic site","Nonantigenic site"))+
  scale_y_continuous(limits=c(0,0.5))
Figure2B

save_plot("./results/Figures/Figure2B.pdf", Figure2B,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2B.pdf")

Figure2B_stats<-wilcox.test(
  Figure2B_data$freq.var[which(Figure2B_data$Antigenic!="None")],
  Figure2B_data$freq.var[which(Figure2B_data$Antigenic=="None")],
  alternative = "greater")

write_to_summary("Antigenic freq vs. Nonantigenic freq:",
                 Figure2B_stats$p.value)

# --------------------------------- Figure 2C ---------------------------------
#   The frequency of NS mutations in HA overtime at the global level
# -----------------------------------------------------------------------------

# We found this one twice.
nextflu %>% select(-HA1.307R.1)%>% gather(mutation,frequency,-x) ->nextflu.l 

nextflu.l %>% mutate(mutation=gsub("(.*)\\.(.*)","\\1:\\2",mutation))->nextflu.l


minor_HA.ns %>% filter(pcr_result=="A/H3N2", !is.na(Antigenic)) ->H3_minor

H3_minor %>% subset(H3_name %in% nextflu.l$mutation & pcr_result=="A/H3N2",
                    select = c(collect,SPECID,H3_name,freq.var,Antigenic,season)) %>%
  mutate(collect=decimal_date(collect))->collection_points

print(knitr::kable(collection_points))


nextflu.l<-inner_join(nextflu.l,collection_points,by=c("mutation" = "H3_name"))

# Filter to just variants that reach above 5% in or after time of collection.
# Verify the 95% cut is a good idea to remove 307R from figure
 
interesting<- nextflu.l%>% group_by(mutation) %>% filter(x>collect) %>%
  summarize(above=any(frequency>0.05 & frequency<0.95))


nextflu.l<-filter(nextflu.l,mutation %in% interesting$mutation[interesting$above==T])
collection_points<-filter(collection_points,H3_name %in% interesting$mutation[interesting$above==T])

nextflu.p<-ggplot()+
  geom_line(data = nextflu.l,
            aes(x=x,y=frequency,color=mutation,alpha = x>collect),size=1.2)+
  scale_alpha_discrete(range = c(0.2,1))+#theme(legend.position = "none")+
  ylab("Global Frequency")+xlab("Year")+
  scale_x_continuous(breaks = seq(2005,2018,by=2))+
  geom_vline(xintercept =collection_points$collect,linetype=2,alpha=0.5)+
  scale_color_manual(values =
                       c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                         "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Allele")

save_plot("./results/Figures/Figure2C_legend.pdf", nextflu.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2C_legend.pdf")

nextflu.p<-direct.label(nextflu.p,c("last.qp"))+theme(legend.position = "none")

write.csv(nextflu.l,"./results/Figures/data/Figure2C_data_frequency.csv") 
write.csv(collection_points,"./results/Figures/data/Figure2C_data_collection_times.csv") 

save_plot("./results/Figures/Figure2C_label.pdf", nextflu.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2C_label.pdf")

# --------------------------------- Figure 2D ---------------------------------
#   The distribution of sampling times for longitudinal sample pairs
# -----------------------------------------------------------------------------
# Remove mixed infections
intra<- intra %>% filter(!(SPECID2 %in% c("HS1530","MH8137","MH8390")) &
                           !(SPECID1 %in% c("HS1530","MH8137","MH8390")))
intra %>% mutate(DPS1 = collect1-onset,DPS2 = collect2-onset) ->intra
#intra<-subset(intra,freq1<0.5) # Just to make sure


intra %>% group_by(ENROLLID) %>% 
  summarize(onset = unique(onset),
            within_host_time = unique(within_host_time),
            DPS1 = unique(DPS1),
            DPS2 = unique(DPS2)) -> intra_meta

intra_meta<-intra_meta[order(intra_meta$DPS1,intra_meta$DPS2,decreasing = T),]
intra_meta <- intra_meta %>% 
  mutate(DPS2 = ifelse(DPS1==DPS2,yes = DPS2+0.3,no = DPS2)) 
intra_meta$sort_order<-1:nrow(intra_meta)
fig_2D<-ggplot(intra_meta,aes(x = DPS1,xend=DPS2,y = sort_order,yend=sort_order))+
  geom_segment(color = cbPalette[1])+ylab("")+
  xlab("Day post symptom onset")+ 
  theme(axis.line.y=element_blank(),axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  geom_point(aes(y=sort_order,x = DPS1),color=cbPalette[1])+
  geom_point(aes(y = sort_order,x = DPS2),color=cbPalette[1]) + 
  scale_x_continuous(breaks = -2:6) 

fig_2D

save_plot("./results/Figures/Figure2D.pdf", fig_2D,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure2D.pdf")

write.csv(rename(intra_meta,day.post.sympotom.onset1=DPS1,
                 day.post.sympotom.onset2=DPS2),"./results/Figures/data/Figure2D.csv")

# --------------------------------- Figure 2E ---------------------------------
#   The within host dynamics in longitudinal sample pairs
# -----------------------------------------------------------------------------

# Just to make sure
intra<-filter(intra,freq1<0.5) 

intra<-intra %>% mutate(Endpoint="Persistent") %>%
  mutate(Endpoint = if_else(freq1==0,"Arisen",Endpoint)) %>%
  mutate(Endpoint = if_else(freq2==0,"Lost",Endpoint))

intra$Endpoint<-factor(intra$Endpoint,levels = c("Persistent","Arisen","Lost"),ordered = T)

intra.plot<-ggplot(intra,aes(x=as.factor(within_host_time),
                                 y=freq2-freq1,
                                 fill = Endpoint))+
  geom_quasirandom(pch=21,color='black',size=2)+
  scale_fill_manual(values=cbPalette[c(1,3,5)],name="")+
  facet_wrap(~class)+
  xlab("Time within host (days)")+ylab("Change in frequency")

  intra.plot

save_plot("./results/Figures/Figure2E.pdf", intra.plot,
          base_aspect_ratio = 1.3,
          base_width = 10)
embed_fonts("./results/Figures/Figure2E.pdf")

write.csv(select(intra,mutation,ENROLLID,DPS1,DPS2,freq1,freq2),
          "./results/Figures/data/Figure2E.csv")


