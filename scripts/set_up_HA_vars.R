HA<-read.csv("../results/HA_nonsynom.csv",stringsAsFactors = F)

H3<-subset(HA,subtype=="H3")

require(tidyverse)
require(magrittr)

H3 %>% subset(select=c(ENROLLID,SPECID,freq.var,mutation,AA_pos,H3_pos,H3_seg,PDB_4HMG,Antigenic))

H3 %>% group_by(H3_seg) %>% summarize(paste(H3_pos,collapse="+")) ->x
