# Choose sample Ids for each sample in meta_snv_qual and add number of NS and S mutations for each segement.

require(magrittr)
require(tidyverse)


meta <-read_csv("data/processed/secondary/meta_snv_qual.csv")


coverage<-read_csv("data/processed/secondary/average_coverages.csv")

# We only want one sample even when some were sequenced in duplicate.
coverage_best <- coverage %>% group_by(LAURING_ID) %>%
  filter(coverage==max(coverage))
  

meta<-left_join(meta,coverage_best,by = "LAURING_ID")


write.csv(meta,"./data/processed/secondary/meta_for_ns.s_calc.csv")

# NS and S counts for genes
# 
qual<-read_csv("./data/processed/secondary/qual.snv.csv",
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               )) # read in quality variant calls from all 
qual<-filter(qual,freq.var<0.5)

qual_OR<- qual %>% select(Class,OR) %>% rowwise%>% mutate(Class = gsub("'","",Class),
                                                 Class = gsub("[","",Class,fixed = T),
                                                Class =gsub("]","",Class),
                                                OR = gsub("'","",OR),
                                                OR = gsub("[","",OR,fixed = T),
                                                OR =gsub("]","",OR))

qual_OR<-qual_OR %>% separate(Class,c("Class1","Class2"),sep =",") %>% separate(OR,c("OR1","OR2"),sep = ',') 

# we loose some of the OR here.
qual_out <- tibble(Class = c(qual_OR$Class1,qual_OR$Class2[!is.na(qual_OR$Class2)]),
                       OR = c(qual_OR$OR1,qual_OR$OR2[!is.na(qual_OR$OR2)]))

qual_out<- qual_out %>% filter(!is.na(Class))


gene_counts <- qual_out %>% group_by(OR) %>%
  summarise(NS=length(which(Class=="Nonsyn")),
            S = length(which(Class=="Syn")))



write.csv(gene_counts,"./data/processed/secondary/mut_count_gene.csv")
