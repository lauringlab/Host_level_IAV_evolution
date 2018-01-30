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

qual_OR<-qual_OR %>% separate(Class,c("Class1","Class2"),sep =", ") %>% separate(OR,c("OR1","OR2"),sep = ', ') 

# we loose some of the OR here.
qual_out <- tibble(Class = c(qual_OR$Class1,qual_OR$Class2[!is.na(qual_OR$Class2)]),
                       OR = c(qual_OR$OR1,qual_OR$OR2[!is.na(qual_OR$OR2)]))


gene_counts <- qual_out %>% group_by(OR) %>%
  summarise(NS=length(which(Class=="Nonsyn")),
            S = length(which(Class=="Syn")))


sites<-read_csv("./data/processed/secondary/NS_S_site.csv")

sites<-rename(sites,"OR"="X1","NS_sites" = "NS","S_sites"="S")
sites<-select(sites,OR,NS_sites,S_sites)

####### REMOVE UNWANTED GENES NOW ###########
######
######
all_data <- left_join(gene_counts,sites)

all_data <- all_data %>% mutate(pn= NS/NS_sites,ps = S/S_sites,
                                dn = -(3/4)*log(1-(4*pn)/3),
                                ds = -(3/4)*log(1-(4*ps)/3),
                                DnDs = dn/ds)

# conconical OR

c_OR<-all_data %>% filter(!(OR %in% c("PA-X","PB1-F2")))

c_sum <- c_OR[,c(2:5)] %>% colSums  
c_sum<-  as.data.frame(c_sum)
c_sum<-as.data.frame(t(c_sum))
c_sum <- c_sum %>% mutate(pn= NS/NS_sites,ps = S/S_sites,
                          dn = -(3/4)*log(1-(4*pn)/3),
                          ds = -(3/4)*log(1-(4*ps)/3),
                          DnDs = dn/ds)
