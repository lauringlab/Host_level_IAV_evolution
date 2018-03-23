require(magrittr)
require(tidyverse)
require(HIVEr)
if(!("package:tidyverse" %in% search())){
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(readr)
  require(rlang)
}
options(warn=1)
# Here we calculate the average coverage for each sequenced sample. As an aside this coverage does not 
# include reads that had gaps at the position in question.
# 
# 
message(
  paste0(
    "##############################################################################\n",
    "############################## Processing coverage ###########################\n",
    "##############################################################################\n"))


#' Read in the coverage files 

cov.files<-c("./data/processed/HK_1/all.coverage.csv",
             "./data/processed/HK_2/all.coverage.csv",
             "./data/processed/HK_6/all.coverage.csv",
             "./data/processed/HK_7/all.coverage.csv",
             "./data/processed/HK_8/all.coverage.csv",
             "./data/processed/cali09/all.coverage.csv",
             "./data/processed/cali09_2/all.coverage.csv",
             "./data/processed/victoria/all.coverage.csv",
             "./data/processed/victoria_2/all.coverage.csv",
             "./data/processed/perth/all.coverage.csv",
             "./data/processed/perth_2/all.coverage.csv")

cov<-read_rbind(cov.files,n=5,cols = list(Id=col_character()))


cov_sample<- cov %>% group_by(run,Id) %>%
  summarize(coverage = mean(coverage))

# I want to elaborate on the Id column. We can do that here.
LAURING_ID_LOOKUP<-data_frame(Id=unique(cov_sample$Id))

LAURING_ID_LOOKUP<-LAURING_ID_LOOKUP %>% tidyr::separate(Id,into=c("LAURING_ID","dup"),
                                                          sep="_",remove=F,fill="right")

# Because I want to be sure
if(any(is.na(LAURING_ID_LOOKUP$LAURING_ID))){
  stop("Something went wrong in the separate command check the fill")
}
cov_sample<-dplyr::left_join(cov_sample,LAURING_ID_LOOKUP,by="Id")

write.csv(cov_sample,file="./data/processed/secondary/average_coverages.csv")

