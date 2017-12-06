require(tidyverse)
require(knitr)
require(grid)
require(doMC)
require(magrittr)
doMC::registerDoMC(cores=4)
####### Write to summary results file ######
write_to_summary<-function(line_pattern,value){
  file = readLines("./results.table.tsv")
  line_pattern_regex = paste0("^",line_pattern)
  line = grep(line_pattern_regex,file)
  file[line] = paste0(line_pattern,"\t",value)
  writeLines(file,"./results.table.tsv")
}
######## Coverage #########
#setwd("./notebook/")
slideFunct <- function(data, window, step){ #dapted from http://coleoguy.blogspot.com/2014/04/sliding-window-analysis.html
  coverage<-data$coverage
  concat.pos<-data$concat.pos
  # c.pos<-data$concat.pos
  total <- length(coverage)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- data.frame(mean=rep(F,times=length(spots)),concat.pos=rep(F,times=length(spots)))
  for(i in 1:length(spots)){
    result$mean[i] <- mean(coverage[spots[i]:(spots[i]+window)])
    result$concat.pos[i]<-mean(concat.pos[spots[i]:(spots[i]+window)])
    #first = (spots[i]-1)*window+1 + min(concat.pos)
    #last = first + window
    #result$concat.pos[i]<-mean(c(first,last))
  }
  return(result)
}

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
    adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=100))) %>%
            mutate(ends = ifelse(starts+200<last,starts+200,last)) %>%
                     rowwise() %>%
                     mutate(concat.pos = mean(c(starts,ends)))  -> setup
  ddply(cov.df,~Id+run,slide,setup) -> cov.slid.df

 # cov.slid.df<-ddply(cov.df,~Id+chr,function(x) slide(x,setup)
  x.labels<-ddply(cov.slid.df,~chr,plyr::summarize,concat.pos=concat.pos[which(abs(concat.pos-mean(concat.pos))==(min(abs(concat.pos-mean(concat.pos)))))]) # No idea if this will get what I want but heres to hoping!
  
  #of course sometimes there are 2 good choices I'll take the first one
  x.labels<-ddply(x.labels,~chr,function(x) return(x[1,]))
  
  #x.labels$chr[x.labels$chr %in% c("NR","N_A")]<-"NA"
  
  cov.plot<-ggplot(cov.slid.df, #subset(cov.slid.df,!(Sample%in%c("90","91","93"))),
                   mapping=aes(x=as.factor(concat.pos),
                               y=mean))+geom_boxplot(fill="white")
  
  cov.plot<-cov.plot+ggtitle(title)+ylab("Read depth")+scale_x_discrete(labels = x.labels$chr,breaks=x.labels$concat.pos)+xlab("Concatenated Genome Position")
  cov.plot<-cov.plot+theme(axis.title.y = element_text(vjust=1.2))
  cov.plot<-cov.plot+theme(legend.position="none")  
  return(cov.plot)
}

######## Procesing #########

# proofed JT = 3-24-2017
read_rbind<-function(list,list2=NULL){
  out<-data.frame()
  for(i in 1:length(list)){
    print(paste0("reading in ",list[i]))
    
    x<-read.csv(list[i],stringsAsFactors = F)
    if(is.null(list2)==F){
      x$run=list2[i]
      print(paste0("appending run column: ",list2[i]))
    }
    out<-rbind(x,out)
  }
  return(out)
}


##### join duplicates/ high quality #####

# Sift_dups, join_dups, and quality proofed by JT 3/24/2017
sift_dups<-function(df){

  try(if(nrow(df)>2) stop(paste0("This mutation is found more than once in this sample : ", unique(df$LAURING_ID), " - ", unique(df$mutation))))
  if(nrow(df)==2){ #it's found in both duplicates
    df<-mutate(df,cov=cov.tst.fw+cov.tst.bw)
    higher_qual<-subset(df,cov==max(df$cov))
    if(nrow(higher_qual)>1){ # same cov
      higher_qual<-higher_qual[1,]
    }
    return(subset(higher_qual,select=-c(cov))) # only returns this if the mutation was found twice.
  } 
}  


join_dups.2<-function(df){
	pos.df<-ddply(df,~mutation,sift_dups)
	if(nrow(pos.df)==0){
		sample<-unique(df$LAURING_ID)
		pos <- unique(df$pos)
		chr <- unique(df$chr)
		warning('All bases were removed from sample at position chr:pos')
	}
	return(pos.df)
}

join_dups<-function(df){
 ddply(df,~chr+pos,join_dups.2,.parallel = TRUE)
}

quality<-function(df){ # this now is written to handle 1 isolate (with possibly 2 sequencing runs at a time)
  #good<-subset(df,gc_ul>=1e5)
  # check only one run/sample
  runs<-ddply(df,~LAURING_ID,summarize,runs=length(unique(run))) # LAURING_ID is unique to the sample. How many times was each sample sequenced? duplicates were sequenced in separate runs
  
  if(runs$runs==2 & unique(df$gc_ul)>1e3){
  	if(unique(df$gc_ul)>1e5){
		print(paste0("Samples ",unique(df$LAURING_ID)), " were sequenced twice even though the titer waw ",  unique(df$gc_ul),".It was  treated as duplicates in the analysis")
 	} 
 	dups_good<-join_dups(df) # We only want those sequenced twice with high enough titers.
	return(dups_good) 
 }else if(runs$runs==1 & unique(df$gc_ul)>1e5){
  	return(df)  
  }
  
}

# 
# equal_compare_helper<-function(position){ # take in all the variants found in a sample pair at a given posistion and correct for differences in infering.Meaning the reference base was called here but only in one sample because in the other sample there wasn't a minor allele. Each position with a minor variant should have at least 2 variants,in both samples. after this is done. A minor varaint and the reference. If a position has novariants in either sample then it is excluded from this analysis. In this case only the reference base was present in both samples.  If only major variants are present then the major allele is fixed and there are no other minor variants to infer. 
#   #print(position)
#   pos_sum.1=sum(position$freq1) # What is the sum of all variants at this site in the first sample
#   pos_sum.2=sum(position$freq2) # What is the sum of all variants at this site in the second sample
#   stopifnot(pos_sum.1+pos_sum.2>0)#"No variant at this position - something terrible happend somewhere"
#   
#   
#   if(pos_sum.1==0){ # there isn't a variant in the first sample.The reference is fixed there. the sum of nothing is 0. Try it.  sum() =0.  sum(c()) =0
#     
#     x<-which(position$ref==position$var) # Where are the infered calls
#     stopifnot(length(x)<2) # should be 1 or 0. 
#     if (length(x)==1){ # The reference has been infered in the second sample but there are no varaints in the first so it was not infered
#       position$freq1[x]<-1 # The reference base is fixed here
#     } else if(length(x)==0){ # the reference was not infered in the second sample and no variants here. Add the reference now. There must have been another base fixed in the second sample.
#       extra<-mutate(position,var=ref,mutation=paste0(chr,"_",ref,pos,var),freq1=1,freq2=0) # adding a line with the reference base fixed in sample 1 but not sample 2
#       position<-rbind(position,extra)
#       #print("need to add a few here")
#     }
#     # Do it again if there are no varaints in the second sample.
#   } else if(pos_sum.2==0){ # there isn't a variant in the second sample. The reference is fixed. also I told you so. (line 196)
#     x<-which(position$ref==position$var)
#     stopifnot(length(x)<2) # should be 1 or 0
#     if (length(x)==1){ # The reference has been infered in the first sample but there are no varaints in the second so it was not infered
#       position$freq2[x]<-1
#     } else if(length(x)==0){ # the reference was not infered in the first sample and no variants here. Add the reference now
#       extra<-mutate(position,var=ref,mutation=paste0(chr,"_",ref,pos,var),freq1=0,freq2=1) # adding a line with the reference base fixed in sample 2 but not sample 1
#       position<-rbind(position,extra)
#     }
#   }
#   
#   if(nrow(position)>1){ # if after all this there is only 1 row in the position then the variant base is fixed in both samples and not interesting to us. We don't look at all the sites that have the reference base only in both.
#     return(position)
#   }else{
#     return(position[F,])
#   }
# }
# 
# 
# equal_compare<-function(data){      ## There are cases where we have infered a majority variant in one of the pairs, but in the other it is not infered since the minority variant was not found. Here I'll add the fixed variants. The work is done in the helper function
#   
#   #out<-ddply(data,~SPECID.1+SPECID.2+chr+pos,equal_compare_helper)
#   out<- data %>% group_by(SPECID.1,SPECID.2,chr,pos)%>% do(equal_compare_helper(.))
#   out<-ungroup(out)
#   return(out)
# }
# 
# get_freqs<-function(pairs,snv){ # take in a data frame of pairs of Ids. Only 1 pair, and a list of snv calls. and output the comparison between the 2. each iSNV and its frequency in both samples
#   stopifnot(nrow(pairs)==1) # Verify only 1 pair here
#   snv<-subset(snv,SPECID %in% c(as.character(pairs$SPECID.1),as.character(pairs$SPECID.2))) # just need the snv in the samples we're looking at.
#   if(nrow(snv)>0){ # There are mutations.
#     mut_table <- dcast(snv, mutation ~ SPECID, value.var = "freq.var") # a data frame with mutation down the first row and then frequency in either sample in the next 2.
#     mut_table[is.na(mut_table)] <- 0 # replace NA with 0. It wasn't found brotha!
#     mut_table$SPECID.1<-pairs$SPECID.1 # add column with first sample ID
#     mut_table$SPECID.2<-pairs$SPECID.2 # add column with second sample ID
#     names(mut_table)[which(names(mut_table)==as.character(pairs$SPECID.1))]<-'freq1' # rename this column as the frequency in the first sample 
#     names(mut_table)[which(names(mut_table)==as.character(pairs$SPECID.2))]<-'freq2' # dido
#     pat="([A-Z]+[1-2]?)_([A|C|T|G]{1})([0-9]+)([A|C|T|G]{1})" # string for getting data from mutation name
#     
#     # This function can only be run on samples that qualified for variant identification. 
#     # If no variants were found in the sample then the SPECID will be missing from mut_table column and so
#     # freq1 or freq2 will be missing since nothing was found we set that to 0 here and add the column. 
#     # equal compare will replace these cases with the reference at 1.
#     if(!('freq1' %in% names(mut_table))){
#       mut_table$freq1=0
#     }
#     if(!('freq2' %in% names(mut_table))){
#       mut_table$freq2=0
#     }
#     mut_table=mutate(mut_table,chr=as.character(sub(pat,"\\1",mutation)),ref = sub(pat,"\\2",mutation),pos=as.numeric(sub(pat,"\\3",mutation)),var=sub(pat,"\\4",mutation))  # add chr pos ref and var from mutation string
#     #print(mut_table)
#     all.freq<-equal_compare(mut_table) # fill in differences based on infered. In this case we are left with only sites that are polymorphic in one or between the 2 samples
#     
#     if(nrow(all.freq)>0){ # Some differences exists
#       return(all.freq)
#     }
#     else{ # some snv were initially present (before frequency filter) but after taking into account differences in inference the samples are the same.
#       #return(data.frame(mutation=NA,freq1=NA,freq2=NA,SPECID.1=pairs$SPECID.1,SPECID.2=pairs$SPECID.2,chr=NA,ref=NA,"pos"=NA,var=NA))
#       return(mut_table[F,])
#     }
#   }
#   else{ # No variants found in either sample
#     #return(data.frame(mutation=NA,freq1=NA,freq2=NA,SPECID.1=pairs$SPECID.1,SPECID.2=pairs$SPECID.2,chr=NA,ref=NA,"pos"=NA,var=NA))
#     return(mut_table[F,])
#     }
# }  







