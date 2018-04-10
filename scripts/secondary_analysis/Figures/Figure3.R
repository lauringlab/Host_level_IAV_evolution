require(tidyverse)
require(HIVEr)
cbPalette<-wesanderson::wes_palette("Zissou") # Set up figures
require(cowplot)
require(extrafont)
require(VGAM)
require(lubridate)
require(doMC)
doMC::registerDoMC(cores=4)
options(warn=1)
# --------------------------------- Functions ---------------------------------
write_to_summary<-function(line_pattern,value){
  file = readLines("./results/results.table.tsv")
  line_pattern_regex = paste0("^",line_pattern)
  line = grep(line_pattern_regex,file)
  file[line] = paste0(line_pattern,"\t",value)
  writeLines(file,"./results/results.table.tsv")
}

# These are for making the tranmission plot
one_for_each_plotting<-function(df){
  df$offset_x.1<-min(df$donset) #This puts the points at the onset 
  df$offset_x.2<-max(df$donset) # it used to be offset_x 
  df$offset_y.1<-min(df$offset_y[df$donset==min(df$donset)])  # This needs work
  df$offset_y.2<-max(df$offset_y[df$donset==max(df$donset)])
  return(df[1,])
}


make_transmission_plot<-function(x,type="all"){ # type should be all or valid or invalid.
  stopifnot(type %in% c('all','valid','invalid'))
  #x<-plyr::adply(x,1,longform_pairs) # Each row is a pair so for each pair split the row into 2 rows 
  x<-longform_pairs(x)
  x<-subset(x,select=c(HOUSE_ID,pair_id,ENROLLID,pcr_result,season,onset,gc_ul,sequenced,sequenced_pair,titer_pair,snv_qualified_pair,valid,quality_distance)) # just select these columns for each row
  x<-plyr::ddply(x,~HOUSE_ID,function(x){ # This will be useful for ording the plotting - at least that's the hope. For each house get the first onset date and set the donset as the distance of each case from that date.
    min_onset<-min(x$onset)
    mutate(x,min_onset=min_onset,donset=onset-min_onset)
  })
  
  
  x<-x[order(x$min_onset,x$HOUSE_ID,decreasing = T),] # I told you it would be useful - based on the first onset in each house.
  HOUSE_order<-unique(x$HOUSE_ID) # This is the order we want the houses to be in on the plot
  
  
  x<-plyr::ddply(x,~HOUSE_ID,function(z){
    z$y_Id=which(HOUSE_order==unique(z$HOUSE_ID))*2 # where does this house fall in the order
    return(z)})
  
  
  # I don't really use this any more but it is useful if you want to offset points with the same x value by shifting them on the x axis
  x<- plyr::ddply(x, ~onset+y_Id, function(x){ # this applys an offset for different individuals who are sick on the same day
    if(length(unique(x$ENROLLID))>1){ # if there are mulitple people here
      ENROLLIDs<-sort(unique(x$ENROLLID),decreasing=T)
      x<-mutate(x,sort_order=match(ENROLLID,ENROLLIDs)-1,offset_x=as.numeric(donset+sort_order/5))
    }
    else{
      x<-mutate(x,offset_x=as.numeric(donset))
    }
    return(x)
  })
  
  x<- plyr::ddply(x, ~onset+y_Id, function(x){ # this applys an offset for different individuals who are sick on the same day
    if(length(unique(x$ENROLLID))>1){ # if there are mulitple people here
      ENROLLIDs<-sort(unique(x$ENROLLID),decreasing=T)
      x<-mutate(x,sort_order=match(ENROLLID,ENROLLIDs)-1,offset_y=y_Id+((-1*length(ENROLLIDs)+4)*sort_order-1))
    }
    else{
      x<-mutate(x,offset_y=y_Id)
    }
    return(x)
  })
  
  
  # Curves and such
  
  # if the onset date is the same for both add a a curve and aline - this will be interpretted as 2 curves later
  
  x$line=T # start with all lines
  x$curve=F
  x<-plyr::ddply(x,~pair_id,function(x){
    if(length(unique(x$onset))==1){ # This pair is sick on the same day # and they are the only ones in the house
      x$curve=T 
    }
    return(x)})
  
  # If there is point between the two in the pair add a curve and no line
  
  # This is a mess of a block of code. We pass each house grab the y axis. 
  #Then we look at each pair in the house. If the points are separated on 
  #the x axis with another point in between then we get rid of the line 
  #connecting and replace it with a curve
  x<-plyr::ddply(x,~HOUSE_ID,function(u){
    offset_y<-sort(unique(u$offset_y))
    plyr::ddply(u,~pair_id,function(y,offset_y){
      diff<-which(offset_y %in% y$offset_y)
      #print(diff)
      #stopifnot(length(diff)==2)
      if(length(diff)==2 & abs(diff[1]-diff[2])!=1){ # ie there is a point between these two
        #print(abs(diff[1]-diff[2]))
        y$curve=T
        y$line=F
      }
      #print(y)
      return(y)
    },offset_y)
  })
  
  seasons = plyr::ddply(x,~season,summarize,ymin=min(y_Id)-0.25,ymax=max(y_Id)+0.25,middle=mean(y_Id))
  
  if(type=="valid"){
    x<-subset(x,valid==T)
  }
  else if(type=="invalid"){
    x<-subset(x,valid==F)
    x$quality_distance=1 # makes all lines solid otherwise dashed lines represnet those with too big L1-norm
    
  }else if(type=="all"){
    x$quality_distance=1 # makes all lines solid
  }
  x_useful<-subset(x,snv_qualified_pair==T )
  x_lost<-subset(x,snv_qualified_pair==F) # they were lost either to not sequencing or poor titer or poor sequencing
  
  x_lines_useful<-subset(x_useful,line==T & curve==F)
  x_curves_over_useful<-subset(x_useful,line==F & curve==T)
  x_curves_both_useful<-subset(x_useful,line==T & curve==T)
  
  
  x_lines_useful<-plyr::ddply(x_lines_useful,~pair_id,one_for_each_plotting) # back to short form, 1 row/ pair for the lines.  offset_x is the onsets here. 1 is the first case 2 is the second
  x_curves_over_useful<-plyr::ddply(x_curves_over_useful,~pair_id,one_for_each_plotting)
  x_curves_both_useful<-plyr::ddply(x_curves_both_useful,~pair_id,one_for_each_plotting)
  
  arrow_length= 0.005
  
  transmission_plot_useful<-ggplot()+ylab("")+xlab("Onset (relative to index case)")+geom_segment(data=seasons,aes(x=-0.5,xend=-0.5,y=ymax,yend=ymin,group=season))+geom_text(data=seasons,aes(x=-1.3,y=middle,label=season,angle=0))
  if(nrow(x_lines_useful)>0){
    transmission_plot_useful<-transmission_plot_useful+geom_segment(data=x_lines_useful,aes(x=offset_x.1,xend=offset_x.2,y=offset_y.1,yend=offset_y.2,group=pair_id,linetype=factor((-2*quality_distance)+3)), # this makes the lines solid if the quality distance ==T and dashed if not
                                                                    arrow = arrow(length = unit(arrow_length, "npc")),
                                                                    color=cbPalette[1])
  }
  if(nrow(x_curves_over_useful)>0){
    transmission_plot_useful<-transmission_plot_useful+geom_curve(data=x_curves_over_useful,aes(x=offset_x.1,xend=offset_x.2,y=offset_y.1,yend=offset_y.2,group=pair_id,linetype=factor((-2*quality_distance)+3)),
                                                                  curvature = -0.1,
                                                                  arrow = arrow(length = unit(arrow_length, "npc")),
                                                                  color=cbPalette[1])
  }
  if(nrow(x_curves_both_useful)>0){
    transmission_plot_useful<-transmission_plot_useful+geom_curve(data=x_curves_both_useful,aes(x=offset_x.1,xend=offset_x.2,y=offset_y.1,yend=offset_y.2,group=pair_id,linetype=factor((-2*quality_distance)+3)),
                                                                  curvature = -0.4,
                                                                  arrow = arrow(length = unit(arrow_length, "npc")),
                                                                  #linetype=2,
                                                                  color=cbPalette[1])+
      geom_curve(data=x_curves_both_useful,aes(x=offset_x.2,xend=offset_x.1,y=offset_y.2,yend=offset_y.1,group=pair_id,linetype=factor((-2*quality_distance)+3)),
                 curvature = -0.4,
                 arrow = arrow(length = unit(arrow_length, "npc")),
                 #linetype=2,
                 color=cbPalette[1])
  }
  transmission_plot_useful<-transmission_plot_useful+geom_point(data=x_useful,aes(x=as.numeric(donset),y=offset_y),size=0.5)
  if(nrow(x_lost)>0){ # if at least some the data is removed then we'll add the differently colored lines
    x_lines_lost<-subset(x_lost,line==T & curve==F)
    x_curves_over_lost<-subset(x_lost,line==F & curve==T)
    x_curves_both_lost<-subset(x_lost,line==T & curve==T)
    
    # Now we just need one x axis point for each onset in the pair
    x_lines_lost<-plyr::ddply(x_lines_lost,~pair_id,one_for_each_plotting)
    x_curves_over_lost<-plyr::ddply(x_curves_over_lost,~pair_id,one_for_each_plotting)
    x_curves_both_lost<-plyr::ddply(x_curves_both_lost,~pair_id,one_for_each_plotting)
    
    transmission_plot_all<-transmission_plot_useful+geom_point(data=x_lost,aes(x=donset,y=offset_y),size=0.5)+
      geom_segment(data=x_lines_lost,aes(x=offset_x.1,xend=offset_x.2,y=offset_y.1,yend=offset_y.2,group=pair_id),
                   arrow = arrow(length = unit(arrow_length, "npc")),
                   color=cbPalette[4])
    if(nrow(x_curves_over_lost)>0){
      transmission_plot_all<-transmission_plot_all+geom_curve(data=x_curves_over_lost,aes(x=offset_x.1,xend=offset_x.2,y=offset_y.1,yend=offset_y.2,group=pair_id),
                                                              curvature = -0.1,
                                                              arrow = arrow(length = unit(arrow_length, "npc")),
                                                              color=cbPalette[4])
    }
    if(nrow(x_curves_both_lost)>0){
      transmission_plot_all<- transmission_plot_all+geom_curve(data=x_curves_both_lost,aes(x=offset_x.1,xend=offset_x.2,y=offset_y.1,yend=offset_y.2,group=pair_id),
                                                               curvature = -0.4,
                                                               arrow = arrow(length = unit(arrow_length, "npc")),
                                                               linetype=2,
                                                               color=cbPalette[4])+
        geom_curve(data=x_curves_both_lost,aes(x=offset_x.2,xend=offset_x.1,y=offset_y.2,yend=offset_y.1,group=pair_id),
                   curvature = -0.4,
                   arrow = arrow(length = unit(arrow_length, "npc")),
                   linetype=2,
                   color=cbPalette[4])+
        ylab("")+xlab("Onset (relative to index case)")
    }
  }  
  else{
    transmission_plot_all<-transmission_plot_useful
  }
  return(transmission_plot_all+theme(legend.position = 'none',axis.line.y=element_blank(),axis.ticks.y = element_blank()) +scale_y_continuous(breaks=c())+scale_x_continuous(breaks=c(0:11),limits = c(-1.5,11)))
}
# PA probability
Pt_PA<-function(x,l,max_Nb){
  s<-0
  for(i in 1:max_Nb){
    c<- ((1-x)^i) *  (l^i)/((exp(l)-1)*factorial(i) )
    s<-s+c
  }
  return(1-s)
}

Pt_BetaBin<-function(x,l,max_Nb){
  s<-0
  for(i in 1:max_Nb){
    prob_not_T<- L.Nb.beta(v_r=1,v_d=(1-x),Nb=i,gc_ul=1e4,
                           threshold=0.02,accuracy_stringent)
    prob_Nb<-(l^i)/((exp(l)-1)*factorial(i) )
    c<-prob_not_T*prob_Nb
    s<-s+c
  }
  return(1-s)
}

Pt_BetaBin<-Vectorize(Pt_BetaBin,vectorize.args = "x")


straight_sum <-function(data,max_nb){
  Nb<-data$Nb[which(data$LL==max(data$LL))] # Get the  max lambda
  good_range<-subset(data,LL> (max(LL)-1.92)) # get the bottlenecks that fall in this region the 95% confidence intereval
  lower<-good_range$Nb[1]
  upper <- good_range$Nb[nrow(good_range)]
  if(length(Nb)>1){
    Nb=max(Nb)
  }
  if(Nb==max_nb){
    return(tibble(Nb=NA,lower_95=NA,
                  upper_95=NA))
  }else{
    return(tibble(Nb=Nb,lower_95=lower,
                  upper_95=upper))
  }
}
# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------

# A
possible_pairs.dist<-read_csv("./data/processed/secondary/possible.pairs.dist.csv")
#B
all_pairs.tp <- read_csv("./data/processed/secondary/transmission_pairs.csv")
#C
trans_freq <- read_csv("./data/processed/secondary/trans_freq.csv")
#D
meta<-read_csv("./data/reference/all_meta.sequence_success.csv")

trans_freq.comp<-read_csv("./data/processed/secondary/transmission_pairs_freq.poly.donor.csv")
#E
accuracy_stringent<-read.csv("./data/reference/accuracy_stringent.csv",stringsAsFactors = F)


# Supplemental
intra<-read_csv("./data/processed/secondary/Intrahost_all.csv")
qual<-read_csv("./data/processed/secondary/qual.snv.csv",
               col_types = list(
                 ENROLLID= col_character(),
                 SPECID = col_character(),
                 LAURING_ID = col_character(),
                 Id = col_character()
               ))

# --------------------------------- Figure 3A ---------------------------------
#   L1 -norm genetic distance between transmission pairs
# -----------------------------------------------------------------------------
# get the percentiles for the community pairs
cutoffs<-as.data.frame(
  quantile(possible_pairs.dist$L1_norm[possible_pairs.dist$Household==F],
           probs = seq(0,1,0.05))) 
names(cutoffs)<-"L1_norm"
cutoffs$threshold<-seq(0,1,0.05)
cutoffs<-cutoffs %>% rowwise() %>%
  mutate(
    valid_pairs=nrow(possible_pairs.dist[(possible_pairs.dist$valid==T & 
                                            possible_pairs.dist$L1_norm<L1_norm),])
    )


#ggplot(cutoffs,aes(y=valid_pairs,x=threshold))+geom_point()
print(knitr::kable(cutoffs))

# adding whether or not a pair is a valid pair
possible_pairs.dist<-mutate(possible_pairs.dist,
                            quality_distance=
                              L1_norm<cutoffs$L1_norm[cutoffs$threshold==0.05])


figure3A_data<-subset(possible_pairs.dist,valid==T | Household==F)
figure3A_data<-figure3A_data %>%
  select(season,SPECID1,SPECID2,L1_norm,valid,Household)
l1norm.p_tp<-ggplot(figure3A_data,
                    aes(x=L1_norm,fill=as.factor((valid-1)*-1),y=..ncount..))+
  geom_histogram(color='white',binwidth = 7.5,boundary=0,position = 'dodge')+
  scale_fill_manual(name="",labels=c("Household transmision","Community pairs"),values = cbPalette[c(1,4)])+
  xlab("L1 Norm")+ylab("Normalized Count")+
  theme(legend.position = c(0.5,0.5))+
  geom_segment(aes(x=15,xend =cutoffs$L1_norm[cutoffs$threshold==0.05],y=0,yend=1),
               linetype=2,color=cbPalette[5],size=0.3)
l1norm.p_tp

write.csv(figure3A_data,"./results/Figures/data/Figure3A.csv")

save_plot("./results/Figures/Figure3A.pdf", l1norm.p_tp,
          base_aspect_ratio = 1.3,
          base_height = 6)
embed_fonts("./results/Figures/Figure3A.pdf")

# --------------------------------- Figure 3B ---------------------------------
#   Visulizing tranmision pairs
#   ransmission rules. These apply to all cases where 2 individuals are sick within
#    the same household within a week of eachtoher (difference in date of 
#    onset <= 7 days).In the event of multiple possible donors we assume the donor 
#    is the individual with symptom onset neast to the recipeient. The donor and 
#    recipeint are never have symptoms on the same day unless they are the only cases 
#    in the house. In this case we will randomize the pair and estimate a bottleneck 
#    both ways. We'll take a look at the pairs we are refering to and then get
#     a list of all those that quailify.
# -----------------------------------------------------------------------------
# Add distance cut off data to all_pairs.tp

all_pairs.tp<-merge(all_pairs.tp,subset(possible_pairs.dist,
                                        select=c(ENROLLID1,ENROLLID2,quality_distance))
                    ,all.x = T)
all_pairs.tp$quality_distance[all_pairs.tp$snv_qualified_pair==F]<-NA

# all_seasons<-make_transmission_plot(all_pairs.tp,type='valid')
# 
# all_seasons

all_seasons_useful<-make_transmission_plot(
  subset(all_pairs.tp,snv_qualified_pair==T & valid==T),
  type='valid')
all_seasons_useful

write.csv(subset(all_pairs.tp,snv_qualified_pair==T & valid==T,
                 select=c(HOUSE_ID,ENROLLID1,ENROLLID2,
                          onset1,onset2,quality_distance,double)),
          "./results/Figures/data/Figure3B.csv")

save_plot("./results/Figures/Figure3B.pdf", all_seasons_useful,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure3B.pdf")

# --------------------------------- Figure 3C ---------------------------------
#   The frequency of variants in transmission pairs. AKA the TV plot
# -----------------------------------------------------------------------------


trans_freq.p<-ggplot(trans_freq,aes(x=freq1,y=freq2))+
  geom_point()+xlab("Frequency in donor")+
  ylab("Frequency in recipient")
trans_freq.p

save_plot("./results/Figures/Figure3C.pdf", trans_freq.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure3C.pdf")

write.csv(x = select(trans_freq,SPECID1,SPECID2,freq1,freq2),
          "./results/Figures/data/Figure3C.csv")

# --------------------------------- Figure 3D ---------------------------------
#   The fit of the PA transmission model
#   The error bars are as in Sobel Lenard et al assuming a binomial distribution
#   with sample size of iSNV in window and probability of success estimated as
#   proportion transmitted.
# -----------------------------------------------------------------------------

# Add gc_ul
trans_freq.comp<-mutate(trans_freq.comp,gc_ul1 = meta$gc_ul[match(SPECID1,meta$SPECID)],
                        gc_ul2 = meta$gc_ul[match(SPECID2,meta$SPECID)])
# Fit the model
pa_total_fit<-trans_fit(trans_freq.comp,Nb_max=100,model="PA",
                        threshold=NULL,acc=NULL,pair_id)

counts<-trans_freq.comp %>% group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1>0 & freq1<0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants,
         weight_factor = 1)

#pa_total_fit<-pa_total_fit %>% filter(pair_id!=184)
require(bbmle)


# zero_truncated poisson
zdpois_fit<-dist_prob_wrapper(ddist = "dzpois",params = "lambda")
# 
# Running with LL_D_given_l gives an l of 1.15 which is
# very close to the 1.12 that was estimated earlier.


dzpois_model_fit<-bbmle::mle2(minuslogl = zdpois_fit,start = list(lambda = 1),
                data = list(data = pa_total_fit,
                            weight = counts))
 
# nb_fit<-dist_prob_wrapper("dposnegbin","size,mu")
# 
# # Just to get an idea of where the max is.
# size = seq(1,20,1)
# mu = seq(0.01,.95,0.05)
# nsize = rep(size,each = length(mu))
# nmu = rep(mu,times = length(size))
# x<-purrr::map2(nsize,nmu,.f = function(x,y) nb_fit(data = pa_total_fit,weight = counts,size = x,mu = y))
# 
# nb_output<-tibble(size = nsize,mu = nmu,negLL = unlist(x))
# 
# ggplot(nb_output,aes(x = size,y= mu,fill = -negLL,z = -negLL))+
#   geom_tile()+geom_contour(binwidth = 1)+
#   geom_point(data = filter(nb_output,negLL ==min(negLL)))
# 
# 
# 
# nb_model_fit<-bbmle::mle2(minuslogl = nb_fit,start = list(size=100,prob=0.5),
#                               data = list(data = pa_total_fit,
#                                           weight = counts),
#                           method = "Nelder-Mead",
#                           hessian = TRUE,
#                           control=list(trace=TRUE))#, maxit=5000))

mean_zpois<-function(l) l/(1-exp(-1*l))
con_int<-confint(dzpois_model_fit)
write_to_summary("P-A Nb:",mean_zpois(dzpois_model_fit@coef))
write_to_summary("P-A lambda:",dzpois_model_fit@coef)
write_to_summary("P-A CI:",paste(bbmle::confint(dzpois_model_fit),collapse = "-"))


prob_above_5 <-1-sum(dzpois(c(1,2,3,4,5),dzpois_model_fit@coef))
write_to_summary("P-A prob >5",prob_above_5)

# Plot the fit
# Sliding window with window of w step of step
w<-0.05
step = 0.025
windows<-tibble(s=seq(0.02,1-w,by=step),end=s+w)

out <- windows %>% rowwise() %>%
  mutate(iSNV = nrow(filter(trans_freq.comp,freq1>=s,freq1<end)),
         transmitted = nrow(filter(trans_freq.comp,freq1>=s,freq1<end,found==T)),
         freq = mean(trans_freq.comp$freq1[which(trans_freq.comp$freq1>=s & trans_freq.comp$freq1<end)]),
         prob = transmitted/iSNV,
         error_bottom = qbinom(c(0.025),iSNV,prob)/iSNV,
         error_top = qbinom(c(0.975),iSNV,prob)/iSNV,
         many = iSNV>5)

model<-tibble(s = seq(0,1,0.01))
model<- mutate(model,prob = Pt_PA(s,dzpois_model_fit@coef,100),
               lower = Pt_PA(s,con_int[1],100),
               upper = Pt_PA(s,con_int[2],100))

window_data.p<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model,aes(x=s,y=prob),color=cbPalette[5])+
  geom_ribbon(data=model,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=cbPalette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,
             aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))

#save_plot("./results/Figures/Figure3D_nodots.pdf", window_data.p,
#          base_aspect_ratio = 1.3)
#embed_fonts("./results/Figures/Figure3D_nodots.pdf")



save_plot("./results/Figures/Figure3D.pdf", window_data.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure3D.pdf")

write.csv(x = select(trans_freq.comp,SPECID1,SPECID2,freq1,freq2),
          "./results/Figures/data/Figure3D_points.csv")
write.csv(x = model,
          "./results/Figures/data/Figure3D_model.csv")
#bottle neck of 10

model_10<-tibble(s = seq(0,1,0.01))
model_10<- mutate(model,prob = Pt_PA(s,10,100))
window_data.p_10<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model,aes(x=s,y=prob),color=cbPalette[5])+
  geom_ribbon(data=model,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=cbPalette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,
             aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))+
  geom_line(data=model_10,aes(x=s,y=prob),color=cbPalette[1])
save_plot("./results/Figures/Figure3D_10.pdf", window_data.p_10,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure3D_10.pdf")

write.csv(x = model_10,
          "./results/Figures/data/Figure3D_model_10.csv")

# --------------------------------- Figure 3E ---------------------------------
#   The fit of the BetaBin transmission model
# -----------------------------------------------------------------------------
# Add gc_ul
trans_freq.comp<-mutate(trans_freq.comp,gc_ul1 = meta$gc_ul[match(SPECID1,meta$SPECID)],
                        gc_ul2 = meta$gc_ul[match(SPECID2,meta$SPECID)])
# Fit the model
beta_total_fit<-trans_fit(subset(trans_freq.comp,freq1<0.5),
                          Nb_max=100,model="BetaBin",
                          threshold=0.02,acc=accuracy_stringent,
                          pair_id)

zdpois_fit<-dist_prob_wrapper(ddist = "dzpois",params = "lambda")

counts<-trans_freq.comp %>% group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1>0 & freq1<0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants,
         weight_factor = 1)

dzpois_model_fit_bb<-bbmle::mle2(minuslogl = zdpois_fit,start = list(lambda = 1),
                              data = list(data = beta_total_fit,
                                          weight = counts))
conf_int_BB<-bbmle::confint(dzpois_model_fit_bb)
summary(dzpois_model_fit_bb)

write_to_summary("BB Nb:",mean_zpois(dzpois_model_fit_bb@coef))
write_to_summary("BB lambda:",dzpois_model_fit_bb@coef)
write_to_summary("BB CI:",paste(conf_int_BB,collapse = "-"))

# Plot the fit
# Sliding window with window of w step of step from above

model_betaBin<-tibble(s = seq(0,1,0.01))
model_betaBin<- model_betaBin %>% rowwise() %>%
  mutate(prob = Pt_BetaBin(s,dzpois_model_fit_bb@coef,100),
               lower = Pt_BetaBin(s,conf_int_BB[1],100),
               upper = Pt_BetaBin(s,conf_int_BB[2],100))

window_data_bb.p<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model_betaBin,aes(x=s,y=prob),color=cbPalette[5])+
  geom_ribbon(data=model_betaBin,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=cbPalette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,
             aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))

# save_plot("./results/Figures/Figure3E_nodots.pdf", window_data_bb.p,
#          base_aspect_ratio = 1.3)
# embed_fonts("./results/Figures/Figure3E_nodots.pdf")



save_plot("./results/Figures/Figure3E.pdf", window_data_bb.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure3E.pdf")





write.csv(x = select(trans_freq.comp,SPECID1,SPECID2,freq1,freq2),
          "./results/Figures/data/Figure3E_points.csv")
write.csv(x = model_betaBin,
          "./results/Figures/data/Figure3E_model.csv")


model_10_bb<-tibble(s = seq(0,1,0.01))
model_10_bb<- model_10_bb%>% rowwise() %>% 
  mutate(prob = Pt_BetaBin(s,10,100))
window_data.p_10<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model,aes(x=s,y=prob),color=cbPalette[5])+
  geom_ribbon(data=model,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=cbPalette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,
             aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))+
  geom_line(data=model_10_bb,aes(x=s,y=prob),color=cbPalette[1])

save_plot("./results/Figures/Figure3E_10.pdf", window_data.p_10,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure3E_10.pdf")

write.csv(x = model_10_bb,
          "./results/Figures/data/Figure3E_model_10.csv")

# --------------------------------- Supplemental Table ---------------------------------
#   The fit of the binomial model for each pair.
# -----------------------------------------------------------------------------


max_nb<-200
beta_Nb<-trans_fit(subset(trans_freq.comp,freq1<0.5),
                   Nb_max=max_nb,model="BetaBin",
                   threshold=0.02,acc=accuracy_stringent,
                   pair_id)


beta_Nb_sum<- beta_Nb %>% group_by(pair_id) %>% do(straight_sum(.,max_nb))
beta_Nb_sum<-beta_Nb_sum[order(beta_Nb_sum$Nb,decreasing = T),]



beta_t<-beta_Nb_sum %>% mutate(
  CI = paste0(lower_95,"-",upper_95),
  lambda = paste0(Nb," (",CI,")")
) %>%
  select(pair_id,Nb,pair_id,CI)
beta_t<-left_join(beta_t,
                  select(trans_freq.comp,pair_id,ENROLLID1,ENROLLID2,
                         collect1,collect2,transmission,SPECID1,SPECID2))%>%
  distinct(pair_id, .keep_all = TRUE) %>%
  dplyr::rename(donor_sample = collect1, recipient_sample=collect2,
                estimated_transmission_date=transmission) 

snv_data<- trans_freq.comp %>% group_by(pair_id,pcr_result) %>%
  summarize(minority_isnv=length(which(freq1<0.5)),
            transmitted_minority_isnv =length(which(freq1<0.5&found==T)))
beta_t<-left_join(beta_t,snv_data)
beta_t <- mutate(beta_t,
                 within_host_time=abs(donor_sample-estimated_transmission_date)+
                   abs(estimated_transmission_date-recipient_sample),
                 proportion_transmitted = transmitted_minority_isnv/minority_isnv)

# Adding the ages. 
ages<- read_csv("./data/reference/HIVE_ages_by_season.csv",
                col_types = cols(STUDY_ID=col_character()))%>%
  mutate(DOB = parse_date_time(DOB,c("db!y","%m/%d/%Y"))) %>%
  mutate(DOB=if_else(condition = (as.numeric(as.POSIXct(today())-DOB)/365)<AGEYR,
                     true = DOB-years(100),false=DOB))


beta_t<-left_join(beta_t,select(ages,STUDY_ID,DOB),by=c("ENROLLID1"="STUDY_ID")) %>% 
  mutate(Donor_age=
           as.numeric(as.POSIXct(estimated_transmission_date)-DOB)/365.2425) %>% 
  select(-DOB) %>%
  left_join(.,select(ages,STUDY_ID,DOB),by=c("ENROLLID2"="STUDY_ID")) %>%
  mutate(Recipient_age=
           as.numeric(as.POSIXct(estimated_transmission_date)-DOB)/365.2425) %>% 
  select(-DOB)
out_beta_t<-beta_t %>% ungroup() %>%
  select(Nb,CI,Subtype=pcr_result,donor_sample,recipient_sample,estimated_transmission_date,
         Donor_age,Recipient_age,minority_isnv,transmitted_minority_isnv,ENROLLID1,ENROLLID2,SPECID1,SPECID2)
out_beta_t$Nb[is.na(out_beta_t$Nb)]<-">200"
write.csv(out_beta_t,"./data/processed/secondary/beta_bottlenecks_by_pair.csv")


# --------------------------------- Supplemental Table ---------------------------------
#   The fit of the binomial model for every possible sample pairing.
# --------------------------------------------------------------------------------------


intra %>% mutate(DPS1 = collect1-onset,DPS2 = collect2-onset) ->intra
intra<-subset(intra,freq1<0.5) # Just to make sure

trans<-trans_freq %>% filter(freq1<0.5) %>%
  mutate(Endpoint="Persistent") %>%
  mutate(Endpoint = if_else(freq1==0,"Arisen",Endpoint)) %>%
  mutate(Endpoint = if_else(freq2==0,"Lost",Endpoint))
trans$Endpoint<-factor(trans$Endpoint,levels = c("Persistent","Arisen","Lost"),ordered = T)


# Now I will work to get longitudinal pairs for people in the transmission subset. 
# The first step is to subset the transmission data into just the meta data for each pair
# and rename the columns according to whom the data comes from. The data starts as 1 row
# for each mutation so the distinct calapses this.

trans_long<- trans %>% select(HOUSE_ID,Donor_ENROLLID=ENROLLID1,
                              Recipient_ENROLLID=ENROLLID2,
                              onset1,onset2,
                              transmission) %>% distinct()

# Now we will get the SPECID for the donor by joining with the meta data
# we are only interested in samples that qualified for snv identification.
# NB : left join includes all combinations so if there are two SPECID in meta, 
# which is present as two rows, then both are kept, again as two rows. Key columns 
# added and the SPECID are the spread to 1 row with _home and _clinic columns.
#  Each time we rename the added SPECID from the meta data according to
#  whose ENROLLID we are using to join.
trans_long<-trans_long %>%left_join(filter(meta,snv_qualified==T),
                                    by=c("Donor_ENROLLID"="ENROLLID"))%>% 
  select(HOUSE_ID=HOUSE_ID.x,Donor_ENROLLID,
         Recipient_ENROLLID,
         onset1,onset2,
         transmission,Donor_SPECID=SPECID) %>%
  mutate(Donor_sample=if_else(condition = grepl("HS",Donor_SPECID),
                              true = "Donor_home",
                              false = "Donor_clinic")) %>%
  spread(key = Donor_sample,value = Donor_SPECID)


# Now we do the same for the recipient.
trans_long<- left_join(trans_long,filter(meta,snv_qualified==T),
                       by=c("Recipient_ENROLLID"="ENROLLID")) %>%
  select(HOUSE_ID=HOUSE_ID.x,Donor_ENROLLID,
         Recipient_ENROLLID,
         onset1,onset2,
         transmission,Donor_home,Donor_clinic,
         Recipient_SPECID=SPECID) %>% 
  mutate(Recipient_sample=if_else(condition = grepl("HS",Recipient_SPECID),
                                  true = "Recipient_home",
                                  false = "Recipient_clinic")) %>%
  spread(key = Recipient_sample,value = Recipient_SPECID)

# Now we really need a long format. We want Donor_ENROLLID Recipient_ENROLLID SPECID1 SPECID2 
# and key columns for each. Then we can use the get freqs command and fit all the combinations.
# Doing this as two commands (one for the donor, and one for the recipient) means we get 
# all combinations Donor clinic - recipeint home donor clinic recipient clinic ect. 

trans_long<- trans_long %>% gather(key = "Donor_location",value = "SPECID1",Donor_home,Donor_clinic) %>%
  gather(key = "Recipient",value = "SPECID2",Recipient_home,Recipient_clinic) %>%
  filter(!(is.na(SPECID1) | is.na(SPECID2))) %>% 
  dplyr::rename("ENROLLID1"="Donor_ENROLLID",
                "ENROLLID2" = "Recipient_ENROLLID")

trans_long$pair_id<-1:length(trans_long$HOUSE_ID)

# Get the frequencies

all_trans_freq<-plyr::adply(trans_long,1,function(x){
  get_freqs(c(x$SPECID1,x$SPECID2),qual)},
  .parallel = T)

# Reduce to sites that are polymorphic in the donor.
all_trans_freq.comp<-polish_freq(all_trans_freq,freq1,0.02)
all_trans_freq.comp$found=all_trans_freq.comp$freq2>0.02 # was it found in the second sample

# Add gc_ul
all_trans_freq.comp <-mutate(all_trans_freq.comp,gc_ul1 = meta$gc_ul[match(SPECID1,meta$SPECID)],
       gc_ul2 = meta$gc_ul[match(SPECID2,meta$SPECID)])
# add collection dates
all_trans_freq.comp <- mutate(all_trans_freq.comp,collect1 = meta$collect[match(SPECID1,meta$SPECID)],
       collect2 = meta$collect[match(SPECID2,meta$SPECID)])

# Fit beta binomial model


max_nb<-200
beta_Nb_all<-trans_fit(subset(all_trans_freq.comp,freq1<0.5),
                   Nb_max=max_nb,model="BetaBin",
                   threshold=0.02,acc=accuracy_stringent,
                   pair_id)


beta_Nb_sum_all<- beta_Nb_all %>% group_by(pair_id) %>% do(straight_sum(.,max_nb))
beta_Nb_sum_all<-beta_Nb_sum_all[order(beta_Nb_sum_all$Nb,decreasing = T),]



beta_t_all<-beta_Nb_sum_all %>% mutate(
  CI = paste0(lower_95,"-",upper_95),
  lambda = paste0(Nb," (",CI,")")
) %>%
  select(pair_id,Nb,pair_id,CI) 
beta_t_all<-left_join(beta_t_all,
                  select(all_trans_freq.comp,pair_id,SPECID1,SPECID2,ENROLLID1,ENROLLID2,
                         collect1,collect2,transmission))%>%
  distinct(pair_id, .keep_all = TRUE) %>%
  dplyr::rename(donor_sample = collect1, recipient_sample=collect2,
                estimated_transmission_date=transmission) 

snv_data<- all_trans_freq.comp %>% group_by(pair_id,pcr_result) %>%
  summarize(minority_isnv=length(which(freq1<0.5)),
            transmitted_minority_isnv =length(which(freq1<0.5&found==T)))
beta_t_all<-left_join(beta_t_all,snv_data)

beta_t_all<-beta_t_all[order(c(beta_t_all$ENROLLID1, beta_t_all$ENROLLID2)),]%>%
  filter(!is.na(pair_id))


# Adding the ages.
ages<- read_csv("./data/reference/HIVE_ages_by_season.csv",
                col_types = cols(STUDY_ID=col_character()))%>%
  mutate(DOB = parse_date_time(DOB,c("db!y","%m/%d/%Y"))) %>%
  mutate(DOB=if_else(condition = (as.numeric(as.POSIXct(today())-DOB)/365)<AGEYR,
                     true = DOB-years(100),false=DOB))

beta_t_all<-left_join(beta_t_all,select(ages,STUDY_ID,DOB),by=c("ENROLLID1"="STUDY_ID")) %>%
  mutate(Donor_age=
           as.numeric(as.POSIXct(estimated_transmission_date)-DOB)/365.2425) %>%
  select(-DOB) %>%
  left_join(.,select(ages,STUDY_ID,DOB),by=c("ENROLLID2"="STUDY_ID")) %>%
  mutate(Recipient_age=
           as.numeric(as.POSIXct(estimated_transmission_date)-DOB)/365.2425) %>%
  select(-DOB)
out_beta_t_all<-beta_t_all %>% ungroup() %>%
  select(Nb,CI,Subtype=pcr_result,donor_sample,recipient_sample,estimated_transmission_date,
         Donor_age,Recipient_age,minority_isnv,transmitted_minority_isnv,ENROLLID1,ENROLLID2,
         SPECID1,SPECID2)
out_beta_t_all$Nb[is.na(out_beta_t_all$Nb)]<-">200"

# mark which samples were used in the analysis
out_beta_t_all<- out_beta_t_all %>% mutate(Used = ifelse(paste(SPECID1,SPECID2) %in% paste(out_beta_t$SPECID1,out_beta_t$SPECID2),
                                                 "X","-"))
write.csv(out_beta_t_all,"./data/processed/secondary/beta_bottlenecks_by_pair_all_combinations.csv")


# Test used samples match the first etimate
original <-out_beta_t
used <- select(filter(out_beta_t_all,Used=="X"),-Used)

# for some reason this sorts add NA
original <- original[order(c(original$SPECID1,original$SPECID2)),] %>% filter(!is.na(SPECID1))
used <- used[order(c(used$SPECID1,used$SPECID2)),] %>% filter(!is.na(SPECID1))

dim(anti_join(used,original))==c(0,14)

