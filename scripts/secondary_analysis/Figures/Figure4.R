require(tidyverse)
require(extrafont)
require(cowplot)
require(HIVEr)
cbPalette<-wesanderson::wes_palette("Zissou") # Set up figures
# --------------------------------- Functions ---------------------------------
# PA probability
Pt_PA<-function(x,l,max_Nb){
  s<-0
  for(i in 1:max_Nb){
    c<- ((1-x)^i) *  (l^i)/((exp(l)-1)*factorial(i) )
    s<-s+c
  }
  return(1-s)
}
# --------------------------------- Data Files --------------------------------
#   Read in the csv files used in the data analysis below
# -----------------------------------------------------------------------------

no_cut.trans_freq<-read_csv("./data/processed/secondary/no_cut_trans_freq.csv", 
                            col_types = list(
                              ENROLLID1= col_character(),
                              ENROLLID2= col_character(),
                              SPECID1 = col_character(),
                              SPECID2 = col_character(),
                              pair_id = col_double()))
no_cut.trans_freq.comp<-read_csv("./data/processed/secondary/no_cut_transmission_pairs_freq.poly.donor.csv")
# --------------------------------- Figure 4A ---------------------------------
#   The frequency of iSNV between transmission pairs with the original pipeline
# -----------------------------------------------------------------------------

no_cut_no_infer<- no_cut.trans_freq %>% group_by(SPECID1,SPECID2,pair_id,chr,pos) %>%
  mutate(minor_infer = ifelse( 
    # There is one site where the minor allele is not the reference. Then it goes away. 
    # This gives a length of 0 and the minor is not infered so it's False. Otherwise
    # we set the minor_infered column to the output of the comparison.
    length(which(ref==var)==which(freq1==min(freq1)))==0,
    F,
    (which(ref==var)==which(freq1==min(freq1))&min(freq1)<0.1)
  )) %>%
  filter(minor_infer==F)


nc_ni_tv_plot<-ggplot(no_cut_no_infer,aes(freq1,freq2))+geom_point()
nc_ni_tv_plot

save_plot("./results/Figures/Figure4A.pdf", nc_ni_tv_plot,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure4A.pdf")

write.csv(x = select(no_cut_no_infer,SPECID1,SPECID2,freq1,freq2),
          "./results/Figures/data/Figure4A.csv")

# --------------------------------- Figure 4B ---------------------------------
#   The transmission bottleneck in the absence of a frequency threshold.
# -----------------------------------------------------------------------------

# Fit the model
no_cut_no_infer.comp<- no_cut.trans_freq.comp %>% group_by(SPECID1,SPECID2,pair_id,chr,pos) %>%
mutate(minor_infer = ifelse( 
  # There is one site where the minor allele is not the reference. Then it goes away. 
  # This gives a length of 0 and the minor is not infered so it's False. Otherwise
  # we set the minor_infered column to the output of the comparison.
  length(which(ref==var)==which(freq1==min(freq1)))==0,
  F,
  (which(ref==var)==which(freq1==min(freq1))&min(freq1)<0.1)
)) %>%
  filter(minor_infer==F)

nc_ni.pa_total_fit<-suppressWarnings(trans_fit(no_cut_no_infer.comp,l=seq(0.01,10,0.01),
                                               Nb_max=100,model="PA",
                                               threshold=NULL,acc=NULL))
nc_ni.pa_fit_sum<-model_summary(nc_ni.pa_total_fit)
print(nc_ni.pa_fit_sum)


# Plot the fit
# Sliding window with window of w step of step
w<-0.05
step = 0.025
windows<-tibble(s=seq(0.02,1-w,by=step),end=s+w)

out <- windows %>% rowwise() %>%
  mutate(iSNV = nrow(filter(no_cut_no_infer.comp,freq1>=s,freq1<end)),
         transmitted = nrow(filter(no_cut_no_infer.comp,freq1>=s,freq1<end,found==T)),
         freq = mean(no_cut_no_infer.comp$freq1[which(no_cut_no_infer.comp$freq1>=s & no_cut_no_infer.comp$freq1<end)]),
         prob = transmitted/iSNV,
         error_bottom = qbinom(c(0.025),iSNV,prob)/iSNV,
         error_top = qbinom(c(0.975),iSNV,prob)/iSNV,
         many = iSNV>5)

model<-tibble(s = seq(0,1,0.01))
model<- mutate(model,prob = Pt_PA(s,nc_ni.pa_fit_sum$lambda,100),
               lower = Pt_PA(s,nc_ni.pa_fit_sum$lower_95,100),
               upper = Pt_PA(s,nc_ni.pa_fit_sum$upper_95,100))

window_data.p<-ggplot()+geom_point(data = out,aes(x=freq,y=prob,alpha=many))+
  geom_errorbar(data=out,aes(x=freq,ymin=error_bottom,ymax=error_top,alpha=many))+
  geom_line(data=model,aes(x=s,y=prob),color=cbPalette[5])+
  geom_ribbon(data=model,aes(x=s,ymin=lower,ymax=upper),alpha=0.5,fill=cbPalette[5])+
  scale_alpha_manual(values=c(0,1))+theme(legend.position = 'none')+
  xlab("Frequency in Donor")+ylab("Probability of transmission")+
  geom_point(data=trans_freq.comp,aes(x=freq1,y=as.numeric(found)+(as.numeric(found)-0.5)/10),alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25))



save_plot("./results/Figures/Figure4B.pdf", window_data.p,
          base_aspect_ratio = 1.3)
embed_fonts("./results/Figures/Figure4B.pdf")

write.csv(x = select(no_cut_no_infer.comp,SPECID1,SPECID2,freq1,freq2),
          "./results/Figures/data/Figure4B_points.csv")
write.csv(x = model,
          "./results/Figures/data/Figure4B_model.csv")


