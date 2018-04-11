require(tidyverse)
require(HIVEr)
cbPalette<-wesanderson::wes_palette("Zissou") # Set up figures
require(cowplot)
require(extrafont)
require(VGAM)
require(lubridate)
require(bbmle)

# Data files
trans_freq <- read_csv("./data/processed/secondary/trans_freq.csv")
trans_freq.comp<-read_csv("./data/processed/secondary/transmission_pairs_freq.poly.donor.csv")
accuracy_stringent<-read.csv("./data/reference/accuracy_stringent.csv",stringsAsFactors = F)
meta<-read_csv("./data/reference/all_meta.sequence_success.csv")


# Add gc_ul
trans_freq.comp<-mutate(trans_freq.comp,gc_ul1 = meta$gc_ul[match(SPECID1,meta$SPECID)],
                        gc_ul2 = meta$gc_ul[match(SPECID2,meta$SPECID)])

beta_total_fit<-trans_fit(subset(trans_freq.comp,freq1<0.5),
                          Nb_max=1000,model="BetaBin",
                          threshold=0.02,acc=accuracy_stringent,
                          pair_id)
write.csv(x = beta_total_fit,
          "./data/processed/secondary/beta_LL_nb1000.csv")


counts<-trans_freq.comp %>% group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1>0 & freq1<0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants,
         weight_factor = weight_factor_kk)#1)


# log normal
# This function discretizes the lognormal such that we
# only draw integers as in our bottleneck. This smooths
# the likelihood surface and means we must interpret the distribution
# accordingly.
dnormalize_lognormal<-function(x,meanlog,sdlog,maxnb){
  # Test we are going far enough out
#  stopifnot(qlnorm(0.999,meanlog,sdlog)<maxnb)
  pdf<-sum(dlnorm(1:maxnb,meanlog = meanlog,sdlog=sdlog))
  dlnorm(x,meanlog = meanlog,sdlog=sdlog)/pdf
}



lgnorm_fit<-dist_prob_wrapper("dnormalize_lognormal","meanlog,sdlog,maxnb")

# Plot to get an idea 
meanlog = seq(-5,5,0.5)
sdlog = seq(0.1,5,0.5)
nmeanlog = rep(meanlog,each = length(sdlog))
nsdlog = rep(sdlog,times = length(meanlog))
x<-purrr::map2(nmeanlog,nsdlog,.f = function(x,y) lgnorm_fit(data = beta_total_fit,weight = counts,meanlog = x,sdlog = y,maxnb=500))

nl_output<-tibble(meanlog = nmeanlog,nsdlog = nsdlog,negLL = unlist(x))

ggplot(nl_output,aes(x = meanlog,y= nsdlog,z = -negLL))+
  #geom_tile()+
  geom_contour(binwidth = 2)+
  geom_point(data = filter(nl_output,negLL ==min(negLL,na.rm = T)))+ 
  ylab("Standard deviation")+
  xlab("Mean")

nlog_model_fit_1000<-bbmle::mle2(minuslogl = lgnorm_fit,start = list(meanlog=-3.50,sdlog=2.10),
                          data = list(data = beta_total_fit,
                                      weight = counts,
                                      maxnb=1000),
                       # method="L-BFGS-B",
                        #lower=c(meanlog=0, sdlog=0),
                          control=list(trace=TRUE,
                                       maxit=5000))

#nlog_model_fit
# Call:
#   bbmle::mle2(minuslogl = lgnorm_fit, start = list(meanlog = -3.5, 
#                                                    sdlog = 2.1), data = list(data = beta_total_fit, weight = counts, 
#                                                                              maxnb = 200), control = list(trace = TRUE, maxit = 5000))
# 
# Coefficients:
#   meanlog      sdlog 
# -38.276714   5.972561 
# 
# Log-likelihood: -37.68 
# Maximum likelihood estimation
# 
# Call:
#   bbmle::mle2(minuslogl = lgnorm_fit, start = list(meanlog = -3.5, 
#                                                    sdlog = 2.1), data = list(data = beta_total_fit, weight = counts, 
#                                                                              maxnb = 200), control = list(trace = TRUE, maxit = 5000))
# 
# Coefficients:
#   Estimate Std. Error z value  Pr(z)
# meanlog -38.2767   164.3652 -0.2329 0.8159
# sdlog     5.9726    12.3918  0.4820 0.6298

#-2 log L: 75.36913 




# final plot

# are we going far enough to the left.
# How does the mean 95 and 99th quantile change using
# these parameters.
mean_customlnorm<-function(maxnb,meanlog,sdlog){
  sum(dnormalize_lognormal(seq(1,maxnb,1),
                meanlog,sdlog,maxnb)*c(1:maxnb))
}

qcustomlnorm<-function(p,maxnb,meanlog,sdlog){
  stopifnot(p<1)
  densities<-dnormalize_lognormal(seq(1,maxnb,1), meanlog,sdlog,maxnb)
  cdf<-0
  i=0
  while(cdf<p){
    i=i+1
    cdf = cdf+densities[i]
  }
  return(i)
}

size<-tibble(maxnb = seq(100,20000,by=100)) %>% rowwise() %>%
  mutate(mean = mean_customlnorm(maxnb,nlog_model_fit_1000@coef[1],nlog_model_fit_1000@coef[2]),
         q95 = qcustomlnorm(0.95,maxnb,nlog_model_fit_1000@coef[1],nlog_model_fit_1000@coef[2]),
         q99 = qcustomlnorm(0.99,maxnb,nlog_model_fit_1000@coef[1],nlog_model_fit_1000@coef[2]),
         q97 = qcustomlnorm(0.975,maxnb,nlog_model_fit_1000@coef[1],nlog_model_fit_1000@coef[2])
  )

slong<-size %>% gather(key = "parameter",value = "value",mean,q95,q97,q99)

require(plotly)
interesting<-ggplot(slong,aes(x=maxnb,y = value,color=parameter))+geom_line()

ggplotly(interesting)
plot(seq(0,20,1),dnormalize_lognormal(seq(0,20,1),
                                      nlog_model_fit@coef[1],nlog_model_fit@coef[2],200),
     xlab = "Bottleneck size",
     ylab = "PMF",
     type = "s")
points(seq(0,20,1),dnormalize_lognormal(seq(0,20,1),
                                       nlog_model_fit_1000@coef[1],nlog_model_fit_1000@coef[2],1000),
     xlab = "Bottleneck size",
     ylab = "PMF",
     type = "s",
     col = "red")
