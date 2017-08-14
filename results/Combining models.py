
# coding: utf-8

# In[11]:

from IPython.display import HTML
import warnings
#warnings.filterwarnings('ignore')
HTML('''<script>
  function code_toggle() {
    if (code_shown){
      $('div.input').hide('500');
      $('#toggleButton').val('Show Code')
    } else {
      $('div.input').show('500');
      $('#toggleButton').val('Hide Code')
    }
    code_shown = !code_shown
  }

  $( document ).ready(function(){
    code_shown=false;
    $('div.input').hide()
  });
</script>
<form action="javascript:code_toggle()"><input type="submit" id="toggleButton" value="Show Code"></form>''')


# In[1]:

from ipywidgets import interact
import numpy as np
import pandas as pd
from bokeh.io import push_notebook, show, output_notebook,hplot
from bokeh.plotting import figure
import numpy as np
import mpmath as mp
import copy 
from matplotlib import pyplot as plt
from scipy.integrate import quad 
from scipy.optimize import minimize
from bokeh.charts import Line
from bokeh.models.layouts import Row
output_notebook()
get_ipython().magic(u'matplotlib inline')



# ## Reading in the accuracy table

# In[2]:

acc_stringent = pd.read_csv("../data/reference/accuracy_stringent.csv")

acc_stringent


#  ## Equation set up
#  
#  Let p(p,x,t,N) be the time dependent probability of a variant at x after t generations when the initial frequency was p, and the effective population size is N.
#  
#  From Kimura 1955 we have,
#  
#  $$
# p(p,x,t,N) = \sum_{i=1}^{\infty}pqi(i+1)(2i+1)F(1-i,i+2,2,p) \times F(1-i,i+2,2,x) e^{-[i(i+1)/2N]t}
# $$
# 
# Where $q=1-p$ and $F$ is the hypergeometric function.
#  
# The code is below.

# In[3]:

def ith_term(i,p,t,x,N): # proofed JT 5/22/17
    q=1-p
    first = p*q*i*(i+1)*((2*i)+1)    
    geometric_1= mp.hyp2f1(1-i,i+2,2,p,zeroprec=10) # verify this the correct function (it's the only hypergeometric I could find with 4 variables) - This is correct 
    geometric_2= mp.hyp2f1(1-i,i+2,2,x,zeroprec = 10)
    exponent= i*(i+1)*t/(2*N) # 4N in the book
    out=first*geometric_1*geometric_2*np.exp(-1*exponent)
    return(float(out))

def non_fixed(p,x,t,N,sensitivity=False,*args,**kwargs):
    gc_ul = kwargs.get('gc_ul', None) # If these  variables are supplied get them. They are used to account for sensitivity. sensivity must be set to True
    acc = kwargs.get('acc', None)
    ith=[ith_term(i=1,p=p,x=x,t=t,N=N),ith_term(i=2,p=p,x=x,t=t,N=N)] # The first 2 terms are done to begin
    i = 3
#    while (ith[i-3]-ith[i-2])>1e-5: # Continue until the difference between the last 2 terms is less than 1e-5
#    while ith[i-2]>1e-3: 
    while i<=50:
        ith.append(ith_term(i=i,p=p,x=x,t=t,N=N))
        i+=1
    #print(ith)
    perfect_detect = np.sum(ith)


    if sensitivity == False :
        return(perfect_detect)
    
    else:# this is the probability of the variant being found where it was given the sensitivity. 
            # This is not used in the fitting of the model. All of these variants are found. For each N that we try this term doesn't change. So it is a constant not dependent on N and so doesn't affect the estimate.
            # It is useful though in getting a pdf that sums to 1 in the plots. If we don't have perfect sensitivity for the lost variants then we should treat these the same.
        acc_gc=10**(np.floor(np.log10(gc_ul))) # round down to the nearest log10 
        if acc_gc>1e5: # set for the max. We assume we do not gain sensitivity above this cut off. We probably do, but we don't have data on it so this is more conservative
            acc_gc=1e5

        ## Here we assume the accuracy of each range is the same as the smaller range
        if x<0.05 and x>0.02:
            sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==0.02),'sensitivity']
            sense = sense.iloc[0]
            prob_detect = perfect_detect*sense
        elif x<0.1 and x>0.05:
            sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==0.05),'sensitivity']
            sense = sense.iloc[0]
            prob_detect = perfect_detect*sense
        else :
            prob_detect = perfect_detect # We assume perfect detection above 10%

        return(prob_detect) 


# The time depenent probabilty of a variant being lost  or not observed at generation t is given by the sum of the probabity that it is lost by generation t p(p,0,t,N) and the probability that it is not detected due to the limit of detection or low sensitivity to detect rare variants P(below_cut)+P(p,x,t,N)\*P(not_detected,x). 
# 
# p(not_detect) = P(other_allele_fixed)+p(below_threshold)+p(present_but_overlooked)
# 
# 
# Variants that are fixed (there are non in the data set but its good to handel anyway)
# 
# p(other_allele_not_detected) = p(fixed)+p(other_below_theshold)+p(other_present_but_overlooked)
# 
#  
# $$
# p(p,0,t,N) = q +\sum_{i=1}^{\infty}(2i+1)pq(-1)^i F(1-i,i+2,2,q) e^{-[i(i+1)/2N]t}
# $$
#  
# Where q is defined as above. (Note : this is simply the probability of fixation for a variant at initial frequency q.
# $$
# P(\text{below_cut}) = \int_0^{0.02} p(p,x,t,N)dx
# $$
#  
# Note the limit of detection is 0.02.
#  
# $$
# P(\text{present_but_not_detected}) = \sum_{f_e}^{[0.02,0.05,0.10)} \big(\text{FNR}|\text{Titer}_r,f_e) \int_{f_e}^{f_e+1} p(p,x,t,N)dx
# $$
#  
# Where $(\text{FNR}|\text{Titer}_r,f_e)$ is the false negative rate given the frequency and the sample titer.

# In[4]:

def ith_term_fixed(i,p,t,N):# proofed JT 5/22/17
    first = (2*i+1)*p*(1-p)*(-1)**i
    geometric = mp.hyp2f1(1-i,i+2,2,p,zeroprec=10)
    exponent= i*(i+1)*t/(2*N) # 4N in the book
    out = first*geometric*np.exp(-1*exponent)
    return(float(out))

def below_cut(p,t,N):
    return quad(lambda x : non_fixed(p,x,t,N),0,0.02)[0] # proofed JT 5/22/17


def just_missed(p,t,N,gc_ul,acc): # This accounts for the variants that are present but we don't detect them 
    acc_gc=10**(np.floor(np.log10(gc_ul))) # again round down to the nearest log10
    if acc_gc>1e5: # set for the max
        acc_gc=1e5
    uncert_term=[]
    f=[0.02,0.05,0.10]
    for i in range(0,(len(f)-1)):
        sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==f[i]),'sensitivity']
        uncert=1-sense.iloc[0]
        #print(uncert)
    # The prob the variant is missed because it is between f[i] and f[i+1] given the sample size 
        uncert_term.append(quad(lambda x : non_fixed(p,x,t,N),f[i],f[i+1])[0]*uncert)
    #print(uncert_term)
    return(np.sum(uncert_term))

def boundaries(p,t,N,final,gc_ul=1e5,sensitivity=False,*args,**kwargs):
    acc = kwargs.get('acc', None)
    
    #if final !=0 or final !=1:
    #    raise(ValueError,"Please select 0 or 1 as final frequency")
    if final==0:
        fixed_freq=1-p   # set for loss. The probabilty the other allele is fixed
    elif final ==1:
        fixed_freq = p    # In this case this is the frequency of the allele we want to fix

    ith=[ith_term_fixed(i=1,p=fixed_freq,t=t,N=N),ith_term_fixed(i=2,p=fixed_freq,t=t,N=N)] # The first 2 terms are done to begin 
    i = 3
#    while (ith[i-3]-ith[i-2])>1e-5: # Continue until the difference between the last 2 terms is less than 1e-5
#    while ith[i-2]>1e-3: 
    while i<=50:
        ith.append(ith_term_fixed(i=i,p=fixed_freq,t=t,N=N))
        i+=1
    fixed = fixed_freq+np.sum(ith) # from the equation above
    #print(fixed)
    if sensitivity == False:
            return(fixed)

        
    elif sensitivity == True:
        if final ==0 : 
            lost_p = p # this is the frequency of the variant we want to observe loss 
        elif final==1:
            lost_p = 1-p # in this case we want the loss of the other allele
        below_threshold = below_cut(p=lost_p,t=t,N=N)
        missed= just_missed(p=lost_p,t=t,N=N,gc_ul=gc_ul,acc=acc)
        lost = below_threshold+missed
        return(lost+fixed)


# In[5]:

def pdf(p,x,t,N,gc_ul=1e5,sensitivity = False, acc=acc_stringent):
    if x <1 and x>0 :
        return(non_fixed(p=p,x=x,t=t,N=N,sensitivity=sensitivity,gc_ul=gc_ul,acc=acc))
    else:
        return(boundaries(p=p,final=x,N=N,t=t,sensitivity=sensitivity,gc_ul=gc_ul,acc=acc))
               


def likelihood(n,data,generation,acc):
    """
    This function takes in a diffusion rate parameter and outputs the negative log likelihood
    of the parameter.
    generation is generation time in hours.
    """
    local_intra=copy.deepcopy(data)
    
    local_intra["generations"] = local_intra.within_host_time*24/generation # convert days to generations
    
    local_intra["log_like"] = local_intra.apply(lambda row: -1*np.log(pdf(p=row["freq1"],x=row["freq2"],N=n,t=row["generations"],sensitivity = True,gc_ul=row["gc_ul.2"],acc=acc)), axis=1)
    #return(local_intra)
    return local_intra.log_like.sum()


# In[6]:

def not_detected(mu,Ne,t,gc_ul,acc=acc_stringent):
 # p0
   naught = 1-quad(lambda x :(2*mu*Ne/x) * np.exp(-1*2*Ne*x/t),0.001,1)[0]  # probability it is at 0
 # below cut
   below = 1-quad(lambda x :(2*mu*Ne/x) * np.exp(-1*2*Ne*x/t),0.02,1)[0]  # below detection
 
 # missed

   acc_gc=10**(np.floor(log10(gc_ul))) # This will round the gc down to the nearest log as we have always done to be conservative
   uncert_term=[]
   f=[0.02,0.05,0.10]
   for i in range(0,(len(f)-1)):
       sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==f[i]),'sensitivity']
       uncert=1-sense.iloc[0]
       #print(uncert)
   # The prob the variant is missed because it is between f[i] and f[i+1] given the sample size 
       integrand = quad(lambda x :(2*mu*Ne/x) * np.exp(-1*2*Ne*x/t) ,f[i],f[i+1])[0]
       uncert_term[i]=integrand*uncert # the probability it is present * the probability of not seeing it
 
   missed = np.sum(uncert_term)
 
   return(naught+below+missed) # the probability of being present and missed or not present

def g_ft(x,mu,Ne,t,gc_ul,sensitivity,acc=acc_stringent): # x is the frequency anything less than 0.001 is 0
   if x>0:
       return (2*mu*Ne/x) * np.exp(-1*2*Ne*x/t) 
   elif x==0 and sensitivity==True:
       return not_detected(mu,Ne,t,gc_ul,acc)
   elif x == 0  and sensitivity==False:
       return 1-quad(lambda x :(2*mu*Ne/x) * np.exp(-1*2*Ne*x/t),0.001,1)[0]  # probability it is at 0


# In[7]:

def conditional_path(xobs,xexp,tobs,texp,mu,Ne,gc_ul,acc): # we assume x0 = 0 at t0
    #This is the probability of observing a variant at frequency x1 at t1 given it was found at x2 at t2 and 0 at t0
    # xobs is subject to our sensitivity limits. xexp is not.
    drift_t = np.abs(tobs-texp)
    if tobs > texp : # then we are interested in the path from 0 to xobs.
        if xexp > 0: 
            return g_ft(xexp,mu,Ne,texp,gc_ul,False,acc)*pdf(xexp,xobs,drift_t,Ne,gc_ul,True,acc) # Sensitivity analysis since x2 is subject to our observed limitations
        elif xexp ==0 :
            return g_ft(xexp,mu,Ne,texp,gc_ul,False,acc)*g_ft(xobs,mu,Ne,drift_t,gc_ul,True,acc)
    elif tobs<texp : 
        return pdf(xobs,xexp,drift_t,Ne,gc_ul,False,acc)


# Let's look at the pdf for frequency at 4 generations (1 day) when we observe a variant at 5% on day 2. (8 generations)

# In[8]:

intermediate_freqs = np.arange(0,1.01,0.01)
prob = [conditional_path(0.50,xexp,12,9,4e-6,31,1e5,acc_stringent) for xexp in intermediate_freqs]


# In[23]:

def normalize(l):
    denom = max(l)
    return [x/denom for x in l]


alpha_list = np.array(normalize(prob))

path10 = [ [0,x] for x in intermediate_freqs]
time10 = [ [0,9] for x in intermediate_freqs]

path20 =[ [x,0.5] for x in intermediate_freqs]
time20 = [ [9,12] for x in intermediate_freqs]

p = figure(plot_width=400, plot_height=400)

first = p.multi_line(time10,path10,alpha = 0.1,line_width = 1)#,line_alpha = alpha_list)
second = p.multi_line(time20,path20,alpha = 0.1,line_width = 1)#,line_alpha = alpha_list)

mid = p.circle(9,intermediate_freqs, size=alpha_list*5, alpha=alpha_list)

l = figure(plot_width=400, plot_height=400)

dist = l.line(intermediate_freqs,prob)

both  = hplot(p,l)


# In[19]:

def update(xobs=0.5,tobs=10,test=9):
    # Calculate data
    intermediate_freqs = np.arange(0,1.01,0.01)
    prob = [conditional_path(xobs,xest,tobs,test,4e-6,31,1e5,acc_stringent) for xest in intermediate_freqs]
    # plot lines
    alpha_list = np.array(normalize(prob))

    path1 = [ [0,x] for x in intermediate_freqs]
    time1 = [ [0,test] for x in intermediate_freqs]

    path2 =[ [x,xobs] for x in intermediate_freqs]
    time2 = [ [test,tobs] for x in intermediate_freqs]

    first.data_source.data["x"] = time1
    first.data_source.data["y"] = path1
    second.data_source.data["x"] = time2
    second.data_source.data["y"] = path2
    

    mid.data_source.data['fill_alpha']=alpha_list
    mid.data_source.data["size"] = alpha_list*5
    mid.data_source.data["x"] = test*np.ones(len(intermediate_freqs))
    
    dist.data_source.data["y"] = prob
    push_notebook()


# In[24]:

show(both,notebook_handle=True)


# In[25]:

interact(update,xobs=(0,1,0.1),tobs = (0,28),test = (0,28))


# In[ ]:



