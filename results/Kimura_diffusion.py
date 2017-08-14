
# coding: utf-8

# In[1]:

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


# In[2]:

from ipywidgets import interact
import numpy as np
import pandas as pd
from bokeh.io import push_notebook, show, output_notebook
from bokeh.plotting import figure
import numpy as np
import mpmath as mp
import copy 
from matplotlib import pyplot as plt
from scipy.integrate import quad 
from scipy.optimize import minimize
import re

output_notebook()
get_ipython().magic(u'matplotlib inline')



# In[3]:

def write_to_summary(line_pattern,value):
    lines = []
    rline_pat = re.compile('^'+line_pattern+".*")
    with open("./results.table.tsv",'r') as results:
        for line in results:
            if rline_pat.match(line):
                line = line_pattern +"\t" + str(value)+"\n"
            lines.append(line)
            
    with open("./results.table.tsv",'w') as output:
        for l in lines:
            output.write(l)


# ## Reading in the accuracy table

# In[4]:

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

# In[5]:

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

# In[6]:

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


# In[7]:

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


# I am assuming one generation is 6 hours. And I am fitting for only samples that were taken at least one day appart.

# In[8]:

intra=pd.read_csv("./Intrahost_initially_present.csv")
#intra.loc[intra.within_host_time == 0, "within_host_time"] = 0.2 # assume about 5 hours a day passed between samples
#intra["generations"] = intra.within_host_time*24/10

intra_minor=intra.loc[intra.freq1<0.5]
intra_minor=intra_minor.loc[intra_minor.within_host_time>0]


intra_minor.count()


# In[9]:

print(intra_minor.head())
intra_minor.loc[intra_minor.donor_class=="Nonsynonymous",]


# In[10]:

LL = np.arange(1,50,1) # These are the effective population sizes
likes = []
for d in LL:
    print( 'working with: '+ str(d))
    likes.append(likelihood(d,intra_minor,6,acc_stringent))


# In[11]:

max_likes=[-1*x for x in likes] # convert back to positive log likelihood.
with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL[19:49],max_likes[19:49])
    ax.plot(LL[max_likes.index(max(max_likes))], max(max_likes), 'ro')
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Population")

Ne=LL[max_likes.index(max(max_likes))]
print(Ne)

write_to_summary("Diffusion model (6) Ne:",Ne)


# ### Save output for plots

# In[12]:

out_LL =[]
for i in range(0,len(LL)):
        out_LL.append({"Ne":LL[i],"LL":max_likes[i]})
out_pd = pd.DataFrame(out_LL)

out_pd.to_csv("./intrahost_change_in_freq.LL.csv")


# In[13]:

cutoff=max(max_likes)-1.92


# In[14]:

above_cut=[x for x in max_likes if x >cutoff ]


# The 95% confidence interval is 

# In[15]:

CI = [LL[max_likes.index(above_cut[0])],LL[max_likes.index(above_cut[-1])]] # get the first and last Ne sizes in the CI range

write_to_summary("Diffusion model (6) CI:",CI)


# In[16]:

LL12 = np.arange(1,50,1) # These are the effective population sizes
likes12 = []
for d in LL12:
    print( 'working with: '+ str(d))
    likes12.append(likelihood(d,intra_minor,12,acc_stringent))


# In[17]:

max_likes12=[-1*x for x in likes12] # convert back to positive log likelihood.
with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL[10:49],max_likes[10:49])
    ax.plot(LL[max_likes.index(max(max_likes))], max(max_likes), 'ro')
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Population")

Ne=LL[max_likes12.index(max(max_likes12))]
print(Ne)
write_to_summary("Diffusion model (12) Ne:",Ne)


# In[18]:

cutoff=max(max_likes12)-1.92

above_cut=[x for x in max_likes12 if x >cutoff ]


# In[19]:

CI =[LL[max_likes12.index(above_cut[0])],LL[max_likes12.index(above_cut[-1])]] # get the first and last Ne sizes in the CI range
write_to_summary("Diffusion model (12) CI:",CI)


# # Looking at Nonsynonymous and Synonymous mutations
# 
# ## Nonsynonymous first

# In[20]:

intra_nonsynon = intra_minor.loc[intra_minor.donor_class=="Nonsynonymous",]

LL = np.arange(1,50,1) # These are the effective population sizes
likes_nonsynom = []
for d in LL:
    print( 'working with: '+ str(d))
    likes_nonsynom.append(likelihood(d,intra_nonsynon,6,acc_stringent))


# In[21]:

max_likes_nonsynom=[-1*x for x in likes_nonsynom]
with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL[12:40],max_likes_nonsynom[12:40])
    ax.plot(LL[max_likes_nonsynom.index(max(max_likes_nonsynom))], max(max_likes_nonsynom), 'ro')
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Population")

Ne=LL[max_likes_nonsynom.index(max(max_likes_nonsynom))]
print(Ne)


# ## Synonymous

# In[22]:

intra_synon = intra_minor.loc[intra_minor.donor_class=="Synonymous",]

LL = np.arange(1,50,1) # These are the effective population sizes
likes_synom = []
for d in LL:
    print( 'working with: '+ str(d))
    likes_synom.append(likelihood(d,intra_synon,6,acc_stringent))


# In[23]:

max_likes_nsynon=[-1*x for x in likes_synom]
with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL[20:49],max_likes_nsynon[20:49])
    ax.plot(LL[max_likes_nsynon.index(max(max_likes_nsynon))], max(max_likes_nsynon), 'ro')
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Population")

Ne=LL[max_likes_nsynon.index(max(max_likes_nsynon))]
print(Ne)


# # Comparision of All, Nonsynonymous and Synonymous fits

# In[24]:

with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL[15:49],max_likes_nsynon[15:49])
    ax.plot(LL[max_likes_nsynon.index(max(max_likes_nsynon))], max(max_likes_nsynon), 'ro')
    
    ax.plot(LL[15:49],max_likes_nonsynom[15:49])
    ax.plot(LL[max_likes_nonsynom.index(max(max_likes_nonsynom))], max(max_likes_nonsynom), 'ro')
    
    ax.plot(LL[15:49],max_likes[15:49])
    ax.plot(LL[max_likes.index(max(max_likes))], max(max_likes), 'ro')
    
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Ne")



# # Looking at number of terms needed for infinite sums
# 
# Above we used 50. Do the outcomes change if we use 25 or 100
# 
# ## 25
# 

# In[25]:

def non_fixed(p,x,t,N,sensitivity=False,*args,**kwargs):
    gc_ul = kwargs.get('gc_ul', None) # If these  variables are supplied get them. They are used to account for sensitivity. sensivity must be set to True
    acc = kwargs.get('acc', None)
    ith=[ith_term(i=1,p=p,x=x,t=t,N=N),ith_term(i=2,p=p,x=x,t=t,N=N)] # The first 2 terms are done to begin
    i = 3
#    while (ith[i-3]-ith[i-2])>1e-5: # Continue until the difference between the last 2 terms is less than 1e-5
#    while ith[i-2]>1e-3: 
    while i<=25:
        ith.append(ith_term(i=i,p=p,x=x,t=t,N=N))
        i+=1
    #print(ith)
    perfect_detect = np.sum(ith)


    if sensitivity == False :
        return(perfect_detect)
    
    else:# this is the probability of the variant being found where it was given the sensitivity. 
            # This is not used in the fitting of the model. All of these variants are found. For each N that we try this term doesn't change. So it is a constant not dependent on N and so doesn't affect the estimate.
            # It is useful though in getting a pdf that sums to 1 in the plots. If we don't have perfect sensitivity for the lost variants then we should treat these the same.
        acc_gc=10**(np.floor(np.log10(gc_ul)))
        if acc_gc>1e5: # set for the max
            acc_gc=1e5


        if x<0.05 and x>0.02:
            sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==0.02),'sensitivity']
            sense = sense.iloc[0]
            prob_detect = perfect_detect*sense
        elif x<0.1 and x>0.05:
            sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==0.05),'sensitivity']
            sense = sense.iloc[0]
            prob_detect = perfect_detect*sense
        else :
            prob_detect = perfect_detect

        return(prob_detect) 
def boundaries(p,t,N,final,gc_ul=1e5,sensitivity=False,*args,**kwargs):
    acc = kwargs.get('acc', None)
    
    #if final !=0 or final !=1:
    #    raise(ValueError,"Please select 0 or 1 as final frequency")
    if final==0:
        fixed_freq=1-p   # set for loss. The probabilty the other allele is fixed
    elif final ==1:
        fixed_freq = p

    ith=[ith_term_fixed(i=1,p=fixed_freq,t=t,N=N),ith_term_fixed(i=2,p=fixed_freq,t=t,N=N)] # The first 2 terms are done to begin 
    i = 3
#    while (ith[i-3]-ith[i-2])>1e-5: # Continue until the difference between the last 2 terms is less than 1e-5
#    while ith[i-2]>1e-3: 
    while i<=25:
        ith.append(ith_term_fixed(i=i,p=fixed_freq,t=t,N=N))
        i+=1
    fixed = fixed_freq+np.sum(ith)
    #print(fixed)
    if sensitivity == False:
            return(fixed)

        
    elif sensitivity == True:
        if final ==0 : 
            lost_p = p
        elif final==1:
            lost_p = 1-p
        below_threshold = below_cut(p=lost_p,t=t,N=N)
        missed= just_missed(p=lost_p,t=t,N=N,gc_ul=gc_ul,acc=acc)
        lost = below_threshold+missed
        return(lost+fixed)


# In[26]:

LL = np.arange(1,50,1)
likes_25 = []
for d in LL:
    print( 'working with: '+ str(d))
    likes_25.append(likelihood(d,intra_minor,6,acc_stringent))


# In[27]:

max_likes_25=[-1*x for x in likes_25]
with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL,max_likes_25)
    ax.plot(LL[max_likes_25.index(max(max_likes_25))], max(max_likes_25), 'ro')
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Population")

Ne=LL[max_likes_25.index(max(max_likes_25))]
print(Ne)


# In[28]:

with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(max_likes,max_likes_25)
    ax.set_ylabel("25")
    ax.set_xlabel("50")


# In[29]:

[max_likes[i]-max_likes_25[i] for i in range(0,len(max_likes))]


# ## 100

# In[30]:

def non_fixed(p,x,t,N,sensitivity=False,*args,**kwargs):
    gc_ul = kwargs.get('gc_ul', None) # If these  variables are supplied get them. They are used to account for sensitivity. sensivity must be set to True
    acc = kwargs.get('acc', None)
    ith=[ith_term(i=1,p=p,x=x,t=t,N=N),ith_term(i=2,p=p,x=x,t=t,N=N)] # The first 2 terms are done to begin
    i = 3
#    while (ith[i-3]-ith[i-2])>1e-5: # Continue until the difference between the last 2 terms is less than 1e-5
#    while ith[i-2]>1e-3: 
    while i<=100:
        ith.append(ith_term(i=i,p=p,x=x,t=t,N=N))
        i+=1
    #print(ith)
    perfect_detect = np.sum(ith)


    if sensitivity == False :
        return(perfect_detect)
    
    else:# this is the probability of the variant being found where it was given the sensitivity. 
            # This is not used in the fitting of the model. All of these variants are found. For each N that we try this term doesn't change. So it is a constant not dependent on N and so doesn't affect the estimate.
            # It is useful though in getting a pdf that sums to 1 in the plots. If we don't have perfect sensitivity for the lost variants then we should treat these the same.
        acc_gc=10**(np.floor(np.log10(gc_ul)))
        if acc_gc>1e5: # set for the max
            acc_gc=1e5


        if x<0.05 and x>0.02:
            sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==0.02),'sensitivity']
            sense = sense.iloc[0]
            prob_detect = perfect_detect*sense
        elif x<0.1 and x>0.05:
            sense=acc.loc[(acc.gc_ul==acc_gc) & (acc.freq==0.05),'sensitivity']
            sense = sense.iloc[0]
            prob_detect = perfect_detect*sense
        else :
            prob_detect = perfect_detect

        return(prob_detect) 
def boundaries(p,t,N,final,gc_ul=1e5,sensitivity=False,*args,**kwargs):
    acc = kwargs.get('acc', None)
    
    #if final !=0 or final !=1:
    #    raise(ValueError,"Please select 0 or 1 as final frequency")
    if final==0:
        fixed_freq=1-p   # set for loss. The probabilty the other allele is fixed
    elif final ==1:
        fixed_freq = p

    ith=[ith_term_fixed(i=1,p=fixed_freq,t=t,N=N),ith_term_fixed(i=2,p=fixed_freq,t=t,N=N)] # The first 2 terms are done to begin 
    i = 3
#    while (ith[i-3]-ith[i-2])>1e-5: # Continue until the difference between the last 2 terms is less than 1e-5
#    while ith[i-2]>1e-3: 
    while i<=100:
        ith.append(ith_term_fixed(i=i,p=fixed_freq,t=t,N=N))
        i+=1
    fixed = fixed_freq+np.sum(ith)
    #print(fixed)
    if sensitivity == False:
            return(fixed)

        
    elif sensitivity == True:
        if final ==0 : 
            lost_p = p
        elif final==1:
            lost_p = 1-p
        below_threshold = below_cut(p=lost_p,t=t,N=N)
        missed= just_missed(p=lost_p,t=t,N=N,gc_ul=gc_ul,acc=acc)
        lost = below_threshold+missed
        return(lost+fixed)


# In[31]:

LL = np.arange(1,50,1)
likes_100 = []
for d in LL:
    print( 'working with: '+ str(d))
    likes_100.append(likelihood(d,intra_minor,6,acc_stringent))


# In[32]:

max_likes_100=[-1*x for x in likes_100]
with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(LL,max_likes_100)
    ax.plot(LL[max_likes_100.index(max(max_likes_100))], max(max_likes_100), 'ro')
    ax.set_ylabel("Log Likelihood")
    ax.set_xlabel("Population")

Ne=LL[max_likes_100.index(max(max_likes_100))]
print(Ne)


# In[33]:

with plt.style.context('fivethirtyeight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(max_likes,max_likes_100)
    ax.set_ylabel("100")
    ax.set_xlabel("50")


# In[34]:

[max_likes[i]-max_likes_100[i] for i in range(0,len(max_likes))]


# In[35]:

def min_diff(data_list):
    out_list = [np.abs(data_list[i]-data_list[i+1]) for i in range(0,(len(data_list)-1))]
    return(min(out_list))
                                                           


# In[36]:

print("25:%f\n50:%f\n100:%f" %(min_diff(max_likes_25),min_diff(max_likes),min_diff(max_likes_100)))


# In[37]:

x=boundaries(p=0.5,t=100,N=30,final=1,sensitivity = True,gc_ul=1e5,acc=acc_stringent)
x


# 
