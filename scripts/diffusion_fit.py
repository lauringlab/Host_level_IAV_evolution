import argparse
import numpy as np
import pandas as pd
import mpmath as mp
import copy
from scipy.integrate import quad
import os
import click
from joblib import Parallel, delayed
import multiprocessing


parser = argparse.ArgumentParser(description='This script runs the Kirmura diffusion maximum likelihood Ne fitting on all files in a directory. This is being used to ensure the method is not biased. This was developed in the kimura notebook',usage ="")

parser.add_argument('inputDir', metavar='inputDir', nargs='+',
                    help='The directory containing the simulated data')
parser.add_argument('outCsv', metavar='outCsv', nargs='+',
                    help='The output CSV')

acc_stringent = pd.read_csv("../data/reference/accuracy_stringent.csv")

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

    local_intra["log_like"] = local_intra.apply(lambda row: -1*np.log(pdf(p=row["freq1"],x=row["sim"],N=n,t=row["generations"],sensitivity = True,gc_ul=row["gc_ul_2"],acc=acc)), axis=1)
    #return(local_intra)
    return local_intra.log_like.sum()

def ML_fit(data,acc_stringent):
    intra=pd.read_csv(data)

    intra_minor=intra.loc[intra.freq1<0.5]
    intra_minor=intra_minor.loc[intra_minor.within_host_time>0]

    LL = np.arange(1,80,1) # These are the effective population sizes
    likes = []
    with click.progressbar(range(1,(max(LL)+1))) as bar:
        for i in bar:
    #for d in LL:
            likes.append(likelihood(i,intra_minor,6,acc_stringent))
    max_likes=[-1*x for x in likes] # convert back to positive log likelihood.
    Ne=LL[max_likes.index(max(max_likes))]
    return(Ne)

def parelleFit(filename,sim_dir):
    simulation = sim_dir+"/"+filename
    if simulation.endswith(".csv"):
        return(ML_fit(simulation,acc_stringent))

def main():

    args=parser.parse_args()
    sim_dir = os.path.abspath(args.inputDir[0])
########## Read in accruacty table #################
    acc_stringent = pd.read_csv("../data/reference/accuracy_stringent.csv")
######### Fit each simulation #######################

    #c= 0
    num_cores = multiprocessing.cpu_count()

    files = os.listdir(sim_dir)
    sim_fits = Parallel(n_jobs=num_cores)(delayed(parelleFit)(filename,sim_dir) for filename in files)

                #c+=1
############### write output ###########################
    out_LL =[]
    for i in range(0,len(sim_fits)):
        out_LL.append({"Simulation":i,"Ne":sim_fits[i]})
        out_pd = pd.DataFrame(out_LL)

    out_pd.to_csv(args.outCsv[0])

main()
