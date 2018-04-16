
# coding: utf-8

# # Dn Ds on our data
# 
# Interpretting Dn/Ds is challenging. It is not robust on short time scales. That said we see an overaboundance of synonymous mutations and have been asked to correct for the number of NS and S mutations. So we do it, keeping in mind the caveots.
# 
# _NB: check http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/ for reference_

# In[2]:


import numpy as np
import pandas as pd
import copy 
from matplotlib import pyplot as plt
import os
import tempfile
import sys
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import codonalign
from ast import literal_eval
import re
get_ipython().run_line_magic('matplotlib', 'inline')



# In[3]:


def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    'fastafile' should specify the name of a FASTA file.

    This function reads all sequences from the FASTA file.  It returns the
        list 'headers_seqs'.  This list is composed of a seq_record objects.
    """
    seqs =[]
    header = None
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq_record.seq.alphabet=IUPAC.unambiguous_dna
        seqs.append(seq_record)

    return seqs




# In[4]:


test= ReadFASTA("../data/reference/NY.OR.main.fa")


# In[5]:


PB2=test[1]
c = codonalign.CodonSeq(str(PB2.seq))


# In[6]:


c


# In[7]:


cList=codonalign.codonseq._get_codon_list(c)
codonalign.codonseq._count_site_NG86(cList[:-1],k=1)


# # Outine
# 
# 1) Read in meta data and filter to samples that qualified for iSNV identification
# 
# 2) Count the number of iSNV in each sample
#         - Count NS and S
#         - Count by segment
#         
# 3) For each sample for each segment
#         - Get OR.
#         - Calculate NS and S
#     - Add to table of NS and S for each segment
#     
# 4) Calculate Dn/Ds for each segment and total genome.
# 
# ### To think about
# What reading frames are we using - at this point just the cononical
#     NS and S was defined as NS in any OR.
# Should we only use 1 sample per illness? - Probabily but not done at this point.
# 
# We assume each mutation is in a unique codon  - almost certainly true but easily checked.
#    

# ## 1) Meta file
# 
# This work is much easier for me to do in R so we'll use that here

# In[8]:


meta = pd.read_csv("../data/processed/secondary/meta_for_ns.s_calc.csv")


# In[9]:


meta = meta.loc[meta.snv_qualified==True]


# In[10]:


sys.path.append("/Users/jt/lauring_lab_repos/variant_pipeline/scripts/")
from fasta_functions import StripGapsToFirstSequence, Align
def get_meta(SPECID,meta):
    """This function takes a SPECID and data frame with meta data and
    returns a dictionary of the meta data for that sample."""
    ID = meta.loc[meta.SPECID==SPECID,"Id"].unique()[0]
    ID = ID.split(".")[0]
    RUN = meta.loc[meta.SPECID==SPECID,"run"].unique()[0]
    season = meta.loc[meta.SPECID==SPECID,"season"].unique()[0]
    ENROLLID = meta.loc[meta.SPECID==SPECID,"ENROLLID"].unique()[0]
    HOUSE_ID = meta.loc[meta.SPECID==SPECID,"HOUSE_ID"].unique()[0]
    if RUN=="vic":
        RUN="victoria"
    if RUN=="vic_2":
        RUN="victoria_2"
    return({"Id":ID,"run":RUN,"season":season,"enrollid":ENROLLID,"house_id":HOUSE_ID})

def run_to_OR(run):
    
    """This function takes in the name of sequencing run from this study and returns the relative path
    to the fasta files with the OR for that sample."""
    conversion={"perth":"../data/reference/perth.OR.fa",
               "perth_2": "../data/reference/perth.OR.fa",
               "cali09":"../data/reference/cali09.OR.fa",
               "cali09_2":"../data/reference/cali09.OR.fa",
               "victoria":"../data/reference/victoria.OR.fa",
               "victoria_2":"../data/reference/victoria.OR.fa",
               "HK_1":"../data/reference/NY.OR.fa",
               "HK_2":"../data/reference/NY.OR.fa",
               "HK_6":"../data/reference/NY.OR.fa",
               "HK_7":"../data/reference/NY.OR.fa",
               "HK_8":"../data/reference/NY.OR.fa"}
    return(conversion[run])

def trim_to_coding(fasta,SPECID,meta):
    """This funciton trims a sample fasta sequence to the reading frame defined in a separate fasta file
    First the funciton takes in the SPECID and looks up where the sample fasta file will be using the meta data
    available in the meta argument. For each sequence in the reference the function looks for the same sequence name
    in the sample fasta. It then alings this and trims the gaps so we are left with just the OR."""
    
    samp_meta = get_meta(SPECID,meta)
    
    regions = run_to_OR(samp_meta["run"])
    
    coding = ReadFASTA(regions)
    # cycle through fasta 
    OR=[]
    for gene in fasta:
        seg_id=gene.id  
        gene.description = "this is the test sample" + gene.description
        for code in coding:
            code_id = code.id
            if seg_id in code_id:
                #print('working with %s' %seg_id) #and seg_id=="NR":
                code_gene=Align([code,gene],"/Users/jt/muscle3.8.31/")
                code_gene_trimmed = StripGapsToFirstSequence(code_gene)
                code_gene_trimmed.id=code_id
                code_gene_trimmed.name=code_id
                OR.append(code_gene_trimmed)
        
    return(OR)
def get_seq(specid_list,meta_run):
    """This function takes a list of SPECIDs and data frame with meta data for the samples.  It only returns the consensus sequences.
    For each SPECID we get the meta data in dictionary form,g et the consensus fasta file, 
    Trim the new sequences to the coding regions and return a dictionary of the results indexed by SPECID."""
    sequences={}
    for specid in specid_list:
        meta = get_meta(specid,meta_run)
        fa = "../data/processed/"+meta["run"]+"/parsed_fa/"+meta["Id"]+".removed.parsed.fasta"
        seq = ReadFASTA(fa)
    
        for seg in seq:
            seg.name=specid
        
        seg_coding = trim_to_coding(seq,specid,meta_run)
        specid_key = '%s_%s_%s_%s_consensus' % (specid,meta["enrollid"],meta["house_id"],meta["season"])
        sequences[specid_key] = seg_coding 
    return(sequences)

def get_NS_s(seqs):
    counter = {"PB2":[0,0],"PB1":[0,0],"PB1-F2":[0,0],
               "PA":[0,0],"PA-X":[0,0],"HA":[0,0],"NP":[0,0],
               "NR":[0,0],"M1":[0,0],"M2":[0,0],"NS1":[0,0],"NS2":[0,0]}
    #seqs = get_seq(specid_list,meta_run)
    # cycle through the samples
    
    errors = pd.DataFrame(columns=["SPECID", "Segment","Error"])

    for sample in seqs:
        sample_seq = seqs[sample]# a list of seqrecords
        for sequence in sample_seq: # cycle through seqrecords
            seq_name = sequence.name
            assert len(str(sequence.seq)) % 3 == 0, "Sequence length is not a triple number"
            try:
                codon_sequence = codonalign.CodonSeq(str(sequence.seq),gap_char='-')
            except ValueError as ve:
                print("ValueError probabable gap - skipping %s in sample %s Remove any mutations from this OR in count" % ( seq_name, sample))
                errors = errors.append({"SPECID": sample.split('_')[0],
                                "Segment":  seq_name,
                                "Error" : str(ve).split("(")[1].split(")")[0]}, ignore_index=True)
            else:
                cList = codonalign.codonseq._get_codon_list(codon_sequence)
                try:
                    counts = codonalign.codonseq._count_site_NG86(cList[:-1],k=1) # remove the stop codon
                except KeyError as e: 
                    print("Codon Key error. Found %s in OR \n skipping %s in sample %s Remove any mutations from this OR in count" % (e, seq_name, sample))
                    errors = errors.append({"SPECID": sample.split('_')[0],
                                "Segment":  seq_name,
                                "Error" : e}, ignore_index=True)
                else:
                            counter[seq_name][0] += counts[0]
                            counter[seq_name][1] += counts[1]
    return(counter,errors)


# In[11]:


sequences = get_seq(meta.SPECID,meta)


# In[11]:


Ns_s,e = get_NS_s(sequences)


# In[12]:


l=[]
for x in e.SPECID:
    y = meta.pcr_result[meta.SPECID==x]
    l.append(y.asobject[0])
    
e["pcr_result"] = l


# In[13]:


e.loc[e.pcr_result=="A/H3N2"]


# ## PB1-F2 Stop
# 
# Sample MH5300 is H1N1 and has a stop codon in the PB1-F2 OR. This is confirmed in the best Blast alignment and in a note in CY188895.1 . The function should ouput any other errors. I will confirm these and then adjust the NS S counts accrodingly.
# 
# This same mutation is found in many samples. MH2516, MH2527 and MH2942 are the only H3N2 with the mutation.
# 

# In[14]:


H3N2_stop = str(sequences["MH2516_330243_3061_2012-2013_consensus"][2].seq)
H3N2_stop


# In[15]:


codon_sequence = codonalign.CodonSeq(H3N2_stop,gap_char='-')
x = codonalign.codonseq._get_codon_list(codon_sequence)

x.index("TAA")
#len(x)


# CY171005.1 matches this sequence. I think we will ignore PB1-F2 and PA-X

# ## PB2 gap
# 
# One sample has a frame shift deletion of 2 nt near the end of PB1. This is confirmed by looking at the igv viewer. This is a victoria sample 1230. 

# In[16]:


gap = sequences["MH2436_331045_3075_2012-2013_consensus"][1].seq


# In[17]:


str(gap)


# ## NS1 issues
# 
# I have checked the following samples. They have a premeture stop codon that terminates the NS1 AA early.
# MH1782  -3 codons early
# MH0922  -10 codons early

# In[18]:



NS1_stop = sequences["MH0922_320110_2028_2011-2012_consensus"]
#NS1_stop = sequences["MH1782_320468_2117_2011-2012_consensus"]


# In[19]:


str(NS1_stop[10].seq)


# Blast these sequences and identify the mutation.

# In[20]:


Ns_s_df = pd.DataFrame(Ns_s)
Ns_s_df = Ns_s_df.transpose()
Ns_s_df.columns = ["S","NS"]


# In[21]:


Ns_s_df.to_csv("../data/processed/secondary/NS_S_site.csv")


# In[22]:


Ns_s_df


# In[23]:


codonalign.codonseq._count_site_NG86(["CCC"],k=1)


# In[12]:


330663.333333/88742.666667


# In[16]:


codonalign.codonseq._count_site_NG86(cList[:-1])


# In[22]:


codonalign.codonseq._count_site_NG86([cList[1]])


# In[30]:


HA=test[3]


c = codonalign.CodonSeq(str(HA.seq))
HAcodonList=codonalign.codonseq._get_codon_list(c)



# In[31]:


len(HAcodonList)


# In[36]:


codonalign.codonseq._count_site_NG86(HAcodonList[:-1],k=1/3.6)


# In[33]:


1339.666666666669/358.33333333333223


# In[37]:


.75/3.7


# In[59]:


HA=test[4]


c = codonalign.CodonSeq(str(HA.seq))
HAcodonList=codonalign.codonseq._get_codon_list(c)
x = codonalign.codonseq._count_site_NG86(HAcodonList[:-1])


# In[60]:


x[1]/x[0]


# In[61]:


HA


# In[62]:


len(test)

