
# coding: utf-8

# # Dn Ds on our data
# 
# Interpretting Dn/Ds is challenging. It is not robust on short time scales. That said we see an overaboundance of synonymous mutations and have been asked to correct for the number of NS and S mutations. So we do it, keeping in mind the caveots.
# 
# _NB: check http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/ for reference_

# In[33]:


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



# In[7]:


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




# In[8]:


test= ReadFASTA("../data/reference/NY.OR.main.fa")


# In[24]:


PB2=test[1]
c = codonalign.CodonSeq(str(PB2.seq))


# In[27]:


c


# In[44]:


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

# In[64]:


meta = pd.read_csv("../data/processed/secondary/meta_for_ns.s_calc.csv")


# In[65]:


meta = meta.loc[meta.snv_qualified==True]


# In[158]:


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

def get_NS_s(specid_list,meta_run):
    counter = {"PB2":[0,0],"PB1":[0,0],"PB1-F2":[0,0],
               "PA":[0,0],"PA-X":[0,0],"HA":[0,0],"NP":[0,0],
               "NR":[0,0],"M1":[0,0],"M2":[0,0],"NS1":[0,0],"NS2":[0,0]}
    seqs = get_seq(specid_list,meta_run)
    # cycle through the samples
    
    for sample in seqs:
        print sample
        sample_seq = seqs[sample]# a list of seqrecords
        for sequence in sample_seq: # cycle through seqrecords
            seq_name = sequence.name
            assert len(str(sequence.seq)) % 3 == 0, "Sequence length is not a triple number"
            codon_sequence = codonalign.CodonSeq(str(sequence.seq))
            cList = codonalign.codonseq._get_codon_list(codon_sequence)
            print seq_name
            counts = codonalign.codonseq._count_site_NG86(cList[:-1],k=1) # remove the stop codon
            for seg in counter:
                if seg==seq_name:
                        counter[seg][0] += counts[0]
                        counter[seg][1] += counts[1]
    return(counter)


# In[159]:


Ns_s = get_NS_s(meta.SPECID,meta)

