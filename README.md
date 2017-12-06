# Influenza evolution within and between naturally occuring acute infections
This repository holds the analysis used in *citation,2018* it relies heavily on our other repository [variant_pipeline](https://github.com/lauringlab/variant_pipeline).

#Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            # raw and primary data, are not changed once created
    |  |- reference/  # reference fasta files to be used in in sequence alignment
    |  |- raw/         # raw data, will not be altered. In practice this is where the raw fastq files go. They will need to be downloaded from the SRA. 
    |  |- process/     # cleaned data, will not be altered once created. During anlysis this includes intermediate files used in variant calling as well as consensus sequence processing
    |
    |- scripts/           # any programmatic code
    |- results         # all output from workflows and analyses
    |  |- Figures/     # graphs, likely designated for manuscript figures
    |  |- *.Rmd  # exicutable R markdown to make the figures 
    +- Makefile        # executable Makefile for this study
    
  --------
# Dependencies    
The analysis expects the variant_pipeline repository to be your home directory. 
```
    cd ~/
    git clone https://github.com/lauringlab/variant_pipeline.git
```
This pipeline requires [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and the R package [DeepSNV](https://www.bioconductor.org/packages/release/bioc/html/deepSNV.html).

More information on the dependencies needed to run the varaint\_pipline can be found [here](https://github.com/lauringlab/variant_pipeline)

The python code used in this analysis has been tested in version 2.7.11 (although it will probably work in other versions) and relies on the following modules. The version used for each module is also listed.

```
biopython -- 1.67
pandas -- 0.18.1
```

The consensus sequence analysis and antigenic analysis relies on [muscle](http://www.drive5.com/muscle/downloads.htm). It expects the exicutable to be named "muscle". Please change the path to this exicutable in the Makefile.

# Reproducing the analysis

We can use the commands in the Makefile to do everyhing from downloading the fastq files from the SRA all the way through making the figures. There are 4 stages to the analysis. 1) Downloading the fastq files from the SRA, 2) Primary analsysis - Calling iSNV in each sample using the variant_pipeline repository referenced above, 3) Secondary analysis - maniputating the data and running the models 4) making the figures. All the intermediate files needed to run the secondary analysis are included in this repository. 

Please note that in many places we refer to the genomic segment "NA" as "NR". This is to avoid complications in R as "NA" is a special term. 


# NOTE : WORK ADDING THESE STEPS TO THE MAKEFILE IS CURRENTLY IN PROGRESS - (12/6/2017)

## Downloading raw data

Becasue we use a plasmid control to estimate the lane specific error rates used to identify iSNV it is important that the fastq files from the SRA are split into separate directories for each Hiseq lane. This stage makes a file list for each run (if needed) and then downloads the .sra files using wget. The sra files are then converted to fastq files and are renamed to match the sample names used in the rest of the analysis. All of this is achieved my running the "download" phony target.

```
make download
```


## Processing the raw data

The following command launches the analyis pipeline from https://github.com/lauringlab/variant_pipeline. This pipeline is based in [bpipe](http://bpipe-test-documentation.readthedocs.io/en/latest/) and is smart enough to remember what commands have been run in the event of failure. It also logs the commands that were run("./data/processed/\*/commandlog.txt"). These steps are time and memory intensive. These commands output a number of large intermediate files. The important files that are used in down stream analysis are the concatenated variant calls "./data/processed/\*/Variants/all.sum.csv". The concatenated coverage files "./data/processed/\*/deepSNV/all.coverage.csv", and the consensus sequences "./data/processed/\*/parsed_fa/\*.fa". This stage of the analysis can be run using the "primary" phony target, and requires that the fastq files are downloaded and named appropriately.
```
make primary
```


## Secondary analysis

The secondary analysis is broken up into 3 separate stages. Each can be run separately or all three can be run at the same time. 

All stages are run by the phony target "secondary"

```
make secondary
```

