#!/bin/bash
#SBATCH --job-name=fasterq_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# This code uses my accessionlist to download all the SRA fastq data that I am interesed in.
# I manually created my accessionlist, which only represents the SRA data listed as "roots".
# For the purposes of my study, it is the only data that I am interested in.  
#################################################################

# load software
module load parallel/20180122
module load sratoolkit/3.0.1

# My outdirectory is fastq and my code uses my accessionlist to select the data to download. 
OUTDIR=../fastq
    mkdir -p ${OUTDIR}
ACCLIST=../metadata/accessionlist.txt 

# use parallel to download 2 accessions at a time. 
cat $ACCLIST | parallel -j 2 "fasterq-dump -O ${OUTDIR} {}"


