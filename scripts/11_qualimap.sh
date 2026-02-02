#!/bin/bash 
#SBATCH --job-name=qualimap
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# calculate stats on alignments
##################################
# this time we'll use qualimap

# load software--------------------------------------------------------------------------
module load qualimap/2.2.1
module load parallel/20180122

# input, output directories--------------------------------------------------------------

INDIR=../alignments
OUTDIR=../qualimap_reports
mkdir -p $OUTDIR

# gtf annotation is required here
GTF=../genome_athaliana/Arabidopsis_thaliana.TAIR10.62.gtf 

# accession list
ACCLIST=../metadata/accessionlist.txt

# run qualimap in parallel
cat $ACCLIST | \
parallel -j 5 \
    qualimap \
        rnaseq \
        -bam $INDIR/{}.bam \
        -gtf $GTF \
        -outdir $OUTDIR/{} \
        --java-mem-size=2G  
