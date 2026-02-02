#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Index the Genome
#################################################################

# load software
module load hisat2/2.2.1

# input/output directories
OUTDIR=../genome/hisat2_index
mkdir -p $OUTDIR

#genome
GENOME=../genome_athaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

#Run hisat2-build
hisat2-build -p 16 $GENOME $OUTDIR/Arab
