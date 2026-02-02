#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# FASTQC
#################################################################
module load fastqc/0.12.1
module load parallel/20180122

# set input/output directory variables
INDIR=../trimmed_fastq/
REPORTDIR=../trimmed_fastqc/
mkdir -p $REPORTDIR

# run FASTQC in parallel on all .fastq files
ls $INDIR/*.fastq | parallel -j 10 "fastqc --outdir $REPORTDIR {}"


