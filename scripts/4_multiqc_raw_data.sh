#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Aggregate reports using MultiQC
#################################################################

module load MultiQC/1.15

INDIR=../fastqc/
OUTDIR=../multiqc

# run on fastqc output
multiqc -f -o ${OUTDIR} ${INDIR}
