#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=5G
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

module load MultiQC/1.9

STATS=../samtools_stats
QUALIMAP=../qualimap_reports
OUTDIR=../multiqc_alignment

# run on samtool stats output
multiqc -f -o ${OUTDIR}/samtools_stats ${STATS}

# run on qualimap output
multiqc -f -o ${OUTDIR}/qualimap ${QUALIMAP}
