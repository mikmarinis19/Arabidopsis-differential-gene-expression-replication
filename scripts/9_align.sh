#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 3
#SBATCH --mem=15G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[1-18]%10

hostname
date

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.16.1

INDIR=../trimmed_fastq
OUTDIR=../alignments
mkdir -p ${OUTDIR}

#Get genome index
INDEX=../genome/hisat2_index/Arab

#Get accessionlist.txt

ACCLIST=../metadata/accessionlist.txt

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${ACCLIST})

# Run HISAT2 alignment for single-end FASTQ
# The data from my paper is single-end
hisat2 \
    -p 4 \
    -x ${INDEX} \
    -U ${INDIR}/${SAMPLE}_trimmed.fastq \
    | samtools sort -@ 1 -T ${OUTDIR}/${SAMPLE}.tmp -o ${OUTDIR}/${SAMPLE}.bam

# Index the sorted BAM
samtools index ${OUTDIR}/${SAMPLE}.bam

