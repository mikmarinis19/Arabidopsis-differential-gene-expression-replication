#!/bin/bash
#SBATCH --job-name=htseq_count
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
date

module load htseq/0.13.5
module load parallel/20180122

# Input & output directories
# I originally had some trouble with this code, so to troubleshoot I used the full paths
# Eventually the code worked

INDIR=/home/FCAM/mmarinis/final_project/alignments
OUTDIR=/home/FCAM/mmarinis/final_project/counts
mkdir -p $OUTDIR

# Accession list
ACCLIST=/home/FCAM/mmarinis/final_project/metadata/accessionlist.txt

# GTF annotation
GTF=/home/FCAM/mmarinis/final_project/genome_athaliana/Arabidopsis_thaliana.TAIR10.62.gtf

# Run htseq-count in parallel for all samples in the accession list
cat $ACCLIST | parallel -j 10 \
    "htseq-count -s no -r pos -f bam $INDIR/{}.bam $GTF > $OUTDIR/{}.counts"

