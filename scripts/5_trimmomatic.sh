#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[1-18]%10

hostname
date

#################################################################
# Trimmomatic
#################################################################

#!/bin/bash
#SBATCH --job-name=trim
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mim18007@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=1-18%10  # match my 18 samples

module load Trimmomatic/0.39
module load parallel/20180122

# Input/output directories
INDIR=../fastq/fastq_output 
TRIMDIR=../trimmed_fastq
mkdir -p $TRIMDIR

# adapters to trim out
ADAPTERS=/isg/shared/apps/Trimmomatic/0.39/adapters/TruSeq3-SE.fa

# Accession list
ACCLIST=../metadata/accessionlist.txt
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${ACCLIST})

# Run Trimmomatic

SAMPLE=$( sed -n ${SLURM_ARRAY_TASK_ID}p ${ACCLIST} )

java -jar $Trimmomatic SE \
    -threads 4 \
    ${INDIR}/${SAMPLE}.fastq \
    ${TRIMDIR}/${SAMPLE}_trimmed.fastq \
    SLIDINGWINDOW:4:25 MINLEN:45

