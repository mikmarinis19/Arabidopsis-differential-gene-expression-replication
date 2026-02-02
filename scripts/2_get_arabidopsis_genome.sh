#!/bin/bash
#SBATCH --job-name=get_genome_athaliana
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Download genome and annotation for Arabidopsis thaliana (TAIR10) from Ensembl Plants
#################################################################

# load software
module load samtools/1.16.1

# output directory
GENOMEDIR=../genome_athaliana
mkdir -p $GENOMEDIR

# download the genome data
wget -P ${GENOMEDIR} ftp://ftp.ensemblgenomes.org/pub/plants/release-62/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gtf.gz
wget -P ${GENOMEDIR} ftp://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget -P ${GENOMEDIR} ftp://ftp.ensemblgenomes.org/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz
wget -P ${GENOMEDIR} ftp://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

# decompress files
gunzip ${GENOMEDIR}/*.gz

# generate samtools faidx indexes 
samtools faidx ${GENOMEDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
samtools faidx ${GENOMEDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa

