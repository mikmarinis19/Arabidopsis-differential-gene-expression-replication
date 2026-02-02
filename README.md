# Replication of Thermomorphogenesis RNA-seq Analysis in Arabidopsis

## Overview
This repository contains a replication of the RNA-seq analysis from
_Nutrient levels control root growth responses to high ambient temperature in plants_ (citation below).

Thix study investigates thermomorphogenesis in Arabidopsis thaliana, focusing on the roles of
HY5 and NRT1.1 in root growth responses to elevated temperature.

## Biological Background
Thermomorphogenesis refers to plant growth responses to elevated ambient temperature,
including enhanced root and shoot growth, curvature, and early flowering.
Key regulators include:
- **HY5**: a master transcriptional regulator of root thermomorphogenesis
- **NRT1.1**: a nitrate transporter downregulated by HY5 under heat stress

This project examines how loss-of-function mutants (hy5-215 and chl1-5) alter global gene
expression compared to wild-type Col-0 plants.

## Experimental Design
- **Organism**: Arabidopsis thaliana
- **Genotypes**: Col-0, hy5-215, chl1-5
- **Conditions**: 21°C (control) vs 28°C (heat)
- **Replicates**: 3 biological replicates per condition (18 samples total)
- **Data source**: GEO accession GSE262197

## Methods Summary
- Quality control: FastQC, MultiQC
- Trimming: Trimmomatic
- Alignment: HISAT2
- Counting: HTSeq
- Differential expression: edgeR
- Functional analysis: clusterProfiler (GO enrichment)

## Key Results
- Differential expression results closely matched the original study
- DEG counts differed by <35 genes across all genotypes
- PCA and heatmap analyses showed clear clustering by genotype and temperature
- GO enrichment recapitulated key thermomorphogenesis-related pathways

## Repository Structure
- `scripts/` – RNA-seq analysis scripts
- `results/` – figures and tables
- `data/` – processed count matrices
- `docs/` – full written report

## Citation
Lee, S., Showalter, J., Zhang, L. et al. Nutrient levels control root growth responses to high ambient temperature in plants. Nat Commun 15, 4689 (2024). https://doi.org/10.1038/s41467-024-49180-6
