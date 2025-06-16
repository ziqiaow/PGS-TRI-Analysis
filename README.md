# PGS-TRI-Analysis

This repository contains the source code to reproduce analyses in "Estimation of Direct and Indirect Polygenic Effects and Gene-Environment Interactions using Polygenic Scores in Case-Parent Trio Studies."

The official software package with instructions and examples to perform PGS-TRI is in [PGS.TRI](https://github.com/ziqiaow/PGS.TRI/tree/main).

Preprint manuscript is released at [medRxiv](https://www.medrxiv.org/content/10.1101/2024.10.08.24315066v1). 

## File structures

### 1_main_and_extended_figures
The directory contains the codes used to plot main and extended figures.

### 2_simulation
This directory contains the codes used to perform simulation analysis in this manuscript.
#### Simulate PGS based on the Model
This directory contains the codes for simulations using the proposed model framework.
#### Simulation using the UK Biobank Data
This directory contrain the codes for simulations using the UK Biobank data.

### 3_data_applications
This directory contains the codes used to perform data analysis of direct and indirect effect, PGSxE analyses of ASD and OFCs, as well as the multi-trait analyses for ASD in the SPARK consortium and GENEVA data in this manuscript. This also includes the codes for the individual genotypic TDT test for OFCs.

### 4_TWAS_MWAS
This directory contains the codes used to construct the genetic scores of metabolites and gene expression levels and perform data analysis of transcriptome-wide and metabolome-wide analysis using PGS-TRI of ASD and OFCs in the SPARK consortium and GENEVA data in this manuscript.

### 5_R_Github
This directory contains the source codes of the PGS-TRI software [PGS.TRI](https://github.com/ziqiaow/PGS.TRI/tree/main).

### 6_revision1
This directory contains the source codes of the analyses done for revision version 1 06/16/2025. Additional analyses performed including the new 18,383 trios in SPARK v3. PGS standardization using PC projections based on 1000G+HGDP reference data. PGSxPC interactions. Direct and indirect effect estimations for ASD and various pleiotropy traits. GTEx v8 TWAS analysis in brain tissues. Simulation studies of assessing selection bias and simulating 20 generations of assortative mating in phenotypes using snipar. 

# Questions
Please send your questions/suggestions to Ziqiao Wang (zwang389@jhu.edu).
