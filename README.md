# 🦴 SAHZU OS Transcriptomics 🧬
This repository includes codes for RNA-sequencing data preprocessing, analysis, and generating plots for the manuscript entitled "Tumour persistence programmes and immune remodelling characterise osteosarcoma chemotherapy response". This repository is intended for editor and reviewer use during peer review process and is to be revised pending review outcomes. 
1. **Data pre-processing and analyses** - Contains code for scRNA-seq, bulk RNA-seq and spatial transcriptomic data preprocessing and analysis. 
  - **scRNA**: Contains code for data pre-processing, quality control, cell type derivation and annotation, and downstream analyses. 
  - **Bulk tissue RNA**: Contains code for in-house dataset integration and downstream analyses. 
  - **Spatial**: Contains code for integrated analyses for in-house 10x Visium ST data. 
2. **Plot generation** - Contains code for generating both main and supplementary figure panels in the manuscript, if not already provided in the analytic pipeline.

## 📊 Data Availability
All sequencing data have been deposited in the China National Center for Bioinformation Genome Sequence Archive for Human (https://ngdc.cncb.ac.cn/gsa-human/) with the below accession numbers. 
  - **Bulk**: HRA009869 (GSA-Human)
  - **scRNA**: HRA009821 (GSA-Human)
  - **Spatial**: HRA017052 (GSA-Human)

## 📂 Software versions
### R 
R (v4.4.2), Seurat (v.5.2.1), miloR (v.2.3.1), inferCNV (v.1.22.0), UCell (v2.10.1), CytoTRACE2 (v.1.1.0), org.Hs.eg.db (v3.20.0), msigdbr (v10.0.1), clusterProfiler (v4.14.6), Monocle3 (v1.3.7), Slingshot (v2.10.0), SCP (v0.5.6), CellChat (v2.2.0), nichenetr (v2.2.0), sva (v3.54.0), DESeq2 (v1.46.0), limma (v3.62.2), GSVA (v2.0.7), survival (v3.8.3), survminer (v0.5.0), ezcox (v1.0.4), FactoMineR (v2.11), TCGAbiolinks (v2.34.1), CARD (v1.1), and here (v1.0.1). 

### Python
Python (v3.8.20), CellPhoneDB (v5.0.1), pySCENIC (v0.12.1), and Uphyloplot2 (v.2.3)

## 📜 License
The code in this repository is released under the MIT License (see `LICENSE`). Note that the EcoTyper helper files in `01_Data_preprocess_and_analysis/02_scRNA/05_Cross_cell_type_analysis/lib/` remain under the terms of the original EcoTyper repository (see `lib/README.md`).

## 📝 Correspondence
Any queries should be addressed to the corresponding author Prof Zhaoming Ye yezhaoming@zju.edu.cn
