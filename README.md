# Singlecell-brain analysis
Pipeline to preprocess Human (postmortem) brain single cell genomics data

# Project Status (Feb 2024)
- Generate ROSMAP snATAC matrix per cell type (annotated thanks to SEA-AD DLPFC scRNA)
- associate with WGS based genotyping 

# Analysis 
We looking into 3 different single cell brain datasets : 
from SEA-AD , ROSMAP, and Gazestani et al cohort

## SEA-AD
Data QC/Preprocessing : [SEAAD_QC.ipynb](notebooks/SEAAD_QC.ipynb)

Astrocyte Pseudobulk Differential Expression analysis : [Pseudobulk_DESeq2_analysis_of_SEAAD_data.ipynb](notebooks/Pseudobulk_DESeq2_analysis_of_SEAAD_data.ipynb)

## Gazestani et al
Some checks in the dataset: [02-Gazestani_EarlyAD.R](https://github.com/TCW-lab/singlecell-brain/blob/main/scripts/02-Gazestani_EarlyAD.R)


## ROSMAP
snATAC data preprocessing [04-ROSMAP_MIT_ATAC.R](https://github.com/TCW-lab/singlecell-brain/blob/main/scripts/04-ROSMAP_MIT_ATAC.R)

Preprocessing of snRNA/ATAC associated WGS based genotyping [05-ROSMAP_MIT_WGS.R](https://github.com/TCW-lab/singlecell-brain/blob/main/scripts/05-ROSMAP_MIT_WGS.R)


# Some  Documentations

**How to :**

Transform large h5ad dataset to Seurat : [notebooks/01-convert_SEAAD_h5ad_to_SeuratV5](notebooks/01-convert_SEAAD_h5ad_to_SeuratV5.md)

downsample large single cell data : [notebooks/downsampling_singlecell_data.md](notebooks/downsampling_singlecell_data.md)


