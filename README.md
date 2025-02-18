# Singlecell-brain analysis
Pipelines to preprocess Human (postmortem) brain single cell genomics data


# Dataset analysis 
We looking into different single nuclei genomics from postmortem brain : 
from SEA-AD , ROSMAP, and Gazestani et al cohort

## SEA-AD
original and preprocess data can be found at `/projectnb/tcwlab/ShareSpace/SEA-AD`
It contains:
- snRNA-seq for  2 brain regions : MTG and DLPFC
- snATAC-seq data for MTG region

Preprocessing notebooks:
Data presentation : [SEAAD_snRNA_data.ipynb](notebooks/SEAAD_snRNA_data.ipynb)  
Data QC/Preprocessing : [SEAAD_QC.ipynb](notebooks/SEAAD_QC.ipynb)  
Data pseudobulk generation : [SEAAD_pseudobulk.ipynb](notebooks/SEAAD_pseudobulk.ipynb)  


Astrocyte Pseudobulk Differential Expression analysis : [pseudobulk_deseq2.ipynb](notebooks/pseudobulk_deseq2.ipynb)  
Astrocyte fgsea analysis : [fgsea_analysis.ipynb)](https://github.com/TCW-lab/singlecell-brain/blob/main/notebooks/fgsea_analysis.ipynb)


## ROSMAP
orginal (preprocess to be added, bug Alexandre if it is still not done) data can be found at `/projectnb/tcwlab/ShareSpace/ROSMAP`
It contains:
- snRNA-seq for DLPFC for ~400 samples from MIT/Kellis.
- snRNA-seq for 4 brain regions for ~80 samples in `MultiRegion` from MIT/Kellis
- snATAC-seq for ~80 samples in DLPFC)
- snRNA-seq for DLPFC for ~400 samples from Columbia/De Jager (To be added).

  
Preprocessing notebooks/scripts:
snATAC data preprocessing [04-ROSMAP_MIT_ATAC.R](https://github.com/TCW-lab/singlecell-brain/blob/main/scripts/04-ROSMAP_MIT_ATAC.R)

Preprocessing of snRNA/ATAC associated WGS based genotyping [05-ROSMAP_MIT_WGS.R](https://github.com/TCW-lab/singlecell-brain/blob/main/scripts/05-ROSMAP_MIT_WGS.R)


## Gazestani et al
Some checks in the dataset: [02-Gazestani_EarlyAD.R](https://github.com/TCW-lab/singlecell-brain/blob/main/scripts/02-Gazestani_EarlyAD.R)


# Some  Documentations

**How to :**

Transform large h5ad dataset to Seurat : [notebooks/01-convert_SEAAD_h5ad_to_SeuratV5](notebooks/01-convert_SEAAD_h5ad_to_SeuratV5.md)

downsample large single cell data : [notebooks/downsampling_singlecell_data.md](notebooks/downsampling_singlecell_data.md)


## Project Status (Feb 2024)
- Generate ROSMAP snATAC matrix per cell type (annotated thanks to SEA-AD DLPFC scRNA)
- associate with WGS based genotyping 

