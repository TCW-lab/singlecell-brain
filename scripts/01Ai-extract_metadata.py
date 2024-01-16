
# #installation of anndata/scanpy : https://scanpy.readthedocs.io/en/stable/installation.html
# micromamba create -n singlecell
# micromamba activate singlecell
# micromamba install -c conda-forge scanpy python-igraph leidenalg
# pip3 install scanpy

#!micromamba activate singlecell
#!python

import numpy as np
import anndata
import scanpy as sc
#import mudata as md
from scipy import sparse
#import pooch
import matplotlib.pyplot as plt
#import pooch
from scipy.sparse import csr_matrix

# #test on endothelial cells
# #load h5ad data
# datapath='/usr2/postdoc/adpelle1/tcwlab-load/ref-data/SEAAD/DLPFC/endothelial.h5ad'
# adata = anndata.read_h5ad(datapath)
# adata
# adata.layers
# #check metadata
# adata.obs.head() #there is so save it in csv.gz
# adata.obs.to_csv('outputs/01-SEAAD_data/DLPFC/endothelial_metadata.csv.gz',compression='gzip')
# 
# #check matrix
# print(adata.X) 
# print(adata.layers['raw']) 


#strat: get the full h5ad on aws, create one h5ad by cell type, and transform to Seurat.
#get MTG and DLPFC on AWS [on R script]
#create one h5ad by cell type
#follow doc of anndata, subseting https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html
# #DLPFC
# datapath='/projectnb/tcwlab-load/ref-data/SEAAD/DLPFC/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad'
# adata = anndata.read(datapath,backed='r')
# 
# adata
# adata.layers
# 
# #check metadata
# adata.obs.head() #there is so save it in csv.gz
# adata.obs.to_csv('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz',compression='gzip')



#MTG
datapath='/projectnb/tcwlab-load/ref-data/SEAAD/MTG/SEAAD_MTG_RNAseq_final-nuclei.2023-07-19.h5ad'
adata = anndata.read(datapath,backed='r')
adata
adata.layers

#save metadata
adata.obs.head() 
adata.obs.to_csv('outputs/01-SEAAD_data/MTG/all_final_RNAseq_nuclei_metadata.csv.gz',compression='gzip')
