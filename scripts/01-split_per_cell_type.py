
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
#DLPFC
datapath='/projectnb/tcwlab-load/ref-data/SEAAD/DLPFC/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad'
adata = anndata.read_h5ad(datapath)

adata
adata.layers

#check metadata
adata.obs.head() #there is so save it in csv.gz
adata.obs.to_csv('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz',compression='gzip')

# #subset endothelial to test
# adata.obs_keys()
# adata.obs['Subclass'] #SubClass
# adata_ct=adata[adata.obs['Subclass']=='Endothelial',:]
# adata_ct.write('outputs/01-SEAAD_data/DLPFC/endothelial.full.h5ad',compression='gzip')
# adata_ct = anndata.read_h5ad('outputs/01-SEAAD_data/DLPFC/endothelial.full.h5ad')

#save anndata by celltype
for ct in adata.obs['Subclass'].cat.categories:
  celltype=ct.replace(' ','_')
  print(celltype)
  adata_ct=adata[adata.obs['Subclass']==celltype,:]
  filename='outputs/01-SEAAD_data/DLPFC/'+celltype+'.full.h5ad'
  adata_ct.write(filename,compression='gzip')
  


#MTG
datapath='/projectnb/tcwlab-load/ref-data/SEAAD/MTG/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad'
adata = anndata.read_h5ad(datapath)
adata
adata.layers

#save metadata
adata.obs.head() 
adata.obs.to_csv('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz',compression='gzip')

#save anndata by celltype
for ct in adata.obs['Subclass'].cat.categories:
  celltype=ct.replace(' ','_')
  print(ct)

  adata_ct=adata[adata.obs['Subclass']==celltype,:]
  filename='outputs/01-SEAAD_data/MTG/'+celltype+'.full.h5ad'
  adata_ct.write(filename,compression='gzip')



