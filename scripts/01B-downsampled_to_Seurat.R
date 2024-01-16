
#Get downsampled seurat obejct of SEAAD h5ad files 
out<-'outputs/01-SEAAD_data/'
library(data.table)
setDTthreads(0)
library(Seurat) #seurat V4, not V5
#if youalready have v5 install, to get back to v4 version:
# remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(BPCells)
library(rhdf5)

#functions ####


#Analysis####

# #DLPFC
# files_paths=list.files('outputs/01-SEAAD_data/DLPFC',pattern = '.h5ad',full.names = T)
# 
# mtd<-fread('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz')  #metadata containing the celltype/donors annotation. should be extract from the anndata object before, see LINKTOADD
# mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column
# 
# pct_keep=0.2
# cell_group='Subclass'
# min_by_group=500
# max_by_group=2000
# h5ls(files_paths[1]) #to know where the raw counts matrix is stored
# 
# raw_counts_location<-'layers/UMIs'
# 
# brain_downsampled<-Reduce(function(x,y)merge(x,y),lapply(files_paths,function(file){
#   data <- open_matrix_anndata_hdf5(file,group = raw_counts_location) #check that UMIs matrix is stored in layers/UMIs using h5ls(path) or reading the h5ad file on python using anndata
#   mtdf<-mtd[colnames(data),on='cell_id']
#   mtdf[,to_keep:=cell_id%in%head(sample(cell_id,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
#       by=cell_group]
# 
#   dataf<-as(data[,mtdf[(to_keep)]$cell_id],'dgCMatrix')
# 
# 
#   objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[(to_keep)],row.names = 'cell_id'))
# 
#   return(objf)
#   
# }))
# brain_downsampled
# 
# saveRDS(brain_downsampled,file = file.path(out,'SEAAD_DLPFC_43k.rds'))
# 
# table(brain_downsampled$Subclass)



#MTG
files_paths=list.files('outputs/01-SEAAD_data/MTG',pattern = '.h5ad',full.names = T)

mtd<-fread('outputs/01-SEAAD_data/MTG/all_final_RNAseq_nuclei_metadata.csv.gz')  #metadata containing the celltype/donors annotation. should be extract from the anndata object before, see LINKTOADD
mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column

pct_keep=0.2
cell_group='Subclass'
min_by_group=500
max_by_group=2000
h5ls(files_paths[1]) #to know where the raw counts matrix is stored

raw_counts_location<-'layers/UMIs'

brain_downsampled<-Reduce(function(x,y)merge(x,y),lapply(files_paths,function(file){
  data <- open_matrix_anndata_hdf5(file,group = raw_counts_location) #check that UMIs matrix is stored in layers/UMIs using h5ls(path) or reading the h5ad file on python using anndata
  mtdf<-mtd[colnames(data),on='cell_id']
  mtdf[,to_keep:=cell_id%in%head(sample(cell_id,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
       by=cell_group]
  
  dataf<-as(data[,mtdf[(to_keep)]$cell_id],'dgCMatrix')
  
  
  objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[(to_keep)],row.names = 'cell_id'))
  
  return(objf)
  
}))
brain_downsampled

saveRDS(brain_downsampled,file = file.path(out,'SEAAD_DLPFC_43k.rds'))

table(brain_downsampled$Subclass)
