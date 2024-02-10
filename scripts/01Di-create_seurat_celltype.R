
#Get  seurat obejct for each cell type
out<-'outputs/01-SEAAD_data/'
library(data.table)
setDTthreads(0)
library(Seurat) #seurat V4, not V5
#if youalready have v5 install, to get back to v4 version:
# remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(BPCells)
library(rhdf5)
library(stringr)
Sys.setenv(USE_SYSTEM_LIBGIT2=1)


#Analysis####

# #DLPFC
# out1<-file.path(out,'DLPFC')
# 
# files_paths=list.files('outputs/01-SEAAD_data/DLPFC',pattern = '.h5ad',full.names = T)
# 
# mtd<-fread('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz')  #metadata containing the celltype/donors annotation. should be extract from the anndata object before, see LINKTOADD
# mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column
# 
# cell_group='Supertype'
# n_cells_max=200000
# h5ls(files_paths[1]) #to know where the raw counts matrix is stored
# 
# raw_counts_location<-'layers/UMIs'
# 
# for(file in files_paths){
#   message(file)
#   
#   ct<-str_remove(basename(file),'\\.h5ad')
#   data <- open_matrix_anndata_hdf5(file,group = raw_counts_location) #check that UMIs matrix is stored in layers/UMIs using h5ls(path) or reading the h5ad file on python using anndata
#   mtdf<-mtd[colnames(data),on='cell_id']
#   
# 
#   #split data according to subtype if more than n_cells_max
#   mtdf[,cell_group:=lapply(.SD,as.factor),.SDcols=cell_group]
#   group_names<-names(table(mtdf$cell_group))
#   groups<-list()
#   i<-1
#   while(length(group_names)>0){
#     
#     groups[[i]]<-group_names[cumsum(table(mtdf$cell_group)[group_names])<n_cells_max]
#     group_names<-setdiff(group_names,groups[[i]])
#     mtdf[cell_group%in%groups[[i]],cellsplit:=i]
#     i<-i+1
#     
#   }
#   table(mtdf$cellsplit)
#   
#   n_groups<-length(groups)
#   
#   for(i in 1:n_groups){
#     dataf<-as(data[,mtdf[cellsplit==i]$cell_id],'dgCMatrix')
#     
#     
#     objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[cellsplit==i],row.names = 'cell_id'))
#     
#     if(n_groups>1){
#       saveRDS(objf,file.path(out1,paste0(ct,'_set',i,'.rds')))
#       
#     }else{
#       saveRDS(objf,file.path(out1,paste0(ct,'.rds')))
#       
#     }
#     
#   }
# 
#  
#   
# }


#MTG
out1<-file.path(out,'MTG')

files_paths=list.files(out1,pattern = '.h5ad',full.names = T)

mtd<-fread(file.path(out1,'all_final_RNAseq_nuclei_metadata.csv.gz')) #metadata containing the celltype/donors annotation. should be extract from the anndata object before, see LINKTOADD
mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column

cell_group='Supertype'
n_cells_max=200000
h5ls(files_paths[1]) #to know where the raw counts matrix is stored

raw_counts_location<-'layers/UMIs'

for(file in files_paths){
  message(file)
  
  ct<-str_remove(basename(file),'\\.h5ad')
  data <- open_matrix_anndata_hdf5(file,group = raw_counts_location) #check that UMIs matrix is stored in layers/UMIs using h5ls(path) or reading the h5ad file on python using anndata
  mtdf<-mtd[colnames(data),on='cell_id']
  
  #split data according to subtype if more than n_cells_max
  mtdf[,cell_group:=lapply(.SD,as.factor),.SDcols=cell_group]
  group_names<-names(table(mtdf$cell_group))
  groups<-list()
  i<-1
  while(length(group_names)>0){
    
    groups[[i]]<-group_names[cumsum(table(mtdf$cell_group)[group_names])<n_cells_max]
    group_names<-setdiff(group_names,groups[[i]])
    mtdf[cell_group%in%groups[[i]],cellsplit:=i]
    i<-i+1
    
  }
  table(mtdf$cellsplit)
  
  n_groups<-length(groups)
  
  for(i in 1:n_groups){
    dataf<-as(data[,mtdf[cellsplit==i]$cell_id],'dgCMatrix')
    
    
    objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[cellsplit==i],row.names = 'cell_id'))
    
    if(n_groups>1){
      saveRDS(objf,file.path(out1,paste0(ct,'_set',i,'.rds')))
      
    }else{
      saveRDS(objf,file.path(out1,paste0(ct,'.rds')))
      
    }
    
  }
  
  
  
}
