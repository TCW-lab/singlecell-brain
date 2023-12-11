# Create Seurat Object from large H5ad Data
#When a dataset contain more than 300k cells, we cannot load the whole dataset in the cache, we need to work on it from the disk

#remotes::install_github("bnprks/BPCells")
#remotes::install_github('satijalab/azimuth', ref = 'master')

library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(stringr)
library(rhdf5)
library(data.table)


#FUNCTIONS####

GetH5ADCellMetadata<-function(path){
  
  cell.var <- rhdf5::h5readAttributes(file = path, name = 'obs')$`_index`
  cell_ids<-h5read(path, "/obs")[[cell.var]]
  
  cats_lvls<-h5read(path,'/obs/__categories')
  
  cats<-names(cats_lvls)
  
  nums<-setdiff(names(h5read(path, "/obs"))[!str_detect(names(h5read(path, "/obs")),'^__')],cats)
  
  message('cell ids format:',head(cell_ids))
  
  metanum <- do.call(cbind,lapply(nums, 
                                  function(n)h5read(path, paste0("/obs/",n))))
  colnames(metanum)<-nums
  
  
  cats.list<-lapply(names(cats_lvls),function(n)h5read(path, paste0("/obs/",n)))
  
  cats.list<-lapply(1:length(cats.list),function(i){
    if(is.list(cats.list[[i]])){
      new.cat<-paste(names(cats.list[i]),names(cats.list[[i]])[1])
      cats[i]<-new.cat
      return(cats_lvls[[i]][[1]][cats.list[[i]][[1]]+1])
    }else{
      return(cats_lvls[[i]][cats.list[[i]]+1])
    }
  })
  
  names(cats.list)<-cats
  metacat <- do.call(cbind,cats.list)
  colnames(metacat)<-cats
  
 
  metadata<-data.frame(as.data.frame(cbind(metacat,metanum)),row.names=cell_ids)
  message('head metadata:')
  
  print(head(metadata))
  return(metadata)
  
}

## Create Seurat Object for DLPFC SEA-AD single cell data
path='/projectnb/tcwlab-load/ref-data/SEAAD/DLPFC/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad'

#Get the count matrix
h5ls(path) #to check how look like the data
#we want the raw counts matrix which is store in `layers/UMIs`
data <- open_matrix_anndata_hdf5(path,group = 'layers/UMIs')
#convert the matrix to a BP directory format to store the data efficiently
write_matrix_dir(
  mat = data,
  dir = str_replace(path,".h5ad", "_BP")
)
#then open this BP matrix
mat <- open_matrix_dir(dir = str_replace(path,".h5ad", "_BP"))


# Get metadata from the h5ad file

metadata<-GetH5ADCellMetadata(path)
table(metadata$Subclass)
#create Seurat Object
seead_dlpfc <- CreateSeuratObject(counts = mat, meta.data = metadata)

#Save the seurat object
?SaveSeuratRds
SaveSeuratRds(
  object = seead_dlpfc,
  file = str_replace(path,".h5ad", ".rds")
)

#'sketch' a subset of cells
#put in cache representative cells
?SketchData

seead_dlpfc <- NormalizeData(seead_dlpfc)
seead_dlpfc <- FindVariableFeatures(seead_dlpfc)
seead_dlpfc <- SketchData(
  object = seead_dlpfc,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
seead_dlpfc
# switch to analyzing the full dataset (on-disk)
DefaultAssay(seead_dlpfc) <- "RNA"
# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(seead_dlpfc) <- "sketch"
#the following analysis is the same than classicaly, you have a 'sketch' data that can be analyze
#the results clustering  can then be projects in the full dataset
#for more see : https://satijalab.org/seurat/articles/seurat5_sketch_analysis 


## Do the same for the others region [MTG]
path='/projectnb/tcwlab-load/ref-data/SEAAD/MTG/SEAAD_MTG_RNAseq_final-nuclei.2023-07-19.h5ad'

#Get the count matrix
data <- open_matrix_anndata_hdf5(path,group = 'layers/UMIs')

#convert the matrix to a BP directory format to store the data efficiently
write_matrix_dir(
  mat = data,
  dir = str_replace(path,".h5ad", "_BP")
)
#then open this BP matrix
mat <- open_matrix_dir(dir = str_replace(path,".h5ad", "_BP"))


# Get metadata from the h5ad file

metadata<-GetH5ADCellMetadata(path)
table(metadata$Subclass)
#create Seurat Object
seead_mtg <- CreateSeuratObject(counts = mat, meta.data = metadata)
seead_mtg
unique(seead_mtg[[]]$Subclass)
#Save the seurat object
SaveSeuratRds(
  object = seead_mtg,
  file = str_replace(path,".h5ad", ".rds")
)



#we could have merge both in one same dataset simply by giving the list of matrices to CreateSeuratObject() function
#see https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette#load-data-from-multiple-h5ad-files

