
out<-'outputs/01-SEAAD_data'
dir.create(out,recursive = T)


library(Seurat)
library(SeuratData)
library(SeuratDisk)
source('../../utils/r_utils.R')

#Functions####
AnnDataToSeurat<-function(afile,meta.data){
  message("Convert from Scanpy to Seurat...")
  require(data.table)
  require(stringr)
  Convert(afile, dest = "h5seurat", overwrite = TRUE)
  obj <- LoadH5Seurat(str_replace(afile,'h5ad$','h5seurat'),meta.data=F,misc=F)
  obj<-AddMetaData(obj,meta.data[colnames(obj),])
  system(paste('rm',str_replace(afile,'h5ad$','h5seurat')))
  return(obj)
}


#ANALYSIS####
#DLFPC
out1<-fp(out,'DLPFC')
mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))

afiles<-fp(out1,list.files(out1,pattern = '.h5ad'))

lapply(afiles, function(afile){
  print(afile)
  cells<-AnnDataToSeurat(afile,meta.data = data.frame(mtd,row.names ='exp_component_name'))
  cells
  #UMIs to RNA
  cells[['RNA']]<-CreateAssayObject(counts = cells@assays$UMIs@counts)
  cells[['UMIs']]<-NULL
  print(head(cells[[]]))
  saveRDS(cells,afile,'h5ad$','.rds')
})

#MTG
out1<-fp(out,'MTG')
mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))

afiles<-fp(out1,list.files(out1,pattern = '.h5ad'))

lapply(afiles, function(afile){
  print(afile)
  cells<-AnnDataToSeurat(afile,meta.data = data.frame(mtd,row.names ='exp_component_name'))
  cells
  #UMIs to RNA
  cells[['RNA']]<-CreateAssayObject(counts = cells@assays$UMIs@counts)
  cells[['UMIs']]<-NULL
  print(head(cells[[]]))
  saveRDS(cells,afile,'h5ad$','.rds')
})


