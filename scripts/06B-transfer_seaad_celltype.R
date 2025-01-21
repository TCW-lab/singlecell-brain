
source('../../utils/r_utils.R')
library(Seurat)
options(future.globals.maxSize = 1e9)

out<-'outputs/06-ROSMAP_MIT_snRNA/'

dlpfc<-readRDS(file = file.path(out,'ROSMAP_MIT_DLPFC_downsampled_representative.rds'))

rds_files<-list.files('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/',
                      pattern = '.rds',full.names = T)
for(f in rds_files){
  message('reading ',basename(f))
  obj<-readRDS(f)
  if(!'predicted.id'%in%colnames(obj[[]])){
    obj=Seurat::UpdateSeuratObject(obj)
    
    dlpfc.anchors <- FindTransferAnchors(reference = dlpfc, query = obj,
                                         dims = 1:50,
                                         reference.reduction = "pca")
    predictions <- TransferData(anchorset = dlpfc.anchors, 
                                refdata = dlpfc$cell_type,
                                dims = 1:50)
    obj <- AddMetaData(obj, metadata = predictions)
    saveRDS(obj,f)
    
  }
 
  
}
