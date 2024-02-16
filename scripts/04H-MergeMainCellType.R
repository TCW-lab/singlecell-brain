out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

#merge sample object per main cell type, add metadata

mtd<-fread(fp(out,'all_final_ATACseq_nuclei_metadata.csv.gz'))

for(ct in unique(mtd$main_cell_type)){
  meassage(ct)
  brainct_list<-lapply(list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',
                                 full.names = T,recursive = F)[1:2], function(d){
                                   s<-basename(d)
                                   f<-fp(d,'signac_object.rds')
                                   brain<-readRDS(f)
                                   brain<-SplitObject(brain,split.by='main_cell_type')[[ct]]
                                   return(brain)
                                 })
  
  brainct<-merge(brainct_list[[1]],brainct_list[2:length(brainct_list)], merge.dr = c('ref.lsi','ref.umap'),merge.data=FALSE)
  saveRDS(brainct,fp(out,ps(ct,'.rds')))
  
  
  
}
