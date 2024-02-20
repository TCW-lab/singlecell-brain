out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

#merge sample object per main cell type, add metadata

mtd<-fread(fp(out,'all_final_ATACseq_nuclei_metadata.csv.gz'))

for(ct in unique(mtd$main_cell_type)){
  message(ct)
  file_out<-fp(out,ps(ct,'.rds'))
  if(!file.exists(file_out)){
    brainct_list<-lapply(list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',
                                   full.names = T,recursive = F), function(d){
                                     f<-fp(d,'signac_object_annotated_celltype_peaks.rds')
                                     brainct<-SplitObject(readRDS(f),split.by='main_cell_type')[[ct]]
                                     return(brainct)
                                   })
    brainct_list<-brainct_list[!sapply(brainct_list,is.null)]
    brainct<-merge(brainct_list[[1]],brainct_list[2:length(brainct_list)],
                   merge.dr = c('ref.lsi','ref.umap'),merge.data=FALSE)
    saveRDS(brainct,file_out)

  }


}


# for(ct in unique(mtd$main_cell_type)){
#   message(ct)
#   file_out<-fp(out,ps(ct,'.rds'))
#   brainct<-readRDS(file_out)
#   
# 
#     brainct_list<-lapply(list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',
#                                    full.names = T,recursive = F), function(d){
#                                      f<-fp(d,'signac_object_annotated_celltype_peaks.rds')
#                                      brainct<-SplitObject(readRDS(f),split.by='main_cell_type')[[ct]]
#                                      return(brainct)
#                                    })
#     
#     brainct_list<-brainct_list[!sapply(brainct_list,is.null)]
#     brainct.old<-merge(brainct_list[[1]],brainct_list[2:length(brainct_list)], 
#                    merge.dr = c('ref.lsi','ref.umap'),merge.data=FALSE)
#     brainct.old<-brainct.old[,colnames(brainct)]
#     Fragments(brainct)<-Fragments(brainct.old)
#     
#     saveRDS(brainct,file_out)
#     
#   
# }

