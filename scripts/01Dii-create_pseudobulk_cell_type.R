
#Get  seurat obejct for each cell type
out<-'outputs/01-SEAAD_data/'
library(data.table)
setDTthreads(0)
library(Seurat) 
fp<-function(...)file.path(...)
library(stringr)

#Analysis####
sample_column='Donor.ID'


# #DLPFC
# out1<-file.path(out,'DLPFC')
# 
# files_paths=list.files(out1,pattern = '.rds',full.names = T)
# 
# for(file in files_paths){
#   ct<-str_remove(basename(file),'\\.rds')
# 
# 
#   celltype<-readRDS(file)
#   pseudo_mat<-AggregateExpression(celltype,assays = 'RNA',slot = 'count',group.by = sample_column,
#                                   return.seurat = FALSE)
# 
#     fwrite(data.table(pseudo_mat$RNA,keep.rownames = 'gene_id'),file.path(out1,paste0(ct,'_pseudobulk.csv.gz')))
# 
# 
# }


#MTG
out1<-file.path(out,'MTG')

files_paths=list.files(out1,pattern = '.rds',full.names = T)

for(file in files_paths){
  ct<-str_remove(basename(file),'\\.rds')
  
  
  celltype<-readRDS(file)
  pseudo_mat<-AggregateExpression(celltype,assays = 'RNA',slot = 'count',group.by = sample_column,
                                  return.seurat = FALSE)
  
  fwrite(data.table(pseudo_mat$RNA,keep.rownames = 'gene_id'),file.path(out1,paste0(ct,'_pseudobulk.csv.gz')))
  
  
}
