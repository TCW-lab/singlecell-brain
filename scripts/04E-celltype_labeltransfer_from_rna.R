out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(future)

# check the current active plan

# Enable parallelization
plan("multicore", workers = 8)

options(future.globals.maxSize = 10000 * 1024^2) #10GB
plan()



brain<-readRDS(fp(out,'brain_12best_samples_qc.rds'))


#Integration with scRNA-seq data
# Load the pre-processed scRNA-seq data
brain_list<-lapply(list.files('outputs/01-SEAAD_data/DLPFC/QC_small/',pattern = '\\.rds$',full.names = TRUE), function(file){
  return(readRDS(file))
})
sum(sapply(brain_list,ncol)) #128384

brain_rna<-merge(brain_list[[1]],brain_list[2:length(brain_list)])
brain_rna<-NormalizeData(brain_rna)
brain_rna<-FindVariableFeatures(brain_rna,nfeatures = 5000)

transfer.anchors <- FindTransferAnchors(
  reference = brain_rna,
  query = brain,
  reduction = 'cca',
  dims = 1:40
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = brain_rna$cell_type,
  weight.reduction = brain[['lsi']],
  dims = 2:30
)

brain <- AddMetaData(object = brain, metadata = predicted.labels)
plot1 <- DimPlot(brain_rna, group.by = 'cell_type', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(brain, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
ps<-plot1 + plot2
ggsave(fp(out,'brain12_labeltransfer_with_seaad_sn_rna.png'),plot = ps,width = 9,height = 7)

# replace each label with its most likely prediction
for(i in levels(brain)) {
  cells_to_reid <- WhichCells(brain, idents = i)
  newid <- names(which.max(table(brain$predicted.id[cells_to_reid])))
  Idents(brain, cells = cells_to_reid) <- newid
}

brain$cell_type<-Idents(brain)

DimPlot(brain, group.by = 'cell_type', label = TRUE, repel = TRUE) + NoLegend() 
ggsave(fp(out,'brain12_cell_type_label.png'),width = 9,height = 7)

saveRDS(brain,fp(out,'brain_12best_samples_qc.rds'))


