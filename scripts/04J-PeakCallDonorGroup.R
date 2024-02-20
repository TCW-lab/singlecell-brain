
out<-'outputs/04-ROSMAP_MIT_ATAC/'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)


rds_files<-list.files('outputs/04-ROSMAP_MIT_ATAC/',pattern = 'Astro|Exc|Inh|VLMC|Mic|Endo|Oligo|OPC',full.names = TRUE)

mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')


#at celltype level, peak call for each groups of donors with similar chromatin accessibility profile

peaks_list<-Reduce(`c`,lapply(rds_files,function(file){

      ct<-str_remove(basename(file),pattern = '.rds$')
    message('Peak calls for', ct)
    
  celltype<-readRDS(file)
  
  #rm outlier
  celltype<-celltype[,!celltype$outlier]
 
  Idents(celltype)<-'cell_type'
  
  
  peaks_list<-lapply(unique(celltype$cell_type), function(ct){
    message(ct)
    celltypef<-subset(celltype,idents=ct)
    # 1) identify groups of donors with similar chromatin access profile
    #filter donors : at least 20 cells 
    samples_to_keep<-names(which(table(celltypef$libraryID)>20))
    celltypef<-celltypef[,celltypef$libraryID%in%samples_to_keep]
    
    celltype_pseudo<-AggregateExpression(celltypef,slot = 'data',
                                         assays = 'peaksCT',return.seurat = T,group.by = 'libraryID')
    
    celltype_pseudo <- FindTopFeatures(celltype_pseudo, min.cutoff = 'q0',assays = 'peaksCT')
    celltype_pseudo <- RunSVD(object = celltype_pseudo)
    
    
    celltype_pseudo <- FindNeighbors(
      object = celltype_pseudo,
      reduction = 'lsi',
      dims = 1:6
    )
    celltype_pseudo <- FindClusters(
      object = celltype_pseudo,
      algorithm = 3,
      resolution = 1.2,
      verbose = FALSE
    )
    
    #save donor cluster
    mtd[cell_type==ct,donor.ct.group:=celltype_pseudo@meta.data[libraryID,'seurat_clusters']]
    celltypef<-AddMetaData(celltypef,data.frame(mtd[cell_type==ct][,.(cell_id_long,donor.ct.group)],
                                                row.names = 'cell_id_long'))
    
    #2) peak call per groups
    

    peaks<-CallPeaks(celltypef,group.by=c('donor.ct.group'),
                     macs2.path = '/projectnb/tcwlab/LabMember/adpelle1/micromamba/envs/macs2/bin/macs2')
    peaks$peak_called_in<-paste0(ct,peaks$peak_called_in)

    return(peaks)
  })
  
  return(peaks_list)



}))
  
#2) combine peaks and reduce, maintaining ident information
peaks.combined <- Reduce(f = c, x = peaks_list)
peaks <- reduce(x = peaks.combined, with.revmap = TRUE)

peaks$peak_called_in <-sapply(1:length(peaks), function(i){
  datasets <-  peaks.combined$peak_called_in[peaks$revmap[[i]]]
  return(paste(unique(datasets), collapse = ";"))
})

peaks$revmap <- NULL

saveRDS(peaks,fp(out,'perDonorGroups_celltype_peaks.rds'))


