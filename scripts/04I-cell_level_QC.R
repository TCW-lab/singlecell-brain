
out<-'outputs/04-ROSMAP_MIT_ATAC/celltype_qc'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)


rds_files<-list.files('outputs/04-ROSMAP_MIT_ATAC/',pattern = 'Astro|Exc|Inh|VLMC|Mic|Endo|Oligo|OPC',full.names = TRUE)


mtd<-rbindlist(lapply(rds_files,function(file){
    ct<-str_remove(basename(file),pattern = '.rds$')
    message('cells QC of', ct)
    
  #3* IQR? of nCount, nFeature, pct read, tss enrichment etc
  celltype<-readRDS(file)
  head(celltype[[]])
dup_bc<-duplicated(colnames(celltype))
if(sum(dup_bc)>0){
  warning('removing ',sum(dup_bc),'duplicated cells')
  #remove duplicated cell names
  celltype<-celltype[,!duplicated(colnames(celltype))]
  
}


    table(celltype$disease)
  #unique(data.table(celltype@meta.data,keep.rownames = 'cell_id_long')[,.(cell_id_long,libraryID)],by='libraryID')
  VlnPlot(celltype,features =  c('pct_frag_in_peaksCT','nFeature_peaksCT','nCount_peaksCT'),pt.size = 0,split.by = 'disease')
  
  ggsave(fp(out,ps(ct,'qc_cell_metrics_set1.png')),width = 8,height = 6)
  
  VlnPlot(celltype,features =  c('nFeature_peaksCT','nCount_peaksCT'),pt.size = 0,split.by = 'disease',log=T)
  ggsave(fp(out,ps(ct,'qc_cell_metrics_set1_log.png')),width = 8,height = 6)
  
  VlnPlot(celltype,features =  c('TSS.enrichment','nucleosome_signal'),split.by = 'disease',pt.size = 0,log=T)
  ggsave(fp(out,ps(ct,'qc_cell_metrics_set2.png')),width = 8,height = 6)
  
  #ggsave(fp(out1,ps(ct,'qc_cell_metrics.png')),width = 8,height = 6)
  
  boxres<-boxplot.stats(celltype$nCount_peaksCT,coef = 3)
  celltype$nCount_peaksCT.high<-colnames(celltype)%in%names(boxres$out)
  sum(celltype$nCount_peaksCT.high) #74/6908 #2%
  #table(celltype$nCount_peaksCT.outlier,celltype$libraryID)#but a removed from one donor
  
  
  boxres<-boxplot.stats(celltype$nFeature_peaksCT,coef = 3)
  celltype$nFeature_peaksCT.high<-colnames(celltype)%in%names(boxres$out)
  sum(celltype$nFeature_peaksCT.high) #54/6908 #1%
  table(celltype$nFeature_peaksCT.high,celltype$libraryID)
  VlnPlot(celltype,features =  'nFeature_peaksCT',split.by = 'nFeature_peaksCT.high',pt.size = 0,log=T)
  ggsave(fp(out,ps(ct,'qc_nFeature_peaksCT.png')),width = 4,height = 6)
  #can exclude based on that
  celltype$nFeature_peaksCT.outlier<-celltype$nFeature_peaksCT.high
  
  
  
  VlnPlot(celltype,features =  'pct_frag_in_peaksCT',pt.size = 0,log=T)
  boxres<-boxplot.stats(celltype$pct_frag_in_peaksCT,coef = 3)
  celltype$pct_frag_in_peaksCT.high<-colnames(celltype)%in%names(boxres$out)
  sum(celltype$pct_frag_in_peaksCT.high) #0
  table(celltype$pct_frag_in_peaksCT.high,celltype$libraryID)
  
  VlnPlot(celltype,features =  'pct_frag_in_peaksCT',
          split.by = 'pct_frag_in_peaksCT.high',pt.size = 0,log=T)
  #can't exclude based on that
  #celltype$nFeature_peaksCT.outlier<-celltype$nFeature_peaksCT.high
  

  VlnPlot(celltype,features =  'TSS.enrichment',pt.size = 0,log=T)
  boxres<-boxplot.stats(celltype$TSS.enrichment,coef = 3)
  celltype$TSS.enrichment.high<-colnames(celltype)%in%names(boxres$out)
  sum(celltype$TSS.enrichment.high) #496
  table(celltype$TSS.enrichment.high,celltype$libraryID)
  VlnPlot(celltype,features =  'TSS.enrichment',split.by = 'TSS.enrichment.high',pt.size = 0,log=T)
  #can't exclude based on that because only loww tss enrichemnnt is bad
  #do normal outlier detection to detect low score outlier
  boxres<-boxplot.stats(celltype$TSS.enrichment)
  celltype$TSS.enrichment.low<-colnames(celltype)%in%names(boxres$out)&celltype$TSS.enrichment<boxres$stats[2]
  sum(celltype$TSS.enrichment.low) #1
  VlnPlot(celltype,features =  'TSS.enrichment',split.by = 'TSS.enrichment.low',pt.size = 0,log=T)
  ggsave(fp(out,ps(ct,'qc_TSS.enrichment.png')),width = 4,height = 6)
  
  celltype$TSS.enrichment.outlier<-celltype$TSS.enrichment.low
  
  
  
  VlnPlot(celltype,features =  'nucleosome_signal',pt.size = 0,log=T)
  boxres<-boxplot.stats(celltype$nucleosome_signal,coef = 3)
  celltype$nucleosome_signal.high<-colnames(celltype)%in%names(boxres$out)
  sum(celltype$nucleosome_signal.high) #99
  table(celltype$nucleosome_signal.high,celltype$libraryID) #outliers fairly distribute accross donors
  VlnPlot(celltype,features =  'nucleosome_signal',split.by = 'nucleosome_signal.high',pt.size = 0,log=T)
  ggsave(fp(out,ps(ct,'qc_nucleosome_signal.png')),width = 4,height = 6)
  
  #can exclude based on that
  celltype$nucleosome_signal.outlier<-celltype$nucleosome_signal.high
  
  
  #Final outlier is:
  celltype$cell.outlier<-celltype$nucleosome_signal.outlier|celltype$TSS.enrichment.outlier|celltype$nFeature_peaksCT.outlier
  
  
  celltype$outlier<-celltype$donor.outlier|celltype$cell.outlier
  
  message(round(sum(celltype$donor.outlier)/nrow(celltype)*100,digits = 2),'% nuclei flagged as donors outliers')
  message(round(sum(celltype$cell.outlier)/nrow(celltype)*100,digits = 2),'% nuclei flagged as cells outliers' )
  
  message('In total ',round(sum(celltype$outlier)/nrow(celltype)*100,digits = 2),'% nuclei flagged as outliers')
  
  #save the full object
  saveRDS(celltype,file)
  
  return(data.table(celltype@meta.data,keep.rownames = 'cell_id_long'))
}),fill = T)

fwrite(mtd,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')
