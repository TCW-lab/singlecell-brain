source('../../utils/r_utils.R')
library(Seurat)
out<-'outputs/01-SEAAD_data/'

# library(SeuratWrappers)
# library(flexmix)


for(region in c('DLPFC','MTG')){
  message('Within ', region,' data...')
  out1<-fp(out,region)
  out2<-fp(out1,'QC_small')
  dir.create(out2)
  rds_files<-list.files(out1,pattern = '.rds$',full.names = TRUE)
  
  for(file in rds_files){
    ct<-str_remove(basename(file),pattern = '.rds$')
    message('cells QC of', ct)
    
    celltype<-readRDS(file)
    VlnPlot(celltype,features =  c('Fraction.mitochondrial.UMIs','Genes.detected','Number.of.UMIs'),combine = F)
    ggsave(fp(out2,ps(ct,'qc_cell_metrics.png')),width = 8,height = 6)
    
    boxres<-boxplot.stats(celltype$Number.of.UMIs,coef = 3)
    celltype$UMIs.outlier<-colnames(celltype)%in%names(boxres$out)
    VlnPlot(celltype,features =  'Number.of.UMIs',split.by = 'UMIs.outlier')
    ggsave(fp(out2,ps(ct,'qc_nUMIs.png')),width = 4,height = 6)
    
    VlnPlot(celltype,features =  'Number.of.UMIs',split.by = 'UMIs.outlier',group.by = 'Donor.ID',fill.by = 'outlier.cellprop')
    ggsave(fp(out2,ps(ct,'qc_nUMIs_per_donor.png')),width = 10,height = 6)
    
    boxres<-boxplot.stats(celltype$Genes.detected,coef = 3)
    celltype$Genes.outlier<-colnames(celltype)%in%names(boxres$out)
    VlnPlot(celltype,features =  'Genes.detected',split.by = 'Genes.outlier')
    ggsave(fp(out2,ps(ct,'qc_nGenes.png')),width = 4,height = 6)
    
    VlnPlot(celltype,features =  'Genes.detected',split.by = 'Genes.outlier',group.by = 'Donor.ID',fill.by = 'outlier.cellprop')
    ggsave(fp(out2,ps(ct,'qc_nGenes_per_donor.png')),width = 10,height = 6)
    
    #percent.mt
    #try use MIQC
    #note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
    #Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.
    #?filterCells
    # celltype_sce<-as.SingleCellExperiment(celltype)
    # celltype_sce@colData$Genes.detected
    # model <- mixtureModel(celltype_sce,
    #                       subsets_mito_percent = 'Fraction.mitochondrial.UMIs',detected = 'Genes.detected')
    
    # model <- flexmix(Fraction.mitochondrial.UMIs ~ Genes.detected, data = celltype@meta.data, 
    #                  k = 2)
    # 
    # intercept1 <- parameters(model, component = 1)[1]
    # intercept2 <- parameters(model, component = 2)[1]
    # if (intercept1 > intercept2) {
    #   compromised_dist <- 1
    #   intact_dist <- 2
    # }else {
    #   intact_dist <- 1
    #   compromised_dist <- 2
    # }
    # 
    # celltype$Mito.outlier <- post[, compromised_dist] > 0.99
    # VlnPlot(celltype,features =  'Fraction.mitochondrial.UMIs',split.by = 'Mito.outlier')
    # FeatureScatter(celltype,'Fraction.mitochondrial.UMIs','Genes.detected',group.by ='Mito.outlier' )
    # #==> do not use MIQC for scNuc https://github.com/TCW-lab/SingleCell_APOE44/issues/2 
    
    boxres<-boxplot.stats(celltype$Fraction.mitochondrial.UMIs,coef = 3)
    celltype$Mito.high<-colnames(celltype)%in%names(boxres$out)
    VlnPlot(celltype,features =  'Fraction.mitochondrial.UMIs',split.by = 'Mito.high')
    ggsave(fp(out2,ps(ct,'qc_Mito.png')),width = 4,height = 6)
    
    VlnPlot(celltype,features =  'Fraction.mitochondrial.UMIs',split.by = 'Mito.high',group.by = 'Donor.ID',fill.by = 'outlier.cellprop')
    ggsave(fp(out2,ps(ct,'qc_Mito_per_donor.png')),width = 10,height = 6)
    
    #For scNuc, keep hardthreshold of 5% mito
    celltype$Mito.outlier<-celltype$Fraction.mitochondrial.UMIs>0.05
    
    #Final outlier is:
    celltype$donor.outlier<-celltype$outlier #conserved donor outlier metadata compute in previous step
    celltype$cell.outlier<-celltype$Genes.outlier|celltype$UMIs.outlier|celltype$Mito.outlier
    
    celltype$outlier<-celltype$donor.outlier|celltype$cell.outlier
    
    message(round(sum(celltype$donor.outlier)/nrow(celltype),digits = 1),'% cells flagged as donors outliers')
    message(round(sum(celltype$cell.outlier)/nrow(celltype),digits = 1),'% cells flagged as cells outliers' )
    
    message('In total ',round(sum(celltype$outlier)/nrow(celltype),digits = 1),'% cells flagged as outliers')
    
    #save the full object
    saveRDS(celltype,file)
    
    
  }
  
}
