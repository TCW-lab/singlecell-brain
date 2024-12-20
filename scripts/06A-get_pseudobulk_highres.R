
library(Seurat)
library(dplyr)
library(DESeq2)
library(stringr)
library(data.table)

#Load additional metadata
cellDF=read.csv('/projectnb/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Cell_Metadata_ByCell.csv')
rownames(cellDF)=cellDF$cellID
sampleDF=fread('/projectnb/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Cell_Metadata_BySample.csv')
#clean
cols_torm<-sapply(sampleDF, function(x)is.character(x)&any(str_detect(x,'\\,')))
sampleDF<-sampleDF[,.SD,.SDcols=!which(cols_torm)]
sampleDF<-sampleDF[,-'Cells']

cellRemoveBool=read.csv('/projectnb/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Outlier_Matrix.csv',row.names=1)

#Get outliers
finalRemove=apply(cellRemoveBool,1,any)

for(file in list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics',full.names = T)){
  if(substr(basename(file),
            start=nchar(basename(file))-3,stop=nchar(basename(file)))=='.rds'){
    print(basename(file))
    raw=readRDS(file)
    raw=Seurat::UpdateSeuratObject(raw)
    rawMDNames=colnames(raw)
    
    keepCells=!finalRemove[rawMDNames]
    subraw=raw[,keepCells]
    #Append additional metadata to scRNASeq object
    subraw@meta.data[['individualID']]=cellDF[rownames(subraw@meta.data),]$individualID
    Idents(subraw)='cell_type_high_resolution'
    #For each high-res cell type
    for(hres in unique(Idents(subraw))){
      print(hres)
      
      #Pull out cells with high-res cell type
      hresSubraw=subset(subraw,idents = hres)
      bulk=AggregateExpression(hresSubraw, return.seurat = TRUE, slot = "counts", assays = "RNA", group.by = c("individualID"))
      
      
      #Append metadata to pseudobulk Seurat object
      bulkMD=bulk@meta.data
      #n cells per donor
      bulkMD$Cells<-table(hresSubraw$individualID)[rownames(bulkMD)]
      
      bulkMD$individualID=rownames(bulkMD)
      
      bulkMD=merge(bulkMD,sampleDF,by='individualID')
      head(bulkMD)
      rownames(bulkMD)=bulkMD$individualID
      bulk=AddMetaData(bulk,bulkMD)
      
      #Calculate %mitochondrial DNA, and attach to metadata
      bulk[['celltype.percent.mt']]=PercentageFeatureSet(bulk,pattern = '^MT-')
      hresFileName=make.names(gsub('/','-',hres))
      saveRDS(bulk,paste0('/projectnb/tcwlab/LabMember/adpelle1/projects/singlecell-brain/outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_highRes_QC/Seurat_Objects/',hresFileName,'_Pseudobulk.rds'))
      #Memory management: Delete Seurat objects and free unused memory
      remove(hresSubraw)
    }
    
    #Memory management: Delete Seurat objects and free unused memory
    remove(raw)
    remove(subraw)
    gc()
    
  }
}