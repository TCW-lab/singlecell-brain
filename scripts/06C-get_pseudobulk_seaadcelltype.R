
library(Seurat)
library(dplyr)
library(DESeq2)
library(stringr)
library(data.table)
source('../../utils/r_utils.R')

out<-'outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_celltype/'
dir.create(out)
#Load additional metadata
cellDF=read.csv('/projectnb/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Cell_Metadata_ByCell.csv')
rownames(cellDF)=cellDF$cellID
sampleDF=fread('outputs/06-ROSMAP_MIT_snRNA/integrated_clinical_metadata.csv')

cellRemoveBool=read.csv('/projectnb/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Outlier_Matrix.csv',row.names=1)

#Get outliers
finalRemove=apply(cellRemoveBool,1,any)

for(file in list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics',full.names = T,pattern = '.rds')){
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
    Idents(subraw)='predicted.id'
    mainct<-unique(str_extract(unique(subraw$cell_type_high_resolution),'Exc|Ast|Inh|Mic|Fib|^T|CAMs|SMC|OPC|Oli|Per|End'))
    #For each  cell type
    for(hres in unique(Idents(subraw))){
      hresSubraw=subset(subraw,idents = hres)
      
      if(str_detect(hres,paste(mainct,collapse = '|'))&ncol(hresSubraw)>100){
        message('generating pseudobulk for ',hres)
        #Pull out cells with high-res cell type
        bulk=AggregateExpression(hresSubraw, return.seurat = TRUE, slot = "counts", assays = "RNA", group.by = c("individualID"))
        
        #Append metadata to pseudobulk Seurat object
        bulkMD=bulk@meta.data
        #n cells per donor
        bulkMD$Cells<-table(hresSubraw$individualID)[rownames(bulkMD)]
        bulkMD$individualID=rownames(bulkMD)
        
        #add donor level info
        bulkMD=merge(bulkMD,sampleDF,by='individualID')
        head(bulkMD)
        rownames(bulkMD)=bulkMD$individualID
        bulk=AddMetaData(bulk,bulkMD)
        
        #Calculate %mitochondrial DNA, and attach to metadata
        bulk[['celltype.percent.mt']]=PercentageFeatureSet(bulk,pattern = '^MT-')
        hresFileName=make.names(gsub('/','-',hres))
        saveRDS(bulk,fp(out,paste0(hresFileName,'_Pseudobulk.rds')))
        #Memory management: Delete Seurat objects and free unused memory
        remove(hresSubraw)
      }else{
        message('not enough or not the good cells, skipping ', hres)
      }
      
    
      }
    
    #Memory management: Delete Seurat objects and free unused memory
    remove(raw)
    remove(subraw)
    gc()
  }
  

}

