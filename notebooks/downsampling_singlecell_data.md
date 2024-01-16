---
title: "Downsampling single Cell Data R Notebook"
output: html_notebook
---


# Downsampling single Cell Data
```r
library(data.table)
library(Seurat)

```

## From one Seurat object

Here we used the metadata of the object to randomly select 10% of the cells

```r
pct_keep=0.1
mtd<-data.table(obj@meta.data,keep.rownames = 'bc')
mtd[,to_keep:=bc%in%sample(bc,round(.N*pct_keep))]
objf<-obj[,mtd[(to_keep)]$bc]
```

### Downsampling by cell type

If we have some rare cell population, we would like to preserve it in downsampled data. To do that we can add condition : For example take at least 500 cells by cell_type and max 2000 cells,

```r
  pct_keep=0.2
  cell_group='cell_type'
  min_by_group=500
  max_by_group=2000
  #subset div by 10
  mtd<-data.table(obj@meta.data,keep.rownames = 'bc')
  mtd[,to_keep:=bc%in%head(sample(bc,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
      by=cell_group] 
  
  objf<-obj[,mtd[(to_keep)]$bc]
```

## Downsampling from several rds file

If you have several seurat objects, that need to be merge but downsampled before because of size issue, you can use this script

```r
files_paths=list.files('../ref-data/SEAAD/MTG/per_celltype/',pattern = '.rds')
pct_keep=0.2
cell_group='cell_type'
min_by_group=500
max_by_group=2000
  
brain_downsampled<-Reduce(function(x,y)merge(x,y,merge.dr=c('scVI','umap')),lapply(files_paths,function(fp)
  obj<-readRDS(fp)
  mtd<-data.table(obj@meta.data,keep.rownames = 'bc')
  mtd[,to_keep:=bc%in%head(sample(bc,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
      by=cell_group] 
  
  objf<-obj[,mtd[(to_keep)]$bc]
  
  return(objf)
  
))
brain_downsampled
```

## Use build in SeuratV5 function

If you do not have a specific cell annotation (cell type, donors..) to use to as cell groups to downsample the objects, you can use build in SeuratV5 function that downsampled object based on a gene expression profile score to conserve rare cell profile (i.e. cell type) warning, here we need to have first [created a Seurat V5 Object]('notebooks/01-convert_SEAAD_h5ad_to_SeuratV5.md')

```r
rds_path='../ref-data/SEAAD/DLPFC/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.rds'
seead_dlpfc<-readRDS(rds_path)
DefaultAssay(seead_dlpfc)<-'RNA'
seead_dlpfc <- SketchData(
  object = seead_dlpfc,
  ncells = 50000,
  assay = 'RNA',
  method = "LeverageScore",
  sketched.assay = "sketch"
)
seead_dlpfc

# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(seead_dlpfc) <- "sketch"
```
FIXME : it seems that the SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.rds do not contain the count matrix stored correctly, we have to recreate this rds file starting from the BP matrix. 

for more information on Sketch analysis using Seurat, see [this seurat vignette](https://satijalab.org/seurat/articles/seurat5_sketch_analysis)
