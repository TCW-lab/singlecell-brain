---
title: "Downsampling single Cell Data R Notebook"
output: html_notebook
editor_options: 
  chunk_output_type: inline
---

# Downsampling single Cell Data

``` r
library(data.table)
setDTthreads(0)
library(Seurat) 

library(BPCells)
library(rhdf5)
```

## From one Seurat object

Here we used the metadata of the object to randomly select 10% of the cells

``` r
pct_keep=0.1
mtd<-data.table(obj@meta.data,keep.rownames = 'cell_id')
mtd[,to_keep:=cell_id%in%sample(cell_id,round(.N*pct_keep))]
objf<-obj[,mtd[(to_keep)]$cell_id]
```

### Downsampling by cell type

If we have some rare cell population, we would like to preserve it in downsampled data. To do that we can add condition : For example take at least 500 cells by cell_type and max 2000 cells,

``` r
  pct_keep=0.2
  cell_group='cell_type'
  min_by_group=500
  max_by_group=2000
  #subset div by 10
  mtd<-data.table(obj@meta.data,keep.rownames = 'cell_id')
  mtd[,to_keep:=cell_id%in%head(sample(cell_id,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
      by=cell_group] 
  
  objf<-obj[,mtd[(to_keep)]$cell_id]
```

## Downsampling from several rds file

If you have several seurat objects, that need to be merge but downsampled before because of size issue, you can use this script

```r
files_paths=list.files('../ref-data/SEAAD/MTG/per_celltype/',pattern = '.rds',full.names = T)
pct_keep=0.2
cell_group='cell_type'
min_by_group=500
max_by_group=2000
dim_red_to_merge=c('scVI','umap')
  
brain_downsampled<-Reduce(function(x,y)merge(x,y,merge.dr=dim_red_to_merge),lapply(files_paths,function(file){
   obj<-readRDS(file)
  mtd<-data.table(obj@meta.data,keep.rownames = 'cell_id')
  mtd[,to_keep:=cell_id%in%head(sample(cell_id,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
      by=cell_group] 
  
  objf<-obj[,mtd[(to_keep)]$cell_id]
  
  return(objf)
  
}))

brain_downsampled
```

## Get downsampled seurat obejct from several h5ad file

If you have several h5ad objects, that need to be merge and transform into one seurat object but downsampled before because of size issue, you can use this script

```r
#make sure you are using Seurat V4
#if you already have v5 install, to get back to v4 version:
# remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

files_paths=list.files('outputs/01-SEAAD_data/DLPFC',pattern = '.h5ad',full.names = T)

mtd<-fread('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz')  #metadata containing the celltype/donors annotation. should be extract from the anndata object before, see LINKTOADD
mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column

pct_keep=0.2
cell_group='Subclass'
min_by_group=500
max_by_group=2000
h5ls(files_paths[1]) #to know where the raw counts matrix is stored

raw_counts_location<-'layers/UMIs'

brain_downsampled<-Reduce(function(x,y)merge(x,y),lapply(files_paths,function(file){
  data <- open_matrix_anndata_hdf5(file,group = raw_counts_location) #check that UMIs matrix is stored in layers/UMIs using h5ls(path) or reading the h5ad file on python using anndata
  mtdf<-mtd[colnames(data),on='cell_id']
  mtdf[,to_keep:=cell_id%in%head(sample(cell_id,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
      by=cell_group]

  dataf<-as(data[,mtdf[(to_keep)]$cell_id],'dgCMatrix')


  objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[(to_keep)],row.names = 'cell_id'))

  return(objf)
  
}))
brain_downsampled

```

## Use build in SeuratV5 function

If you do not have a specific cell annotation (cell type, donors..) to use to as cell groups to downsample the objects, you can use build in SeuratV5 function that downsampled object based on a gene expression profile score to conserve rare cell profile (i.e. cell type) warning, here we need to have first [created a Seurat V5 Object]('notebooks/01-convert_SEAAD_h5ad_to_SeuratV5.md')

``` r
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
