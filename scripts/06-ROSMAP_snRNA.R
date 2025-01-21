source('../../utils/r_utils.R')
out<-'outputs/06-ROSMAP_MIT_snRNA/'
dir.create(out)
library(Seurat)

#ROSMAT_MITout#ROSMAT_MIT preprocessing####
obj<-readRDS('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Astrocytes.rds')
head(obj[[]])
table(obj$projid)
DimPlot(obj)

DimPlot(obj,group.by = 'projid')+NoLegend()

#metadata
mtd<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_assay_snRNAseq_metadata.csv')
mtds<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv')
mtdc<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Metadata/ROSMAP_clinical.csv')
mtdc2<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/ROSMAP_xqtl_complete_samples_covariates_sex_death_pmi_study.csv')
mtds2<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/ROSMAP_biospecimen_metadata.csv')
setnames(mtdc2,'sample_id','specimenID')
mtdc2<-merge(mtdc2,mtds2[,.(specimenID,individualID)],by='specimenID')



mtd<-merge(mtd[,-'assay'][libraryPreparationMethod=='10x'],mtds[assay=='scrnaSeq'],by='specimenID')
mtd<-merge(mtd,mtdc,by='individualID',all.x=T)

mtd<-merge(mtd,mtdc2,by='individualID',all.x=T)

mtd[,age_death:=ifelse(is.na(age_death.y),age_death.x,age_death.y)]
mtd[,msex:=ifelse(is.na(msex.y),msex.x,msex.y)]
mtd[,pmi:=ifelse(is.na(pmi.y),pmi.x,pmi.y)]
setnames(mtd,'specimenID.y','WGS_ID')
setnames(mtd,'specimenID.x','specimenID')
mtd<-unique(mtd,by='individualID')
fwrite(mtd,fp(out,'integrated_clinical_metadata.csv'))
mtd<-fread(fp(out,'integrated_clinical_metadata.csv'))


#downsampling
files_paths=list.files('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/',pattern = '.rds',full.names = T)

#mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column

pct_keep=0.2
cell_group='cell_type_high_resolution'
min_by_group=200
max_by_group=2000
#h5ls(files_paths[1]) #to know where the raw counts matrix is stored
#raw_counts_location<-'layers/UMIs'

brain_downsampled<-Reduce(function(x,y)merge(x,y),lapply(files_paths,function(file){
  message('reading ',file)
 # data <- open_matrix_anndata_hdf5(file,group = raw_counts_location) #check that UMIs matrix is stored in layers/UMIs using h5ls(path) or reading the h5ad file on python using anndata
 # mtdf<-mtd[colnames(data),on='cell_id']

  
 # dataf<-as(data[,mtdf[(to_keep)]$cell_id],'dgCMatrix')
  
  
 # objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[(to_keep)],row.names = 'cell_id'))
  obj<-readRDS(file)
  obj<-UpdateSeuratObject(obj)
  mtc<-merge(data.table(obj@meta.data,keep.rownames = 'cell_id'),mtd,by='projid')

  
  obj<-AddMetaData(object = obj,data.frame(mtc,row.names = 'cell_id'))
  
  mtc[,to_keep:=cell_id%in%head(sample(cell_id,ifelse(.N*pct_keep>min_by_group,round(.N*pct_keep),min_by_group)),max_by_group),
      by=cell_group]
  
  objf<-obj[,mtc[(to_keep)]$cell_id]
  return(objf)
  
}))
brain_downsampled
brain_downsampled<-NormalizeData(brain_downsampled)
brain_downsampled<-FindVariableFeatures(brain_downsampled)
brain_downsampled<-ScaleData(brain_downsampled)
brain_downsampled<-RunPCA(brain_downsampled)
brain_downsampled<-RunUMAP(brain_downsampled,dims = 1:50)
DimPlot(brain_downsampled,group.by = 'cell_type_high_resolution',label = T)+NoLegend()

saveRDS(brain_downsampled,file = file.path(out,'ROSMAP_MIT_DLPFC_downsampled_representative.rds'))

table(brain_downsampled$cell_type_high_resolution)

brain_downsampled<-readRDS(file = file.path(out,'ROSMAP_MIT_DLPFC_downsampled_representative.rds'))

DimPlot(brain_downsampled,group.by = 'cell_type',label = T)

#Pseudobulk by main cell type
out1<-fp(out,'Pseudobulk_main_celltype')
dir.create(out1)
#for1
obj<-readRDS('~/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Seurat_Objects/Astrocytes_Pseudobulk.rds')
obj<-CreateSeuratObject(obj@assays$RNA@counts,meta.data = obj@meta.data)
obj=RenameCells(obj,new.names = obj$x)
head(obj[[]])
mtd<-data.table(obj@meta.data)

mtdf<-mtd[,.(x,Cells,Reads.Cell, Features.Cell, percent.mt,Astrocytes,Immune_cells,Inhibitory_neurons,
       Oligodendrocytes,OPCs,Vasculature_cells,Excitatory_neurons,projid)]
mtdg<-mtd[,.SD,.SDcols = c('x',paste0('PC',1:20))]
mtd<-merge(mtdf,mtdg)
mtd[,individualID:=x]
obj@meta.data<-data.frame(mtd,row.names = 'x')[colnames(obj),]

saveRDS(obj,fp(out1,'Astrocytes_Pseudobulk.rds'))

#for all
pseudo_files<-list.files('~/tcwlab/LabMember/acandib/APOE_Ab/Pseudobulk_QC/Seurat_Objects/',pattern = '.rds',full.names = T)
pseudo_files<-pseudo_files[!str_detect(basename(pseudo_files),'set[1-3]')]
for(f in pseudo_files){
  message(basename(f))
  obj<-readRDS(f)
  obj<-CreateSeuratObject(obj@assays$RNA@counts,meta.data = obj@meta.data)
  obj=RenameCells(obj,new.names = obj$x)
  head(obj[[]])
  mtd<-data.table(obj@meta.data)
  
  mtdf<-mtd[,.(x,Cells,Reads.Cell, Features.Cell, percent.mt,Astrocytes,Immune_cells,Inhibitory_neurons,
               Oligodendrocytes,OPCs,Vasculature_cells,Excitatory_neurons,projid)]
  mtdg<-mtd[,.SD,.SDcols = c('x',paste0('PC',1:20))]
  mtd<-merge(mtdf,mtdg)
  mtd[,individualID:=x]
  obj@meta.data<-data.frame(mtd,row.names = 'x')[colnames(obj),]
  
  saveRDS(obj,fp(out1,basename(f)))
  
}

readRDS('outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_main_celltype/Astrocytes_Pseudobulk.rds')@meta.data
readRDS('outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_highRes_QC/Seurat_Objects/Ast.CHI3L1_Pseudobulk.rds')@meta.data

#true prop
pseudo_files<-list.files('outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_main_celltype/',pattern = '.rds',full.names = T)

props<-rbindlist(lapply(pseudo_files, function(f){
  ct=basename(f)|>str_remove('_Pseudobulk.rds')
  obj=readRDS(f)
  data.table(individualID=colnames(obj),
             n.cells=obj@meta.data[[ct]],
             cell_type=ct)
  
}),fill = T)
props[,n.cells.sample:=sum(n.cells),by='individualID']
props[,pct.cells:=n.cells/n.cells.sample]
ggplot(props)+geom_boxplot(aes(y=cell_type,x=pct.cells))+theme_bw()

fwrite(props,'outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_main_celltype/true_proportions.csv.gz')



#by subtypes
#get suntype pseudobulk
CreateJobForRfile('/projectnb/tcwlab/LabMember/adpelle1/projects/singlecell-brain/scripts/06A-get_pseudobulk_highres.R',nThreads = 16)
RunQsub('/projectnb/tcwlab/LabMember/adpelle1/projects/singlecell-brain/scripts/06A-get_pseudobulk_highres.R',job_name = 'highressnMIT')

#true prop
pseudo_files<-list.files('outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_highRes_QC/Seurat_Objects/',pattern = '.rds',full.names = T)

props<-rbindlist(lapply(pseudo_files, function(f){
  obj=readRDS(f)
  data.table(individualID=colnames(obj),
             n.cells=obj$Cells,
             cell_type=basename(f)|>str_remove('_Pseudobulk.rds'))
  
}),fill = T)
props[,n.cells.sample:=sum(n.cells),by='individualID']
props[,pct.cells:=n.cells/n.cells.sample]
ggplot(props)+geom_boxplot(aes(y=cell_type,x=pct.cells))+theme_bw()

fwrite(props,'outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_highRes_QC/true_proportions.csv.gz')

# annotate with SEA-AD ref####
dlpfc_ref<-readRDS('outputs/01-SEAAD_data/DLPFC/all_celltype_qc_small_50k.rds')
dlpfc_ref
dlpfc_ref <- NormalizeData(dlpfc_ref)
dlpfc_ref <- FindVariableFeatures(dlpfc_ref)
dlpfc_ref <- ScaleData(dlpfc_ref)
dlpfc_ref <- RunPCA(dlpfc_ref)
dlpfc_ref <- FindNeighbors(dlpfc_ref, dims = 1:50)
dlpfc_ref <- FindClusters(dlpfc_ref)
dlpfc_ref <- RunUMAP(dlpfc_ref,dims = 1:50)

head(dlpfc_ref[[]])
table(dlpfc_ref$cell_type)
DimPlot(dlpfc_ref,group.by='Subclass',label=T)+NoLegend()

dlpfc=readRDS(file = file.path(out,'ROSMAP_MIT_DLPFC_downsampled_representative.rds'))

dlpfc.anchors <- FindTransferAnchors(reference = dlpfc_ref, query = dlpfc, dims = 1:50,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = dlpfc.anchors, refdata = dlpfc_ref$cell_type, dims = 1:50)
dlpfc <- AddMetaData(dlpfc, metadata = predictions)
dlpfc[[]]|>head()
u1<-DimPlot(dlpfc,group.by='predicted.id',label=T)+NoLegend()
u2<-DimPlot(dlpfc,group.by='cell_type_high_resolution',label=T)+
  NoLegend()

u1+u2

#to get unified celltype, for Exc/Inh keep the SEA-AD anno, for others remove the cell 'state' anno (with a gene marker)

dlpfc$cell_type<-ifelse(str_detect(dlpfc$predicted.id,'^Exc|^Inh'),dlpfc$predicted.id,
                        str_remove(dlpfc$cell_type_high_resolution,' [A-Z0-9]+$'))
table(dlpfc$cell_type)
summary(factor(dlpfc$cell_type))

DimPlot(dlpfc,group.by='cell_type',label = T)+
  NoLegend()

table(dlpfc$predicted.id)
table(dlpfc$cell_type)

table(dlpfc$cell_type,dlpfc$cogdx)

dlpfc[[]]
dlpfc[['percent.mt']]<-PercentageFeatureSet(dlpfc,pattern = '^MT-')
FeaturePlot(dlpfc,features = 'nCount_RNA')
FeaturePlot(dlpfc,features = 'percent.mt')
saveRDS(dlpfc,file = file.path(out,'ROSMAP_MIT_DLPFC_downsampled_representative.rds'))


#check good 
mtd<-dlpfc@meta.data|>data.table()
mtd[cell_type|>str_detect('Lamp5')]$cell_type_high_resolution|>table()

u1<-DimPlot(dlpfc,group.by='cell_type',label=T)+NoLegend()
u2<-DimPlot(dlpfc,group.by='cell_type_high_resolution',label=T)+
  NoLegend()

u1+u2

#ok


#annotate others cells
#transfer label
#try by main celltype for 1
inh<-readRDS('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Inhibitory_neurons.rds')
dlpfc.anchors <- FindTransferAnchors(reference = dlpfc, query = inh, dims = 1:50,
                                     reference.reduction = "pca")
predictions <- TransferData(anchorset = dlpfc.anchors, refdata = dlpfc$cell_type, dims = 1:50)
inh <- AddMetaData(inh, metadata = predictions)
saveRDS(inh,'~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Inhibitory_neurons.rds')

inh[[]]|>head()
table(inh$predicted.id)
u1<-DimPlot(inh,group.by='cell_type_high_resolution',label=T)
u2<-DimPlot(inh,group.by='predicted.id',label=T)

u1+u2
ggplot(inh@meta.data)+geom_bar(aes(x=factor(projid),fill=cell_type_high_resolution),position = 'fill')
inh<-AddMetaData(inh,data.frame(merge(data.table(inh@meta.data,keep.rownames = 'cell_id'),
                           fread('outputs/06-ROSMAP_MIT_snRNA/integrated_clinical_metadata.csv'),by='projid'),row.names = 'cell_id'))

ggplot(inh@meta.data)+geom_bar(aes(x=factor(cogdx),fill=cell_type_high_resolution),position = 'dodge')
ggplot(inh@meta.data)+geom_bar(aes(x=factor(cogdx),fill=predicted.id),position = 'dodge')

table(inh$projid,inh$cell_type_high_resolution)

#for all
CreateJobForRfile('scripts/06B-transfer_seaad_celltype.R',nThreads = 28)
RunQsub()

#create psuedobulk by this celltype
CreateJobForRfile('scripts/06C-get_pseudobulk_seaadcelltype.R',nThreads = 28)
RunQsub(wait_for = 'transfer_seaad_')
#check
obj<-readRDS('outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_celltype/Ast_Pseudobulk.rds')

obj[['RNA']]<-as(obj[['RNA']],Class = 'Assay')

#for all
dir='outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_main_celltype/'
files<-list.files(dir,pattern = '.rds')

for(file in fp(dir,files)){
  obj<-readRDS(file)
  
  obj[['RNA']]<-as(obj[['RNA']],Class = 'Assay')
  saveRDS(obj,file)
}

dir='outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_celltype//'
files<-list.files(dir,pattern = '.rds')

for(file in fp(dir,files)){
  obj<-readRDS(file)
  
  obj[['RNA']]<-as(obj[['RNA']],Class = 'Assay')
  saveRDS(obj,file)
}

#true prop
pseudo_files<-list.files('outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_celltype/',pattern = '.rds',full.names = T)

props<-rbindlist(lapply(pseudo_files, function(f){
  obj=readRDS(f)
  data.table(individualID=colnames(obj),
             n.cells=obj$Cells,
             cell_type=basename(f)|>str_remove('_Pseudobulk.rds'))
  
}),fill = T)
props[,n.cells.sample:=sum(n.cells),by='individualID']
props[,pct.cells:=n.cells/n.cells.sample]
ggplot(props)+geom_boxplot(aes(y=cell_type,x=pct.cells))+theme_bw()

fwrite(props,'outputs/06-ROSMAP_MIT_snRNA/Pseudobulk_celltype/true_proportions.csv.gz')

#downsampling based on Seaad celltypes####
pct_keep=0.2
cell_group='predicted.id'
cell_id='cellID'

min_by_group=200
max_by_group=2000


mtd<-fread('../../../acandib/APOE_Ab/Pseudobulk_QC/Cell_Metadata_ByCell.csv')
#add pred id in the metadata
files_paths=list.files('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/',pattern = '.rds',full.names = T)

mtdid<-rbindlist(lapply(files_paths,function(f)readRDS(f)@meta.data|>data.table(keep.rownames =cell_id )))
mtdid
mtd<-merge(mtd,mtdid)
fwrite(mtd,fp(out,'Cell_Metadata_ByCell.csv.gz'))
mtd<-fread(fp(out,'Cell_Metadata_ByCell.csv.gz'))

#mtd[,cell_id:=exp_component_name] #here be sure to save the name of your cells matching with the name in the matrix in the cell_id column

#h5ls(files_paths[1]) #to know where the raw counts matrix is stored
#raw_counts_location<-'layers/UMIs'
# obj<-readRDS('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Astrocytes.rds')
# head(obj[[]])

mtd[,to_keep:=1:.N%in%head(sample(1:.N,ifelse(.N*pct_keep>min_by_group,
                                              round(.N*pct_keep),
                                              min(c(min_by_group,.N))),replace=F),
                           max_by_group),
    by=cell_group]

mtd[(to_keep)]

brain_downsampled<-Reduce(function(x,y)merge(x,y),
                          lapply(files_paths,function(file){
  message('reading ',file)
  # objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[(to_keep)],row.names = 'cell_id'))
  obj<-readRDS(file)
  obj<-UpdateSeuratObject(obj)
  mtc<-merge(data.table(obj@meta.data,keep.rownames = cell_id),
             mtd[,.SD,.SDcols=c(setdiff(colnames(mtd),
                                      colnames(obj@meta.data)),cell_id)],by=cell_id)
  
  obj<-AddMetaData(object = obj,data.frame(mtc,row.names = cell_id))
  
  objf<-obj[,mtc[(to_keep)][[cell_id]]]
  return(objf)
  
}))

brain_downsampled
brain_downsampled<-NormalizeData(brain_downsampled)
brain_downsampled<-FindVariableFeatures(brain_downsampled)
brain_downsampled<-ScaleData(brain_downsampled)
brain_downsampled<-RunPCA(brain_downsampled)
brain_downsampled<-RunUMAP(brain_downsampled,dims = 1:50)
DimPlot(brain_downsampled,group.by = 'predicted.id',label = T)+NoLegend()

saveRDS(brain_downsampled,file = file.path(out,'ROSMAP_MIT_DLPFC_downsampled_representative_seaad_anno.rds'))

table(brain_downsampled$predicted.id)

#downsampling based on main celltype####
pct_keep=0.2
cell_group='main_cell_type'
cell_id='cellID'

min_by_group=500
max_by_group=2000

mtd<-fread(fp(out,'Cell_Metadata_ByCell.csv.gz'))
#add main celltype in the metadata
mtd[,main_cell_type:=str_remove(cell_type,'_set[1-3]+')]

fwrite(mtd,fp(out,'Cell_Metadata_ByCell.csv.gz'))

mtd[,to_keep:=1:.N%in%head(sample(1:.N,ifelse(.N*pct_keep>min_by_group,
                                              round(.N*pct_keep),
                                              min(c(min_by_group,.N))),replace=F),
                           max_by_group),
    by=cell_group]

brain_downsampled<-Reduce(function(x,y)merge(x,y),lapply(files_paths,function(file){
  message('reading ',file)
  # objf<-CreateSeuratObject(dataf,meta.data = data.frame(mtdf[(to_keep)],row.names = 'cell_id'))
  obj<-readRDS(file)
  obj<-UpdateSeuratObject(obj)
  mtc<-merge(data.table(obj@meta.data,keep.rownames = cell_id),
             mtd[,.SD,.SDcols=c(setdiff(colnames(mtd),
                                        colnames(obj@meta.data)),cell_id)],
             by=cell_id)
  
  obj<-AddMetaData(object = obj,data.frame(mtc,row.names = cell_id))
  
  
  
  objf<-obj[,mtc[(to_keep)][[cell_id]]]
  return(objf)
  
}))

brain_downsampled
brain_downsampled<-NormalizeData(brain_downsampled)
brain_downsampled<-FindVariableFeatures(brain_downsampled)
brain_downsampled<-ScaleData(brain_downsampled)
brain_downsampled<-RunPCA(brain_downsampled)
brain_downsampled<-RunUMAP(brain_downsampled,dims = 1:50)
DimPlot(brain_downsampled,group.by = 'main_cell_type',label = T)+NoLegend()

saveRDS(brain_downsampled,
        file = file.path(out,
                         'ROSMAP_MIT_DLPFC_downsampled_representative_maincelltypes.rds'))

brain_downsampled<-readRDS(file = file.path(out,
                         'ROSMAP_MIT_DLPFC_downsampled_representative_maincelltypes.rds'))
table(brain_downsampled$main_cell_type)



#ANNEXE####
#check fgsea results
source('../../utils/r_utils.R')
res<-readRDS('../../../acandib/Project/APOE_Ab/Pseudobulk_QC/DESeq2_Results/fgsea_Results.RDS')

res_gsea<-res[["Astrocytes_Pseudobulk"]][["AD"]][["4-3"]][padj<0.05][order(padj)]


#annotate fgsea results
pathways_info<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
res_gsea<-merge(res_gsea,pathways_info,by=c('pathway'))
table(res_gsea[padj<0.05])

#table
ggplot(res_gsea[padj<0.05])+geom_bar(aes(x=NES>0,fill=subcat))+
  labs(y='number of enriched terms (padj<0.05)')

#toppathways by cell type
source('../../utils/visualisation.R')
res_gsea_top30<-res_gsea[padj<0.05,.SD[head(order(padj),40)],
                             by=.(subcat)]

emmaplot(head(res_gsea_top30[order(pval)],100),label.size=2,cols=c('blue','white','red'))

#only down
emmaplot(head(res_gsea[padj<0.05][order(pval)][NES<0],100),label.size=2,cols=c('blue','white','red'))



#autophagy 
res_gsea[padj<0.05][str_detect(pathway,'AUTO')]
emmaplot(head(res_gsea_top30[cell_type=='Astrocyte'][order(pval)],100),label.size=2,cols=c('blue','white','red'))
#rep to eostradiol..because sex bias ?
table(mtd[ipsc_type=='population'&cell_type=='Astrocyte']$genotype,mtd[ipsc_type=='population'&cell_type=='Astrocyte']$gender)


#get suntype pseudobulk
CreateJobForRfile('scripts/06A-get_pseudobulk_highres.R',nThreads = 16)
RunQsub(job_name = 'highressnMIT')
