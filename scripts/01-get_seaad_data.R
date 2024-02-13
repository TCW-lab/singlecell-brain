
#get all SEAAD data
source('../../utils/r_utils.R')
library(Seurat)
out<-'../singlecell-brain/outputs/01-SEAAD_data'
dir.create(out,recursive = T)


#get Full DLPFC and MTG  data####

#try to get from h5ad entire dataset
#try first see if Seurat conversion works on endothelial
system('curl -o /projectnb/tcwlab-load/ref-data/SEAAD/DLPFC/endothelial.h5ad "https://corpora-data-prod.s3.amazonaws.com/013c6a00-0507-49e8-998a-5022bd65621a/local.h5ad?AWSAccessKeyId=ASIATLYQ5N5XYHYGK57I&Signature=mNweGqhUiE17AFZ5A2PQ3wNBNG0%3D&x-amz-security-token=IQoJb3JpZ2luX2VjENv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJIMEYCIQDE7wdbmVBXMuaOzX7p3sWEDTSgFgsN6UgU0Mubo2U8KQIhAJgFjOljwKQ%2FmC5Fk%2Ft0wCaQhRgnuHyzLFB3JfseciDhKusDCBQQARoMMjMxNDI2ODQ2NTc1IgxPq01YYAAvROa4dn4qyANQtypMxD23b8f44aj1OWk2gsCmNxRuurSwXZqdERVvPFcMQmaKe8fKFo7rWZxGHTb04m6IEKblfCmLtxzTY%2Fv4OsuIETbyJczZVsrbL1XSkUlj%2BODM1WjIYvBve9YV7QpgDOZ6H144vUA%2FhGTc1%2Br%2B1gy0IwqWevi0JsGVM7GQeP%2BRJgroPDAHmqNPHExkpKCCxXT4IiVB0yCBhGzJiD2HOuNrxtFfyUtkjeUxyTmUmDVxm%2BOZBXi7j1hu%2FjUBjhKyqgB2Gilwn%2BcgUxIog0NAGxOTKaV5RgIdOhVAyPQq5g8Bn4h8m%2BBvqe1NEvChd3w3aS1cTtDDs1U%2B2lDHupw%2Bi%2FATnndjTm1fiWVIMlywTsTJJDySpyo%2FJLtEAFdjQucxywXsm8IR%2FrTQ18m1XzU4AlQM1Cmac4ARR25rovrzwk%2BuH2xrwAVhRAcbQsO%2Bma2DlHV4HtXik3WMEH2Xj%2BCec9IPs4USZQld8i27RN9Lb0w9krnl64Cpavj%2FBGGakVl%2FNWE6aqfu%2BrnsCELcGYem8A6J0C6zB5fnjMwC6APtzxdHxE02BPoQ4%2F5iepcurf2Azv6kD2SZd3Ee7uSXtPVBJa9LQSXsRrUwmN2IqgY6pAF2T4hIGViTGPKloDpjIzMMkp20to3X4trvis1CVYsLH7s8NrFpnY4YW8lH94UryaOhT9%2FDmBK0GeSsValTRSwFOYSoJf86dERjp8w23pxfp2EIodvy83gkqnFvjD0lNE3T8Z9416ejxzDJyNDZbNiEU4eCXRs1tyiLRAy28NSZwGXRSsoU8QHkSmAqIJBxlBbRiwMishc4VVTaqctUREVUn28Oxg%3D%3D&Expires=1699456777"')


#strat: get the full h5ad on aws, create one h5ad by cell type, and transform to Seurat.

#1) get MTG and DLPFC on AWS


cmds<-list(DLPFC=paste('wget -q -O','/usr2/postdoc/adpelle1/tcwlab-load/ref-data/SEAAD/DLPFC/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad',
                       'https://sea-ad-single-cell-profiling.s3.amazonaws.com/DLPFC/RNAseq/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad'),
           MTG=paste('wget -q -O','/usr2/postdoc/adpelle1/tcwlab-load/ref-data/SEAAD/MTG/SEAAD_MTG_RNAseq_final-nuclei.2023-07-19.h5ad',
                       'https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2023-05-05.h5ad'))


CreateJobFile(cmds,file = 'scripts/01i-get_full_data.qsub',nThreads = 4,maxHours = 12)
RunQsub('scripts/01i-get_full_data.qsub',job_name = 'getFullSEAAD')

#2) create one anndata (h5ad format) object by cell type
#run 01A
CreateJobForPyFile('scripts/01A-split_per_cell_type.py',micromamba_env = 'singlecell',
                   nThreads = 28,maxHours = 12)
RunQsub('scripts/01A-split_per_cell_type.qsub',job_name = 'splitSEAAD')

#extract the metdata 
CreateJobForPyFile('scripts/01Ai-extract_metadata.py',micromamba_env = 'singlecell',
                   nThreads = 16,maxHours = 12)
RunQsub('scripts/01Ai-extract_metadata.py',job_name = 'extractMTD')


#TBecause o trouble to use the SeuratV5 on disk option to manage large dataset (see 01C) we will downsampled and use the SeuratV4 
#to have a look at the data
#run 01B-
CreateJobForRfile('scripts/01B-downsampled_to_Seurat.R',nThreads = 16,maxHours = 24)
RunQsub('scripts/01B-downsampled_to_Seurat.R',job_name = 'ToSeurat',wait_for = 'extractMTD')

brain_mtg<-readRDS(file = file.path(out,'SEAAD_MTG_43k.rds'))

brain_mtg<-NormalizeData(brain_mtg)
brain_mtg<-FindVariableFeatures(brain_mtg)
brain_mtg<-ScaleData(brain_mtg)
brain_mtg<-RunPCA(brain_mtg)
brain_mtg<-RunUMAP(brain_mtg,dims = 1:50)
DimPlot(brain_mtg,group.by = 'Subclass',label = T)
DimPlot(brain_mtg,group.by = 'Class',label = T)

brain_mtg[['cell_type']]<-sapply(1:ncol(brain_mtg),function(i)ifelse(str_detect(brain_mtg$Class[i],'Glut'),paste0('Exc_',brain_mtg$Subclass[i]),
                                                                     ifelse(str_detect(brain_mtg$Class[i],'GABA'),paste0('Inh_',brain_mtg$Subclass[i]),brain_mtg$Subclass[i])))
DimPlot(brain_mtg,group.by = 'cell_type',label = T,label.size = 3.5,repel = T)+NoLegend()

brain_mtg[['RNA']]<-CreateAssayObject(counts = brain_mtg@assays$RNA@counts)

saveRDS(brain_mtg,file = file.path(out,'SEAAD_MTG_43k.rds'))


brain_dlpfc<-readRDS(file = file.path(out,'SEAAD_DLPFC_43k.rds'))

brain_dlpfc<-NormalizeData(brain_dlpfc)
brain_dlpfc<-FindVariableFeatures(brain_dlpfc)
brain_dlpfc<-ScaleData(brain_dlpfc)
brain_dlpfc<-RunPCA(brain_dlpfc)
brain_dlpfc<-RunUMAP(brain_dlpfc,dims = 1:50)
DimPlot(brain_dlpfc,group.by = 'Subclass',label = T)
DimPlot(brain_dlpfc,group.by = 'Class',label = T)

brain_dlpfc[['cell_type']]<-sapply(1:ncol(brain_dlpfc),function(i)ifelse(str_detect(brain_dlpfc$Class[i],'Glut'),paste0('Exc_',brain_dlpfc$Subclass[i]),
                                                                     ifelse(str_detect(brain_dlpfc$Class[i],'GABA'),paste0('Inh_',brain_dlpfc$Subclass[i]),brain_dlpfc$Subclass[i])))
DimPlot(brain_dlpfc,group.by = 'cell_type',label = T,label.size = 3.5,repel = T)+NoLegend()

brain_dlpfc[['RNA']]<-CreateAssayObject(counts = brain_dlpfc@assays$RNA@counts)

saveRDS(brain_dlpfc,file = file.path(out,'SEAAD_DLPFC_43k.rds'))



#DATA exploration and QCs
#see 01D




