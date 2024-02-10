#get pseudobulk and metadata 
library(data.table)
library(Seurat)
out<-'outputs/03-pseudobulk_data'
dir.create(out)


#ROSMAP_Multiomics
out1<-file.path(out,'MIT_ROSMAP_Multiomics')
dir.create(out1)
astro<-readRDS('../../../acandib/Project/APOE_Ab/Pseudobulk/Astrocytes.rds_Pseudobulk.rds')



astro@meta.data<-data.frame(data.table(astro@meta.data,keep.rownames = 'sample_id')[,.(x,projid, Study, msex.x ,educ,
                                                           race, spanish, apoe_genotype,   age_at_visit_max,
                                                           age_first_ad_dx,age_death.x, cts_mmse30_first_ad_dx, cts_mmse30_lv, pmi.x ,braaksc, ceradsc ,cogdx, dcfdx_lv)],row.names = 'x')


saveRDS(astro,file.path(out1,'astro_pseudobulk_seurat.rds'))

dim(astro@assays$RNA@counts) #33538   427
pseudo_bulk_counts<-astro@assays$RNA@counts

head(astro@meta.data) #33538   427

metadata<-astro@meta.data

mtd<-data.table(metadata,keep.rownames = 'sample_id') #to put in data.table format

