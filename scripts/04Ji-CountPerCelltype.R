out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(future)
n_cores_mc=6
 n_core_future=4
# 
options(future.globals.maxSize = 10000 * 1024^2) #10GB



rds_files<-list.files('outputs/04-ROSMAP_MIT_ATAC/',pattern = 'Astro|Exc|Inh|VLMC|Mic|Endo|Oligo|OPC',full.names = TRUE)

mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')


peaks<-readRDS(fp(out,'perDonorGroups_celltype_peaks.rds'))


mtd_anno<-rbindlist(mclapply(rds_files,function(file){
  ct<-str_remove(basename(file),pattern = '.rds$')
  message('Counting for', ct)
  
  celltype<-readRDS(file)
  
 plan("multicore", workers = n_core_future)
  
 celltype<-AddMetaData(celltype,data.frame(mtd[main_cell_type==ct],row.names = 'cell_id_long'))
 
 peaks_mat<-FeatureMatrix(Fragments(celltype),
                          features =peaks,
                          process_n = 2000,
                          cells=colnames(celltype))
 
 celltype[["peaksDCT"]]<-CreateChromatinAssay(peaks_mat)
 DefaultAssay(celltype)<-'peaksDCT'

 #compute pct_read_in_peak
 celltype$pct_frag_in_peaksDCT<- celltype$nCount_peaksDCT/celltype$n_tot_fragments
 saveRDS(celltype,file)
 
 return(data.table(celltype@meta.data,keep.rownames = 'cell_id_long'))
 
                    
},mc.cores =n_cores_mc ))

mtd<-merge(mtd,mtd_anno[,.(cell_id,libraryID,nCount_peaksDCT,nFeature_peaksDCT,pct_frag_in_peaksDCT,cell_id_long)])

fwrite(mtd,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')

