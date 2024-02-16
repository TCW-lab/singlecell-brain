out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(future)
# n_cores_mc=6
 n_core_future=4
# 
# options(future.globals.maxSize = 10000 * 1024^2) #10GB

#merge sample object per main cell type, add metadata

mtd<-fread(fp(out,'all_final_ATACseq_nuclei_metadata.csv.gz'))

sample_dirs<-list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',
          full.names = T,recursive = F)
peaks<-readRDS(fp(out,'brain_12best_samples_qc_celltype_peaks.rds'))

mtd_anno<-rbindlist(lapply(sample_dirs[1:2],function(d){
  s<-basename(d)
  f<-fp(d,'signac_object.rds')
  message(s)
  
 plan("multicore", workers = n_core_future)
  
 brain<-readRDS(f)
 brain<-AddMetaData(brain,data.frame(mtd[libraryID==s],row.names = 'cell_id'))
 
 brain$cell_id<-colnames(brain)
 brain<-RenameCells(brain,add.cell.id =s)
 
 peaks_mat<-FeatureMatrix(Fragments(brain),
                          features =peaks,
                          process_n = 2000,
                          cells=colnames(brain))
 
 brain[["peaksCT"]]<-CreateChromatinAssay(peaks_mat)
 brain[["peaks"]]<-NULL
 #compute pct_read_in_peak
 brain$pct_frag_in_peaksCT<- brain$nCount_peaksCT/brain$n_tot_fragments
 saveRDS(brain,fp(d,'signac_object_annotated_celltype_peaks.rds'))
 
 return(data.table(brain@meta.data,keep.rownames = 'cell_id_long'))
 
                    
}))

mtd<-merge(mtd,mtd_anno[,.(cell_id,libraryID,nCount_peaksCT,nFeature_peaksCT,cell_id_long)])

fwrite(mtd,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')

