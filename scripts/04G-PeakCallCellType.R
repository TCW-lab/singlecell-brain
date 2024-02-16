out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

#merge sample object per main cell type, add metadata

brain<-readRDS('outputs/04-ROSMAP_MIT_ATAC/brain_12best_samples_qc.rds')
DefaultAssay(brain)<-'peaks'
?CallPeaks
peaks<-CallPeaks(brain,group.by='cell_type',
                 macs2.path = '/projectnb/tcwlab/LabMember/adpelle1/micromamba/envs/macs2/bin/macs2')

saveRDS(peaks,fp(out,'brain_12best_samples_qc_celltype_peaks.rds'))


