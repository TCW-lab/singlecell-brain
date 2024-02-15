out<-'outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(parallel)
n_cores=6

signac_files<-file.path(list.dirs(out,full.names = T,recursive =F ),'signac_object.rds')

annotations<-readRDS('ref-data/EnsDb.Hsapiens.v86_annotations.rds')
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
mts<-fread('ref-data/MIT_ROSMAP_Multiomics/snATAC-seq/samples_metadata.csv')


n_fragments<-fread('outputs/04-ROSMAP_MIT_ATAC/n_fragment_per_barcode_alllibs.csv.gz')
setnames(n_fragments,'N','n_tot_fragments')


mtd_all<-rbindlist(mclapply(signac_files,function(f){
  s<-basename(dirname(f))
  message(s)
  
  brain<-readRDS(f)
  
  #i) annotate with genome annoation
  Annotation(brain) <- annotations
  
  
  
  #ii) add the clinical data / assay metadata
  if(!'libraryID'%in%colnames(brain[[]])){
    brain$libraryID<-brain$donor_id
    brain$donor_id<-NULL
  }

  mtd_cells<-merge(mts,
                   data.table(brain@meta.data[,c('nCount_peaks','nFeature_peaks','libraryID')],
                              keep.rownames = 'cell_id')[,.(cell_id,libraryID)],by='libraryID')
  
  brain<-AddMetaData(brain,data.frame(mtd_cells,row.names = 'cell_id'))
  head(brain[[]])
  
  # iii) add nuclei QC metric
  #add % reads falling in peak
  brain<-AddMetaData(brain,data.frame(n_fragments[libraryID==unique(brain$libraryID)],row.names = 'cell_id'))
  
  brain$pct_reads_in_peaks<-brain$nCount_peaks/brain$n_tot_fragments
  head(brain[[]])
  
  #add Nucleosome signal metrics, TSS enrichment
  brain <- NucleosomeSignal(object = brain)
  brain<-TSSEnrichment(brain)
  brain[[]]
  VlnPlot(brain,c('nCount_peaks','pct_reads_in_peaks','TSS.enrichment','nucleosome_signal'),ncol = 2)
  ggsave(fp(out,s,'nuclei_qc_metric.png'),width = 7,height = 8)
  
  saveRDS(brain,fp(out,s,'signac_object.rds'))
  
  return(data.table(brain@meta.data,keep.rownames = 'cell_id'))
},mc.cores = n_cores))

fwrite(mtd_all,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei.csv.gz')
