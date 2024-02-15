
#Transfer label to others samples/nuclei ####
# compute UMAP and store the UMAP model
out<-'outputs/04-ROSMAP_MIT_ATAC/'
dir.create(out)
source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(future)
options(future.globals.maxSize = 10*1000 * 1024^2) #10GB

n_cores_mc=6
n_core_future=4


mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei.csv.gz')
signac_files<-file.path(list.dirs(out,full.names = T,recursive =F ),'signac_object.rds')


#load ref
brain<-readRDS('outputs/04-ROSMAP_MIT_ATAC/brain_12best_samples_qc.rds')
DefaultAssay(brain)<-'peaks'
brain <- RunUMAP(brain, reduction = "lsi", dims = 2:30, return.model = TRUE)



mtd_anno<-rbindlist(mclapply(signac_files,function(f){
  s<-basename(dirname(f))
  message(s)
  
  plan("multicore", workers = n_core_future)
  

  brain.sample<-readRDS(f)
  
  brain.sample <- RunTFIDF(brain.sample)
  brain.sample <- FindTopFeatures(brain.sample, min.cutoff = 'q0')
  brain.sample <- RunSVD(object = brain.sample)
  
  # find transfer anchors
  transfer.anchors <- FindTransferAnchors(
    reference = brain,
    query = brain.sample,
    reference.reduction = "lsi",
    reduction = "lsiproject",
    dims = 2:30
  )
  
  # map query onto the reference dataset
  brain.sample <- MapQuery(
    anchorset = transfer.anchors,
    reference = brain,
    query = brain.sample,
    refdata = brain$cell_type,
    reference.reduction = "lsi",
    new.reduction.name = "ref.lsi",
    reduction.model = 'umap'
  )
  
  brain.sample$cell_type<-brain.sample$predicted.id
  
  saveRDS(brain.sample,fp(out,s,'signac_object.rds'))
  
  return(data.table(brain@meta.data,keep.rownames = 'cell_id'))
},mc.cores = n_cores_mc))

mtd<-merge(mtd,mtd_anno[,.(cell_id,libraryID,cell_type)])

fwrite(mtd,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')

