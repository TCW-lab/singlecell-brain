out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)


#Generate the well annotated signac objects 

#create first the tabix file for every fragment files
fragment_files<-list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/',pattern = 'fragments\\.tsv\\.gz$',full.names = TRUE)

cmds<-sapply(fragment_files,function(f)paste('tabix -p bed',f))
CreateJobFile(cmds,'scripts/04A-tabix_fragments_files.qsub',modules = 'samtools',nThreads = 16)
RunQsub('scripts/04A-tabix_fragments_files.qsub',job_name = 'tabix')        


#1) produce the peak count matrix from fragments.tsv file for each individual
CreateJobForRfile('scripts/04B-joincall_peak_cell_matrix_per_individual.R',nThreads = 28,maxHours = 72)
RunQsub('scripts/04B-joincall_peak_cell_matrix_per_individual.R',job_name = 'PeakMat')        


#1) annotate and QC nuclei
#i) annotate with metadata/QC stats and 
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(brain) <- annotations


#add the clinical data / assay metabdata
mtd<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_assay_snATACseq_metadata.csv')
mtd<-mtd[,.(libraryID,specimenID,sequencingBatch,assay,platform)]
mtd_bio<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv')[]

mtd<-merge(mtd,mtd_bio[assay=='snATACSeq'][,.(individualID,assay,specimenID,specimenIdSource,tissue)])

# #all donors present in snRNA ?
# inds<-mtd_bio[assay=='snATACSeq']$individualID
# length(unique(inds))#91/92, 1 replicates
# nrow(mtd_bio[assay=='snrnaSeq'][individualID%in%inds]) #91/92


#add clinical
mtd_clin<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Metadata/ROSMAP_clinical.csv')
#mtd_clin<-merge(mtd_bio[assay=='snATACSeq'],mtd_clin,by='individualID')
mtd<-merge(mtd,mtd_clin,all.x=TRUE,by='individualID')
nrow(mtd) #92
brain$libraryID<-brain$donor_id
brain$donor_id<-NULL

mtd_cells<-merge(mtd,
                 data.table(brain@meta.data[,c('nCount_peaks','nFeature_peaks','libraryID')],keep.rownames = 'cell_id')[,.(cell_id,libraryID)],by='libraryID')
brain<-AddMetaData(brain,data.frame(mtd_cells,row.names = 'cell_id'))


# ii) nuclei QC

#add % reads falling in peak
# fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/D19-12513.fragments.tsv.gz')
# fragment_files<-list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/',pattern = 'fragments\\.tsv\\.gz$',full.names = TRUE)
# 
# n_fragments<-rbindlist(lapply(fragment_files, function(file){
#   libID=str_remove(basename(file),'.fragments.tsv.gz')
#   fragments<-fread(file,select = 4,col.names = 'cell_id')
#   return(fragments[,.N,by='cell_id'][,libraryID:=libID][N>500])
# }))
# fwrite(n_fragments,fp(out,'n_fragment_per_barcode_alllibs.csv'))
# 
# table(n_fragments$libraryID)
# 
# fwrite(n_fragments,fp(out,'n_fragment_per_cells_donor.csv'))
CreateJobForRfile('scripts/04C-calculate_N_read_per_barcode_for_cellQC.R',nThreads = 28)
RunQsub('scripts/04C-calculate_N_read_per_barcode_for_cellQC.R',job_name ='nRead')

n_fragments<-fread(fp(out,'n_fragment_per_barcode_alllibs.csv.gz'))
n_fragments[,cell_id:=paste(libraryID,cell_id,sep='_')]
brain<-AddMetaData(brain,n_fragments)
head(brain[[]])


#add Nucleosome signal metrics, TSS enrichment
brain <- NucleosomeSignal(object = brain)
brain<-TSSEnrichment(brain)

#filter nuclei


saveRDS(brain,fp(out,'brain_atac_6_first_samples.rds'))



# 3) annotate the celltypes based on snRNA SEA-AD
#Annotate first high quality samples using RNA annotation transfer label using QC_small SEA-AD DLPFC data to do an ATAC reference,
# then transfer label from this reference for remanining samples

#select the donors/samples : young with good stats (nreads in peak, %reads in peak, TSS enrichment and nucleo signal)

# i ) merge the objects
library(EnsDb.Hsapiens.v86)

brain_list<-lapply(list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',recursive = F), function(dir){
  file<-fp(dir,'signac_object.rds')
  if(file.exists(file))
    readRDS(file)
  else
    NULL
})
libids<-list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',recursive = F,full.names = F)
brain<-merge(brain_list[[1]],brain_list[2:length(brain_list)],add.cell.ids = libids[1:6])


#3) annotate the celltypes based on snRNA SEA-AD


#4)cell type level peak calling
#split object per cell type
#5) pseudobulk peak count creation


#6) test a first experiment

