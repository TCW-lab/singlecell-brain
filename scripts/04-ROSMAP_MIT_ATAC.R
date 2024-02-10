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
RunQsub('scripts/04B-joincall_peak_cell_matrix_per_individual.R',job_name = 'tabix')        


#2) all the different cell types
#first, need annotate and merge the objects

#30 test a first experiment

