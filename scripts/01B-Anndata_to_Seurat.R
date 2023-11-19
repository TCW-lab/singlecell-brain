
out<-'outputs/01-SEAAD_data'
dir.create(out,recursive = T)


library(Seurat)
library(SeuratData)
library(SeuratDisk)
source('../../utils/r_utils.R')


#Functions####



#ANALYSIS####
#DLFPC
out1<-fp(out,'DLPFC')
#mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))

afiles<-fp(out1,list.files(out1,pattern = '.h5ad'))

lapply(afiles, function(afile){
  print(afile)
  sceasy::convertFormat(afile, from="anndata", to="seurat",
                        outFile=str_replace(afile,'h5ad$','rds'))
  
})

#MTG
out1<-fp(out,'MTG')
#mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))

afiles<-fp(out1,list.files(out1,pattern = '.h5ad'))

lapply(afiles, function(afile){
  print(afile)
  sceasy::convertFormat(afile, from="anndata", to="seurat",
                        outFile=str_replace(afile,'h5ad$','rds'))
  
})


