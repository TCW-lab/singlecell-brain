out<-'outputs/04-ROSMAP_MIT_ATAC/'
dir.create(out)


library(Seurat)
library(Signac)
library(GenomicRanges)
library(Matrix)
source('../../utils/r_utils.R')
states<-fread('ref-data/chromHMM18states.csv')


celltype_files<-list.files(out,'Ast|Exc|Inh|Mic|Oligo|OPC|VLMC',full.names = T)
for(file in celltype_files){
  message(file)
  celltype<-readRDS(file)
  for(s in states$state){
    message('state ',s)
    state<-fread(ps('ref-data/ChromHMM_brain_dlpfc/bed_file_per_state/state',s,'.bed.gz'),col.names = c('chr','start','end','state_name'))
    state_gr<-makeGRangesFromDataFrame(state)
    count<-CountsInRegion(celltype,assay = 'peaksCT',regions = state_gr)
    celltype[[paste0('n',states[state==s]$name,'_ChromHMM')]]<-count
  }
    saveRDS(celltype,file)
  
  
}



