out<-'outputs/04-ROSMAP_MIT_ATAC/'
dir.create(out)


library(Seurat)
library(Signac)
library(GenomicRanges)
source('../../utils/r_utils.R')

#Count 
#get 18 state model of DLPFC (E073) STATES FOR EACH 200bp BIN:
#in https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state
#for 1chr
system('wget -O ref-data/ChromHMM_brain_dlpfc/E073_18_core_K27ac_chr1_statebyline.txt.gz https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/STATEBYLINE/E073_18_core_K27ac_chr1_statebyline.txt.gz')

'https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/STATEBYLINE/E073_18_core_K27ac_chr1_statebyline.txt.gz'

chr1<-fread('E073_18_core_K27ac_chr1_statebyline.txt.gz',skip = 1,col.names = 'state')
unique(chr1$state)
#add pos
chr1[,chr:='chr1'][,start:=(0:(.N-1))*200][,end:=(1:.N)*200-1]
#add states anno
states<-fread('ref-data/chromHMM18states.csv')
chr1<-merge(chr1,states)

#for all 
#download them
cmds=sapply(c(1:22,'X','Y'),function(c)paste('wget -O',ps('ref-data/ChromHMM_brain_dlpfc/E073_18_core_K27ac_chr',c,'_statebyline.txt.gz'),
                                             ps('https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/STATEBYLINE/E073_18_core_K27ac_chr',c,'_statebyline.txt.gz')))
sapply(cmds, function(cmd)system(cmd))

#open add pos, anno and bind them
states<-fread('ref-data/chromHMM18states.csv')

chromst<-rbindlist(lapply(paste0('chr',c(1:22,'X','Y')), function(chr){
  chrx<-fread(fp('ref-data/ChromHMM_brain_dlpfc/',ps('E073_18_core_K27ac_',chr,'_statebyline.txt.gz')),skip = 1,col.names = 'state') 
  #add pos
  chrx[,chr:=chr][,start:=(0:(.N-1))*200][,end:=(1:.N)*200-1]
  #add states anno
  chrx<-merge(chrx,states)
  return(chrx)
}))

class(chromst$start)
chromst[chr=='chr1']
chromst[,start:=as.numeric(start)]
fwrite(chromst[order(chr,start)],'ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_200pb_bin_all_chr.csv.gz')

fwrite(chromst[order(chr,start)][,.(chr,start,end,state)],
       'ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_200pb_bin_all_chr.bed.gz',
       sep = '\t',col.names = FALSE,scipen = 999)

#create bed file for each state
out_chromref<-'ref-data/ChromHMM_brain_dlpfc/bed_file_per_state'
dir.create(out_chromref)

for(s in unique(chromst$state)){
  fwrite(chromst[state==s][order(chr,start)][,.(chr,start,end,state)],
         fp(out_chromref,paste0('state',s,'_200pb_bin.bed.gz')),
         sep = '\t',col.names = FALSE,scipen = 999)
  
}
#count by cells for each cell type
chromst<-fread('../APOE4_aging/ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_200pb_bin_all_chr.csv.gz')

#source('../../utils/splicing_deficiency.R')
#for astro####
astro<-readRDS('../singlecell-brain/outputs/04-ROSMAP_MIT_ATAC/Astro.rds')
frag<-Fragments(astro[['peaksCT']])[[1]]
frag@cells
colnames(astro)
frag<-readRDS('../singlecell-brain/outputs/04-ROSMAP_MIT_ATAC/')

unique(astro$libraryID)
colnames(astro[,astro$libraryID=='D19-12513'])
FeatureMatrix()
?FeatureMatrix
Fragments(astro)
# res_overlap<-CountFragmentsBedOverlap(astro,
#                          genomic_regions_file ='ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_200pb_bin_all_chr.bed.gz')
# 
#take a while because of cells extraction from the fragment file, have to optimize it 
#for now try Signac::CountsInRegion
?Signac::CountsInRegion


states_mat<-sapply(unique(chromst$state),function(s){
  message('state ',s)
  state<-fread(ps('../APOE4_aging/ref-data/ChromHMM_brain_dlpfc/bed_file_per_state/state',s,'_200pb_bin.bed.gz'),col.names = c('chr','start','end','state'))
  state_gr<-makeGRangesFromDataFrame(state)
  CountsInRegion(astro,assay = 'peaksCT',regions = state_gr)
  
  
})
colnames(states_mat)<-paste0('state',unique(chromst$state))
astro[['ChromHMM18DLPFC']]<-CreateAssayObject(t(states_mat))
head(astro[[]])
astro<-NormalizeData(astro,assay = 'ChromHMM18DLPFC',normalization.method = 'CLR')
VlnPlot(astro,assay = 'ChromHMM18DLPFC',features = 'state13',group.by = 'disease',log=T,pt.size = 0)
VlnPlot(astro,assay = 'ChromHMM18DLPFC',features = 'state13',group.by = 'apoe_genotype',log=T,pt.size = 0)
saveRDS(astro,'../singlecell-brain/outputs/04-ROSMAP_MIT_ATAC/Astro.rds')

#For all cell type####
celltype_files<-list.files(out,'Ast|Exc|Inh|Mic|Oligo|OPC|VLMC',full.names = T)
for(file in celltype_files){
  message(file)
  celltype<-readRDS(file)
  
  states_mat<-sapply(unique(chromst$state),function(s){
    message('state ',s)
    state<-fread(ps('../APOE4_aging/ref-data/ChromHMM_brain_dlpfc/bed_file_per_state/state',s,'_200pb_bin.bed.gz'),col.names = c('chr','start','end','state'))
    state_gr<-makeGRangesFromDataFrame(state)
    CountsInRegion(celltype,assay = 'peaksCT',regions = state_gr)
    
    
  })
  
  colnames(states_mat)<-paste0('state',unique(chromst$state))
  celltype[['ChromHMM18DLPFC']]<-CreateAssayObject(t(states_mat))
  celltype<-NormalizeData(celltype,assay = 'ChromHMM18DLPFC',normalization.method = 'CLR')
  VlnPlot(celltype,assay = 'ChromHMM18DLPFC',features = 'state13',group.by = 'apoe_genotype',log=T,pt.size = 0)+ggtitle(basename(tools::file_path_sans_ext(file)))
  saveRDS(celltype,file)
  
  
}

