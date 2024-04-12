out<-'outputs/04-ROSMAP_MIT_ATAC/'
dir.create(out)


library(Seurat)
library(Signac)
library(GenomicRanges)
source('../../utils/r_utils.R')

#get 18 state model of DLPFC (E073) STATES FOR EACH 200bp BIN:
#in https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state
#for 1chr
system('wget -O ref-data/ChromHMM_brain_dlpfc/E073_18_core_K27ac_chr1_statebyline.txt.gz https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/STATEBYLINE/E073_18_core_K27ac_chr1_statebyline.txt.gz')
system('wget -O ref-data/ChromHMM_brain_dlpfc/E073_18_core_K27ac_hg38lift_mnemonics.bed.gz https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/E073_18_core_K27ac_hg38lift_mnemonics.bed.gz')



chromst<-fread('ref-data/ChromHMM_brain_dlpfc/E073_18_core_K27ac_hg38lift_mnemonics.bed.gz',col.names = c('chr','start','end','state_name'))

chromst[,state:=str_extract(state_name,'^[0-9]+')]
chromst[,name:=str_remove(state_name,'^[0-9]+_')]
fwrite(chromst,'ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_all_chr.csv.gz')
fwrite(chromst,'ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_all_chr.bed.gz',
       col.names = F,scipen = 999,sep='\t')
unique(chromst$chr)
fwrite(chromst[chr!='ChrM'],'ref-data/ChromHMM_brain_dlpfc/chromHMM_18states_noChrM.bed.gz',
       col.names = F,scipen = 999,sep='\t')

#create bed file for each state
out_chromref<-'ref-data/ChromHMM_brain_dlpfc/bed_file_per_state'
dir.create(out_chromref)

for(s in unique(chromst$state)){
  fwrite(chromst[state==s][order(chr,start)][,.(chr,start,end,state_name)],
         fp(out_chromref,paste0('state',s,'.bed.gz')),
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
#run 04Ki
CreateJobForRfile('../singlecell-brain/scripts/04Ki-CountChromHMM.R',nThreads = 16)
RunQsub('../singlecell-brain/scripts/04Ki-CountChromHMM.R',job_name = 'countCHStates')




