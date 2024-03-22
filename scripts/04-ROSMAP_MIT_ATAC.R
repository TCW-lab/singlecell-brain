out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)


#Generate the well annotated signac objects 

#create first the tabix file for every fragment files
fragment_files<-list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/',pattern = 'fragments\\.tsv\\.gz$',full.names = TRUE)

cmds<-sapply(fragment_files,function(f)paste('tabix -p bed',f))
CreateJobFile(cmds,'scripts/04A-tabix_fragments_files.qsub',modules = 'samtools',nThreads = 16)
RunQsub('scripts/04A-tabix_fragments_files.qsub',job_name = 'tabix')        


#1) produce the peak count matrix from fragments.tsv file for each individual
CreateJobForRfile('scripts/04B-joincall_peak_cell_matrix_per_individual.R',nThreads = 28,maxHours = 72)
RunQsub('scripts/04B-joincall_peak_cell_matrix_per_individual.R',job_name = 'PeakMat')        


#2) annotate and QC nuclei####
#for 1####
brain<-readRDS('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/D19-12513/signac_object.rds')

#i) annotate with genome annoation

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19


saveRDS(annotations,'ref-data/EnsDb.Hsapiens.v86_annotations.rds')
annotations<-readRDS('ref-data/EnsDb.Hsapiens.v86_annotations.rds')

seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(brain) <- annotations


#ii) add the clinical data / assay metabdata
mtd<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_assay_snATACseq_metadata.csv')
mtd<-mtd[,.(libraryID,specimenID,sequencingBatch,assay,platform)]
mtd_bio<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv')[]

mtd<-merge(mtd,mtd_bio[assay=='snATACSeq'][,.(individualID,assay,specimenID,specimenIdSource,tissue)])

mtd_clin<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/Metadata/ROSMAP_clinical.csv')
#mtd_clin<-merge(mtd_bio[assay=='snATACSeq'],mtd_clin,by='individualID')
mtd<-merge(mtd,mtd_clin,all.x=TRUE,by='individualID')
nrow(mtd) #92
#save mtd
fwrite(mtd,'ref-data/MIT_ROSMAP_Multiomics/snATAC-seq/samples_metadata.csv')
mtd<-fread('ref-data/MIT_ROSMAP_Multiomics/snATAC-seq/samples_metadata.csv')

#race?
table(mtd$race)
#  1  2 
# 91  1 
#2 is african

# #all donors present in snRNA ?
inds<-mtd$individualID
length(unique(inds))#91/92, 1 replicates
nrow(mtd_bio[assay=='snrnaSeq'][individualID%in%inds]) #91,, yes


brain$libraryID<-brain$donor_id
brain$donor_id<-NULL

mtd_cells<-merge(mtd,
                 data.table(brain@meta.data[,c('nCount_peaks','nFeature_peaks','libraryID')],
                            keep.rownames = 'cell_id')[,.(cell_id,libraryID)],by='libraryID')

brain<-AddMetaData(brain,data.frame(mtd_cells,row.names = 'cell_id'))
head(brain[[]])

# ii) nuclei QC

#add % reads falling in peak
#for 1
fragments<-fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/D19-12513.fragments.tsv.gz')

n_fragments<-fragments[,.N,by='cell_id'][N>500]

# for all
CreateJobForRfile('scripts/04C-calculate_N_read_per_barcode_for_cellQC.R',nThreads = 28)
RunQsub('scripts/04C-calculate_N_read_per_barcode_for_cellQC.R',job_name ='nRead')

n_fragments<-fread(fp(out,'n_fragment_per_barcode_alllibs.csv.gz'))
setnames(n_fragments,'N','n_tot_fragments')

brain<-AddMetaData(brain,data.frame(n_fragments[libraryID==unique(brain$libraryID)],row.names = 'cell_id'))

brain$pct_reads_in_peaks<-brain$nCount_peaks/brain$n_tot_fragments
head(brain[[]])

#add Nucleosome signal metrics, TSS enrichment
brain <- NucleosomeSignal(object = brain)
brain<-TSSEnrichment(brain)
brain[[]]
VlnPlot(brain,c('nCount_peaks','pct_reads_in_peaks','TSS.enrichment','nucleosome_signal'),ncol = 2)


saveRDS(brain,'outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/D19-12513/signac_object.rds')


#for all and get the combined metadata####
#run 04D
CreateJobForRfile('scripts/04D-nuclei_anno_and_qc.R',nThreads = 28,maxHours = 12)
RunQsub('scripts/04D-nuclei_anno_and_qc.R',job_name = 'NucAnnoQC')

#QC check nreads in peak, %reads in peak, TSS enrichment and nucleo signal
mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei.csv.gz')
mtd
mtd[,disease:=cogdx%in%4:6]
mtd[,age_decade:=paste0(str_extract(age_death,'^[1-9]'),'0')]
mtd[,age_at_death_num:=as.numeric(ifelse(age_death=='90+',91,age_death))]

ggplot(mtd)+geom_violin(aes(x=libraryID,y=nCount_peaks,fill=disease))+
  scale_y_log10()+
  scale_x_discrete(guide = guide_axis(angle = 60))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8))

ggplot(mtd)+
  geom_violin(aes(x=libraryID,y=pct_reads_in_peaks,fill=disease))+
  scale_x_discrete(guide = guide_axis(angle = 60))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8))

ggplot(mtd)+
  geom_violin(aes(x=libraryID,y=TSS.enrichment,fill=disease))+
  scale_x_discrete(guide = guide_axis(angle = 60))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8))

ggplot(mtd)+
  geom_violin(aes(x=libraryID,y=nucleosome_signal,fill=disease))+
  scale_x_discrete(guide = guide_axis(angle = 60))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8))


mtd[,median_in_peaks:=median(nCount_peaks),by='libraryID']
mtd[,avg_pct_in_peaks:=mean(pct_reads_in_peaks),by='libraryID']
mtd[,avg_tss_enrichment:=mean(TSS.enrichment),by='libraryID']
mtd[,avg_nucleosome_signal:=mean(nucleosome_signal),by='libraryID']

mts<-unique(mtd,by='libraryID')
ggplot(mts)+
  geom_boxplot(aes(x=disease,y=log10(median_in_peaks)))+theme_bw()

mts[,outlier.count_in_peaks:=log10(median_in_peaks)<2.5&log10(median_in_peaks)%in%boxplot.stats(log10(median_in_peaks))$out,by='disease']

ggplot(mts)+
  geom_boxplot(aes(x=disease,y=avg_pct_in_peaks))+theme_bw()

ggplot(mts)+
  geom_boxplot(aes(x=disease,y=avg_tss_enrichment))+theme_bw()

ggplot(mts)+
  geom_boxplot(aes(x=disease,y=avg_nucleosome_signal))+theme_bw()

outliers<-mts[(outlier.count_in_peaks)]$libraryID

mtd[,outlier.count_in_peaks:=libraryID%in%outliers]


#select the top 20% high qual samples : young with good stats (nreads in peak, %reads in peak, TSS enrichment and nucleo signal)
ggplot(mts)+
  geom_point(aes(x=age_at_death_num,y=avg_pct_in_peaks,col=disease))+theme_bw()

ggplot(mts)+
  geom_point(aes(x=age_at_death_num,y=median_in_peaks,col=disease))+theme_bw()

ggplot(mts)+
  geom_point(aes(x=age_at_death_num,y=avg_tss_enrichment,col=disease))+theme_bw()

ggplot(mts)+
  geom_point(aes(x=age_at_death_num,y=avg_nucleosome_signal,col=disease))+theme_bw()

#age, n count peaks, tss enrichment
mts[,score_qual:=rowSums(data.frame(age_at_death_num<=quantile(age_at_death_num,0.2),
                                    avg_pct_in_peaks>=quantile(avg_pct_in_peaks,0.8),
                                    median_in_peaks>=quantile(median_in_peaks,0.8),
                                    avg_tss_enrichment>=quantile(avg_tss_enrichment,0.8)))]
mts[order(-score_qual)]
high_qual_samples <-mts[order(-score_qual)][1:round(nrow(mts)*0.2)]$libraryID


mtd[libraryID%in%high_qual_samples] #79k/440k

setnames(mts,'score_qual','donor_score_qual')
mtd<-merge(mtd,mts[,.(libraryID,donor_score_qual)])

fwrite(mtd,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei.csv.gz')
fwrite(mts,'outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_sample_level.csv.gz')


# 3) annotate the celltypes based on snRNA SEA-AD ####
#Annotate top 20% using RNA annotation and transfer label using QC_small SEA-AD DLPFC data to do an ATAC reference

# i ) merge the objects
library(EnsDb.Hsapiens.v86)

brain_list<-lapply(file.path('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',high_qual_samples), function(dir){
  file<-fp(dir,'signac_object.rds')
  return(readRDS(file))
})

brain<-merge(brain_list[[1]],brain_list[2:length(brain_list)],add.cell.ids = high_qual_samples)

#ii) conserve only the high quqlity one
VlnPlot(brain,
        c('nCount_peaks','pct_reads_in_peaks','TSS.enrichment','nucleosome_signal'),ncol = 2)

VlnPlot(brain,
        c('nCount_peaks','pct_reads_in_peaks','TSS.enrichment','nucleosome_signal'),ncol = 2,pt.size = 0)
VlnPlot(brain,
        c('nCount_peaks'),pt.size = 0,log = T)
VlnPlot(brain,
        c('nucleosome_signal'),pt.size = 0,log = T)
VlnPlot(brain,
        c('TSS.enrichment'),pt.size = 0,y.max = 10)

brain <- subset(
  x = brain,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 0.25 & 
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
brain #64651/79203

#iii) normalization and dim reduction
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(object = brain)
DepthCor(brain)

brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)

brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindClusters(
  object = brain,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(object = brain, label = TRUE) + NoLegend()


# iv) compute gene activities
gene.activities <- GeneActivity(brain)

# add the gene activity matrix to the Seurat object as a new assay
brain[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)
DefaultAssay(brain) <- 'RNA'
FeaturePlot(
  object = brain,
  features = c('OLIG1','SLC17A6','PVALB','SST',"GAD2","NEUROD6","RORB","MBP",'TMEM119','MPO','TREM2','AQP4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

saveRDS(brain,fp(out,'brain_12best_samples_qc.rds'))

#v) Integrating with scRNA-seq data
#run 04E
CreateJobForRfile('scripts/04E-celltype_labeltransfer_from_rna.R',nThreads = 36)
RunQsub('scripts/04E-celltype_labeltransfer_from_rna.R',job_name = 'TrasferLabelRNA')

#check good annotation
brain<-readRDS('outputs/04-ROSMAP_MIT_ATAC/brain_12best_samples_qc.rds')
DefaultAssay(brain)<-'peaks'
DimPlot(brain,reduction = 'umap',group.by = 'cell_type',label=T,repel = T)

mtd<-data.table(brain@meta.data,keep.rownames = 'cell_id')
mtd[,main_cell_type:=str_extract(cell_type,'Oligo|Exc|Inh|Astro|Mic|Endo|VLMC|OPC')]
ggplot(mtd)+geom_bar(aes(x=main_cell_type,fill=main_cell_type))

#top left Exc neurons are far more the main cluster of Exc, poor quaility cells?
DimPlot(brain,reduction = 'umap',group.by = c('individualID','sequencingBatch','cogdx'))

DimPlot(brain,reduction = 'umap',cells.highlight = WhichCells(brain,expression = individualID=='R9781891'))

#R9181891 sample have a completely different profile of Exc neur chromatin access
#beacause bad quqlity sample?
FeaturePlot(brain,reduction = 'umap',features =  c('pct_reads_in_peaks'),label=T,repel = T) #no

FeaturePlot(brain,reduction = 'umap',features =  c('prediction.score.Exc_L5.IT'),label=T,repel = T) #well predicted

brain$logNCount_peaks<-log(brain$nCount_peaks)
FeaturePlot(brain,reduction = 'umap',features = 'logNCount_peaks') #no
#FeaturePlot(brain,reduction = 'umap',features = 'nFeature_peaks')

#different brain region??
mtd[individualID=='R9781891']
mtd_bio<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv')
mtd_bio[individualID=='R9781891']#have scRNA from multiple brain region and cnMultiome from  medial frontal cortex instread of PFC
table(mtd_bio$individualID)#not specific to hom
mtd_bio[,tissue.per.donor:=length(unique(tissue)),by='individualID']
unique(mtd_bio[tissue.per.donor==7]$individualID) #take part of the 4 donors with all brain region
#guess : this sample coming from a different brainregion. annotated as PFC but inversion when taking/identfying the frozen sample
#Note . number of snMultiome donors: 
mtd_bio[assay=='snMultiome']#19 multiome
#maybe to use to brdge, but come from a different brain region ( medial frontal cortex)
#so we can transfer data like that

#4) Transfer label to others samples/nuclei ####
#run 04F
CreateJobForRfile('scripts/04F-MapSamples.R',nThreads = 36)
RunQsub('scripts/04F-MapSamples.R',job_name = 'MapSamples',wait_for ='TrasferLabelRNA' )




#5) flag sample outliers, ####
# annotate main celltype 
mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')
mtd

mtd[,main_cell_type:=str_extract(cell_type,'Oligo|Exc|Inh|Astro|Mic|Endo|VLMC|OPC')]
mtd[,main_cell_type:=factor(main_cell_type,levels = c('Exc','Inh','Oligo','Astro','OPC','Mic','VLMC','Endo'))]
table(mtd$main_cell_type)

# Exc    Inh  Oligo  Astro    OPC    Mic   VLMC   Endo 
# 168288  53321 161413  23174  13445  13762   6908      0 
ggplot(mtd)+
  geom_bar(aes(x=main_cell_type,fill=main_cell_type),position='dodge')+theme_bw()


#i)Clinical Status
mts<-unique(mtd,by=c('libraryID'))
#Age, PMI, ethnicity
#flag missing info
table(mts$cogdx)
ggplot(mts)+geom_bar(aes(x=cogdx,fill=as.factor(cogdx)))+theme_bw( )
mts[,cognitive_status:=ifelse(cogdx==1,'NL',ifelse(cogdx%in%c(2,3),'MCI','AD'))]
ggplot(mts)+geom_bar(aes(x=cognitive_status,fill=cognitive_status))+theme_bw( )

#APOE
table(mts$apoe_genotype) #2 44
# 23 24 33 34 44 
#  9  3 62 16  2 
mts[,apoe4:=apoe_genotype%in%c(34,44,24)]
mts[,apoe2:=apoe_genotype%in%c(24,23,22)]

mts[,outlier.apoe:=apoe4&apoe2]

ggplot(mts)+
  geom_bar(aes(x=as.factor(apoe_genotype),fill=cognitive_status))+theme_bw() 

ggplot(mts)+geom_boxplot(aes(x=as.factor(apoe_genotype),y=age_at_death_num,fill=as.factor(apoe_genotype)))+theme_bw()


#age
mts[,age_decade:=paste0(str_extract(age_death,'^[1-9]'),'0')]
mts[,age_at_death_num:=as.numeric(ifelse(age_death=='90+',91,age_death))]

table(mts$age_decade,mts[]$cognitive_status)
  #    AD MCI NL
  # 70  3   3  5
  # 80 24   8 20
  # 90  7   9 13

ggplot(mts)+geom_boxplot(aes(x=cognitive_status,y=age_at_death_num,fill=cognitive_status))+theme_bw()
ggplot(mts)+geom_boxplot(aes(y=age_at_death_num)) #ok


#PMI
ggplot(mts)+geom_boxplot(aes(x=cognitive_status,y=pmi,fill=cognitive_status))+theme_bw() #MCI have bigger PMI, to take into account in analysis

#ethnicity
table(mts$race)
# 1  2 
# 91  1 

mts[,ethnicity:=ifelse(race==1,'White',ifelse(race==2,'Black or African American',NA))]

table(mts$ethnicity)
# Black or African American                     White 
# 1                        91 
ggplot(mts)+
  geom_bar(aes(x=ethnicity,fill=cognitive_status))+theme_bw() 

mts[,outlier.ethnicity:=ethnicity=='Black or African American']

mts[,outlier.clinical.status:=outlier.ethnicity|outlier.apoe]
mts[(outlier.clinical.status)] #4

#add this information to the main metadata
mtd<-merge(mtd,
           mts[,.SD,.SDcols=c('libraryID',setdiff(colnames(mts),colnames(mtd)))],by='libraryID')



#ii) Cellular/Molecular Traits QC

#based on %read in peak
ggplot(mtd)+geom_violin(aes(x=libraryID,y=pct_reads_in_peaks))
mtd[,med.n_tot_fragment:=median(n_tot_fragments),by=.(libraryID)]
mtd[,med.n_tot_fragment.ct:=median(n_tot_fragments),by=.(libraryID,main_cell_type)]

mtd[,avg.pct.read.in.peak:=mean(`pct_reads_in_peaks`),by=c('libraryID')]
mtd[,avg.pct.read.in.peak.ct:=mean(`pct_reads_in_peaks`),by=c('libraryID','main_cell_type')]

mtd[,med.tss.enrich:=median(`pct_reads_in_peaks`),by=c('libraryID')]
mtd[,med.tss.enrich.ct:=median(`pct_reads_in_peaks`),by=c('libraryID','main_cell_type')]

mtd[,med.nucleosome_signal:=median(`nucleosome_signal`),by=c('libraryID')]
mtd[,med.nucleosome_signal.ct:=median(`nucleosome_signal`),by=c('libraryID','main_cell_type')]

mtsc<-unique(mtd,by=c('libraryID','main_cell_type'))
mts<-unique(mtd,by=c('libraryID'))
ggplot(mts)+
  geom_boxplot(aes(x=cognitive_status,y=med.n_tot_fragment,fill=cognitive_status))+theme_bw()

ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=med.n_tot_fragment.ct,fill=cognitive_status))+theme_bw()

ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=med.n_tot_fragment.ct,fill=as.factor(apoe_genotype)))+theme_bw()


ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=avg.pct.read.in.peak.ct,fill=cognitive_status))+theme_bw()

ggplot(mtsc)+
  geom_boxplot(aes(x=as.factor(apoe_genotype),
                   y=avg.pct.read.in.peak.ct,fill=as.factor(apoe_genotype)))+theme_bw()+
  facet_wrap('main_cell_type')


ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=med.tss.enrich.ct,fill=cognitive_status))+theme_bw()

ggplot(mtsc)+
  geom_boxplot(aes(x=as.factor(apoe_genotype),
                   y=med.tss.enrich.ct,fill=as.factor(apoe_genotype)))+theme_bw()+
  facet_wrap('main_cell_type')



ggplot(mtsc)+
  geom_boxplot(aes(x=as.factor(apoe_genotype),
                   y=med.nucleosome_signal.ct,fill=as.factor(apoe_genotype)))+theme_bw()+
  facet_wrap('main_cell_type')

ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=med.nucleosome_signal.ct,fill=cognitive_status))+theme_bw()

#no removal because associated to Dementia


#based on Cellular Distribution
#Cell Distribution
ggplot(mtd)+
  geom_bar(aes(x=main_cell_type,fill=main_cell_type),position='dodge')+
  facet_wrap(~libraryID)+theme_bw()

mtd[,n.cells:=.N,by=.(`libraryID`)]
mtd[,pct.ct:=.N/n.cells,by=.(`libraryID`,main_cell_type)]
mtd[,n.ct:=.N,by=.(`libraryID`,main_cell_type)]
mtd[,n.ct.all:=.N,by=.(main_cell_type)]


mtsc<-unique(mtd,by=c('libraryID','main_cell_type'))

ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=pct.ct,fill=`cognitive_status`),position='dodge')+
  theme_bw()

ggplot(mtsc)+
  geom_boxplot(aes(x=main_cell_type,y=pct.ct,fill=apoe4),position='dodge')+
  theme_bw()


#Summurize Cell type proportions abnormalities at donor level
#=> PCA of the % deviation from IQR
source('../../utils/pca_utils.R')
mtsc[,pct.ct.q25:=quantile(pct.ct,0.25),by='main_cell_type']
mtsc[,pct.ct.q75:=quantile(pct.ct,0.75),by='main_cell_type']

mtsc[,pct.ct.IQR:=pct.ct.q75-pct.ct.q25]
mtsc[pct.ct<pct.ct.q25,pct.from.IQR:=(pct.ct.q25-pct.ct)/pct.ct.IQR]
mtsc[pct.ct>pct.ct.q75,pct.from.IQR:=(pct.ct-pct.ct.q75)/pct.ct.IQR]
mtsc[is.na(pct.from.IQR),pct.from.IQR:=0]
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.003323  0.327975  0.362668 13.114971 

iqr_skew_mat<-dcast(mtsc,libraryID~main_cell_type,value.var = 'pct.from.IQR')
iqr_skew_pca<-RunPca(t(data.frame(iqr_skew_mat,row.names = 'libraryID')),scale = T)
iqr_skew_pca$x

iqr_skew_pca_dt<-PcaPlot(iqr_skew_pca,mts,group.by ='cognitive_status',sample_col = 'libraryID',return_pcs_mtd = T)
ggplot(iqr_skew_pca_dt)+geom_boxplot(aes(x=cognitive_status,y=PC1))

#outlier if sample > 3*IQR of PC1
iqr_skew_pca_dt[,outlier.cellprop:=PC1%in%boxplot.stats(PC1,coef=3)$out]
iqr_skew_pca_dt[(outlier.cellprop)]

PcaPlot(iqr_skew_pca,mts,group.by ='cognitive_status',
        sample_col = 'libraryID',label = iqr_skew_pca_dt$outlier.cellprop)

ggplot(mtd[libraryID%in%iqr_skew_pca_dt[(outlier.cellprop)]$libraryID])+
  geom_bar(aes(x=main_cell_type,fill=main_cell_type),position='dodge')+
  facet_wrap(~libraryID)+theme_bw() #samples with very high Oligo /Very low Exc Neu

iqr_skew_pca_dt[,cell_prop_dev_pc1:=PC1]


mtsc<-merge(mtsc,iqr_skew_pca_dt[,.(libraryID,outlier.cellprop,cell_prop_dev_pc1,outlier.cellular.status)])
mtd<-merge(mtd,iqr_skew_pca_dt[,.(libraryID,outlier.cellprop,cell_prop_dev_pc1,outlier.cellular.status)])
unique(mtd[(outlier.cellprop)],by='libraryID')
      
#based on chromatin profile
#UMAP
umap_coords<-rbindlist(lapply(list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',
                                        full.names = T,recursive = F), function(d){
                                          s<-basename(d)
                                          f<-fp(d,'signac_object.rds')
                                          sign<-readRDS(f)
                                          data.table(sign@reductions$ref.umap@cell.embeddings,keep.rownames = 'cell_id')[,libraryID:=s]
                                          
                                        }))
mtd<-merge(mtd,umap_coords,by=c('cell_id','libraryID'))

ggplot(mtd)+geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=cell_type),size=0.2)+theme_bw()
ggplot(mtd)+geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=main_cell_type),size=0.2)+theme_bw()
ggplot(mtd)+geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=libraryID),size=0.2)+theme_bw()+NoLegend()
table(mtd[refUMAP_1<(-5)&refUMAP_2>(7.5)]$libraryID)
exc.donor.outlier<-names(which(table(mtd[refUMAP_1<(-5)&refUMAP_2>(7.5)]$libraryID)>50))
ggplot(mtd[libraryID%in%exc.donor.outlier&main_cell_type=='Exc'])+
  geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=cell_type),size=0.2)+
  theme_bw()+facet_wrap('libraryID')

#flag as outlier D19-12532 and flag all Exc cells cluster as outlier

ggplot(mtd)+
  geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=cell_type),size=0.2)+
  theme_bw()+geom_hline(yintercept = 7)+geom_vline(xintercept = -5)

table(mtd[refUMAP_1<(-5)&refUMAP_2>7]$cell_type)

mtd[,exc.donor.outlier:=libraryID=='D19-12532']
mtd[,exc.cells.outlier:=refUMAP_1<(-5)&refUMAP_2>7]

#oligo outlier
table(mtd[refUMAP_1>4&refUMAP_2>7.5]$cell_type)
table(mtd[refUMAP_1>4&refUMAP_2>7.5]$libraryID) #mainly D19-125455

ggplot(mtd[libraryID=='D19-125455'&main_cell_type=='Oligo'])+
  geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=cell_type),size=0.2)+
  theme_bw()+facet_wrap('libraryID') #have some 'normal' oligo

ggplot(mtd)+
  geom_point(aes(x=refUMAP_1,y=refUMAP_2,col=cell_type),size=0.2)+
  theme_bw()+facet_wrap('libraryID')+facet_wrap(~libraryID=='D19-125455') #seems pretty normal for other cell type

mtd[libraryID=='D19-125455']
table(mts$spanish)
#do not seems specifc known phenotype. flag as oligo.outlier
mtd[,oligo.donor.outlier:=libraryID=='D19-125455']
mtd[,oligo.cells.outlier:=refUMAP_1>4&refUMAP_2>7.5]


#final donors outliers
mtd[,outlier.cellular.status:=outlier.cellprop|oligo.donor.outlier|exc.donor.outlier]

mtd[,donor.outlier:=outlier.cellular.status|outlier.clinical.status]

mts<-unique(mtd,by=c('libraryID'))
tot.outliers<-mts[(donor.outlier)]$libraryID
length(tot.outliers) #12/92

#save 
fwrite(mtd,fp(out,'all_final_ATACseq_nuclei_metadata.csv.gz'))

mtsc<-UniqueClean(mtd,key_cols=c('libraryID','main_cell_type'),pattern_to_exclude = 'cell_id')

fwrite(mtsc,fp(out,'all_final_ATACseq_nuclei_main_cell_type_level_metadata.csv.gz'))

mts<-UniqueClean(mtsc,key_cols=c('libraryID'),pattern_to_exclude = 'cell_type')

fwrite(mts,fp(out,'all_final_ATACseq_nuclei_sample_level_metadata.csv.gz'))



#6)cell type level peak calling####
#peak call in the good quality top20%nuclei, and produce the faeture matrix for every samples
#to then QC nuclei according to these peaks
#peak call
CreateJobForRfile('scripts/04G-PeakCallCellType.R',nThreads = 28)
RunQsub('scripts/04G-PeakCallCellType.R',job_name = 'CellTypePeak')

#feature matrix for all samples
CreateJobForRfile('scripts/04Gi-CountPerSample.R',nThreads = 36)
RunQsub('scripts/04Gi-CountPerSample.R',job_name = 'CTPeakCount',wait_for ='CellTypePeak' )


#check how many new peak and cell specificity
peaks<-readRDS(fp(out,'brain_12best_samples_qc_celltype_peaks.rds'))
mtd<-fread(fp(out,'metadata_all_nuclei_celltype_annotated.csv.gz'))

peaks_dt<-data.table(as.data.frame(peaks))
peaks_dt[,peak_id:=paste(seqnames,start,end,sep="-")]
peaks_dt[,chr:=seqnames]
peaks_dt<-Reduce(rbind,lapply(str_replace_all(unique(mtd$cell_type)," |/","_"), function(x){
  ocrs<-peaks_dt[str_detect(peak_called_in,x)]
  return(ocrs[,cell_type:=x][,-'peak_called_in'])
}))
peaks_dt[,n.ct:=.N,by="peak_id"]
length(unique(peaks_dt$peak_id))#318278

fwrite(peaks_dt,fp(out,"brain_12best_samples_qc_celltype_peaks.csv.gz"))

peaks_samples<-fread('outputs/04-ROSMAP_MIT_ATAC/initial_peak_calling_by_indivs_tidy.csv.gz')
length(unique(peaks_samples$peak_id))#401697
#more peak in initial call, but more cell type specific peak?
#=> %frag_in_peak before / after
mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')
mtd$pct_frag_in_peaksCT<- mtd$nCount_peaksCT/mtd$n_tot_fragments
mtd2<-rbind(mtd[,.(cell_id,cell_type,pct_reads_in_peaks)][,peak_calling:='per donor'],mtd[,.(cell_id,cell_type,pct_frag_in_peaksCT)][,pct_reads_in_peaks:=pct_frag_in_peaksCT][,peak_calling:='per cell type'],fill=T)
ggplot(mtd2)+geom_boxplot(aes(x=cell_type,y=pct_reads_in_peaks,fill=peak_calling))+theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))

#7) merge sample object per main cell type, ####
#for one
unique(mtd[,.(main_cell_type)])
astro_list<-lapply(list.dirs('outputs/04-ROSMAP_MIT_ATAC/peak_count_matrices/',
                             full.names = T,recursive = F)[1:2], function(d){
                               s<-basename(d)
                               f<-fp(d,'signac_object.rds')
                               brain<-readRDS(f)
                               brain<-AddMetaData(brain,data.frame(mtd[libraryID==s],row.names = 'cell_id'))
                               brain<-RenameCells(brain,add.cell.id =s)
                               brain<-SplitObject(brain,split.by='main_cell_type')[['Astro']]
                               return(brain)
                             })
astro_list
astro<-merge(astro_list[[1]],astro_list[2:length(astro_list)], merge.dr = c('ref.lsi','ref.umap'),merge.data=FALSE)
saveRDS(astro,fp(out,'Astro.rds'))

#for all
CreateJobForRfile('scripts/04H-MergeMainCellType.R',nThreads = 28)
RunQsub('scripts/04H-MergeMainCellType.R',job_name = 'mergeMain')


#8) cell type level QC : 
#for 1
#3* IQR? of nCount, nFeature, pct read, tss enrichment etc
Exc<-readRDS('outputs/04-ROSMAP_MIT_ATAC/Exc.rds')
head(Exc[[]])
VlnPlot(Exc,features =  c('pct_frag_in_peaksCT','nFeature_peaksCT','nCount_peaksCT'),pt.size = 0)
VlnPlot(Exc,features =  c('pct_frag_in_peaksCT','nFeature_peaksCT','nCount_peaksCT'),pt.size = 0,split.by = 'disease')
VlnPlot(Exc,features =  c('nFeature_peaksCT','nCount_peaksCT'),pt.size = 0,split.by = 'disease',log=T)

VlnPlot(Exc,features =  c('TSS.enrichment','nucleosome_signal'),split.by = 'disease',pt.size = 0,log=T)
#ggsave(fp(out1,ps(ct,'qc_cell_metrics.png')),width = 8,height = 6)

boxres<-boxplot.stats(Exc$nCount_peaksCT,coef = 3)
Exc$nCount_peaksCT.high<-colnames(Exc)%in%names(boxres$out)
sum(Exc$nCount_peaksCT.high) #2908/168288 #2%
table(Exc$nCount_peaksCT.outlier,Exc$libraryID)#but a removed from one donor
VlnPlot(Exc,features =  'nCount_peaksCT',split.by = 'nCount_peaksCT.outlier',pt.size = 0,log=T)
VlnPlot(subset(Exc,libraryID=='D19-12521'),features =  'nCount_peaksCT')

#ggsave(fp(out2,ps(ct,'qc_nUMIs.png')),width = 4,height = 6)

boxres<-boxplot.stats(Exc$nFeature_peaksCT,coef = 3)
Exc$nFeature_peaksCT.high<-colnames(Exc)%in%names(boxres$out)
sum(Exc$nFeature_peaksCT.high) #2133/168288 #1%
table(Exc$nFeature_peaksCT.high,Exc$libraryID)
VlnPlot(Exc,features =  'nFeature_peaksCT',split.by = 'nFeature_peaksCT.high',pt.size = 0,log=T)
#can exclude based on that
Exc$nFeature_peaksCT.outlier<-Exc$nFeature_peaksCT.high


VlnPlot(Exc,features =  'pct_frag_in_peaksCT',pt.size = 0,log=T)
boxres<-boxplot.stats(Exc$pct_frag_in_peaksCT,coef = 3)
Exc$pct_frag_in_peaksCT.high<-colnames(Exc)%in%names(boxres$out)
sum(Exc$pct_frag_in_peaksCT.high) #0
table(Exc$pct_frag_in_peaksCT.high,Exc$libraryID)
VlnPlot(Exc,features =  'pct_frag_in_peaksCT',split.by = 'nFeature_peaksCT.high',pt.size = 0,log=T)
#can't exclude based on that
#Exc$nFeature_peaksCT.outlier<-Exc$nFeature_peaksCT.high


VlnPlot(Exc,features =  'TSS.enrichment',pt.size = 0,log=T)
boxres<-boxplot.stats(Exc$TSS.enrichment,coef = 3)
Exc$TSS.enrichment.high<-colnames(Exc)%in%names(boxres$out)
sum(Exc$TSS.enrichment.high) #496
table(Exc$TSS.enrichment.high,Exc$libraryID)
VlnPlot(Exc,features =  'TSS.enrichment',split.by = 'TSS.enrichment.high',pt.size = 0,log=T)
#can't exclude based on that because only loww tss enrichemnnt is bad
#do normal outlier detection to detect low score outlier
boxres<-boxplot.stats(Exc$TSS.enrichment)
Exc$TSS.enrichment.low<-colnames(Exc)%in%names(boxres$out)&Exc$TSS.enrichment<boxres$stats[2]
sum(Exc$TSS.enrichment.low) #29
Exc$TSS.enrichment.outlier<-Exc$TSS.enrichment.low


VlnPlot(Exc,features =  'nucleosome_signal',pt.size = 0,log=T)
boxres<-boxplot.stats(Exc$nucleosome_signal,coef = 3)
Exc$nucleosome_signal.high<-colnames(Exc)%in%names(boxres$out)
sum(Exc$nucleosome_signal.high) #2876
table(Exc$nucleosome_signal.high,Exc$libraryID) #outliers fairly distribute accross donors
VlnPlot(Exc,features =  'nucleosome_signal',split.by = 'nucleosome_signal.high',pt.size = 0,log=T)
#can exclude based on that
Exc$nucleosome_signal.outlier<-Exc$nucleosome_signal.high

#Final outlier is:
Exc$cell.outlier<-Exc$nucleosome_signal.outlier|Exc$TSS.enrichment.outlier|Exc$TSS.enrichment.outlier|Exc$nFeature_peaksCT.outlier


Exc$outlier<-Exc$donor.outlier|Exc$cell.outlier

message(round(sum(Exc$donor.outlier)/nrow(Exc),digits = 2),'% nuclei flagged as donors outliers')
message(round(sum(Exc$cell.outlier)/nrow(Exc),digits = 2),'% nuclei flagged as cells outliers' )

message('In total ',round(sum(Exc$outlier)/nrow(Exc),digits = 2),'% nuclei flagged as outliers')

#save the full object
#saveRDS(exc,file)

#for all
#run 04I
CreateJobForRfile('scripts/04I-cell_level_QC.R',nThreads = 28)
RunQsub('scripts/04I-cell_level_QC.R',wait_for = 'mergeMain',job_name = 'CellQC')


#9)DonorsGroups-cell type level peak calling####
# to increase discovery of subpop specific peak, call per groups of donors (non outliers donors only).
# for each celltype group donors according to KNN/SNN celltype spe PCA of celltype spe peak count matrix
# (group can be different depending of cell type)
#for 1 ct
astro<-readRDS('outputs/04-ROSMAP_MIT_ATAC/Astro.rds')

#rm outlier
astro<-astro[,!astro$outlier]
table(astro$libraryID)

#per cell type 
#1) peak call per cluster

#at least 20 cells per donor
samples_to_keep<-names(which(table(astro$libraryID)>20))
astro<-astro[,astro$libraryID%in%samples_to_keep]

astro_pseudo<-AggregateExpression(astro,slot = 'data',
                                  assays = 'peaksCT',return.seurat = T,group.by = 'libraryID')

mtd<-data.table(astro@meta.data)
mtd[,median_n_tot_fragment:=median(n_tot_fragments),by='libraryID']
mts<-unique(mtd,by='libraryID')

astro_pseudo<-AddMetaData(astro_pseudo,data.frame(mts,row.names = 'libraryID'))

astro_pseudo <- FindTopFeatures(astro_pseudo, min.cutoff = 'q0',assays = 'peaksCT')
astro_pseudo <- RunSVD(object = astro_pseudo)
DepthCor(astro_pseudo)

astro_pseudo <- RunUMAP(
  object = astro_pseudo,
  reduction = 'lsi',
  dims = 1:6
)

astro_pseudo <- FindNeighbors(
  object = astro_pseudo,
  reduction = 'lsi',
  dims = 1:6
)
astro_pseudo <- FindClusters(
  object = astro_pseudo,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(astro_pseudo)
DimPlot(astro_pseudo,group.by = 'cognitive_status')
DimPlot(astro_pseudo,group.by = 'apoe4')
DimPlot(astro_pseudo,group.by = 'apoe2')
FeaturePlot(astro_pseudo, 'avg.pct.read.in.peak.ct')
DimPlot(astro_pseudo,group.by = c('msex'))


FeaturePlot(astro_pseudo, c('n.ct','avg.pct.read.in.peak.ct'))
FeaturePlot(astro_pseudo, c('braaksc','pmi'))

#save donor cluster
mtd[cell_type==ct,donor.ct.group:=astro_pseudo@meta.data[libraryID,'seurat_clusters']]
astro<-AddMetaData(astro,data.frame(mtd[cell_type==ct][,.(cell_id_long,donor.ct.group)],
                                            row.names = 'cell_id_long'))

#peak call per cluster
peaks<-CallPeaks(astro,group.by=c('donor.ct.group'),
                 macs2.path = '/projectnb/tcwlab/LabMember/adpelle1/micromamba/envs/macs2/bin/macs2')
peaks$peak_called_in<-paste0(ct,peaks$peak_called_in)


#for all
#run 04J

CreateJobForRfile('scripts/04J-PeakCallDonorGroup.R',nThreads = 28,maxHours = 48)
RunQsub('scripts/04J-PeakCallDonorGroup.R',job_name = 'PeakDonorGroup',wait_for = 'mergeMain')

peaks<-readRDS(fp(out,'perDonorGroups_celltype_peaks.rds'))
mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')

peaks_dt<-data.table(as.data.frame(peaks))
peaks_dt[,peak_id:=paste(seqnames,start,end,sep="-")]
peaks_dt[,chr:=seqnames]
peaks_dt<-Reduce(rbind,lapply(unique(mtd$cell_type), function(x){
  ocrs<-peaks_dt[str_detect(peak_called_in,x)]
  return(ocrs[,cell_type:=x][,-'peak_called_in'])
}))
peaks_dt[,n.ct:=.N,by="peak_id"]
length(unique(peaks_dt$peak_id))#491503 - 318278 = 173k new peak

fwrite(peaks_dt,fp(out,"perDonorGroups_celltype_peaks.csv.gz"))
peaks_dt<-fread(fp(out,"perDonorGroups_celltype_peaks.csv.gz"))

#celltype spe peak
peaksct_dt<-fread(fp(out,"brain_12best_samples_qc_celltype_peaks.csv.gz"))
table(peaksct_dt$cell_type)
peaks_dt<-rbind(peaks_dt[,call:='PerSampleGroups'],peaksct_dt[,call:='12BestSamples'])
peaks_dt[,cell_type:=str_replace_all(cell_type,' ','_')]
peaks_dt[,cell_type:=str_replace_all(cell_type,'/','_')]

ggplot(peaks_dt)+geom_bar(aes(x=cell_type,fill=call),position = 'dodge')+
  scale_x_discrete(guide = guide_axis(angle = 60))+theme_bw()
peaks_dt[,celltype_spe:=.N==1,by=.(peak_id,call)]
ggplot(peaks_dt[(celltype_spe)])+geom_bar(aes(x=cell_type,fill=call),position = 'dodge')+
  scale_x_discrete(guide = guide_axis(angle = 60))+theme_bw()

#in each celltype, n of sample group spe peaks
peaks_dt<-data.table(as.data.frame(peaks))
peaks_dt[,peak_id:=paste(seqnames,start,end,sep="-")]
peaks_dt[,chr:=seqnames]
peaks_l<-lapply(unique(mtd$cell_type), function(ct){
  ocrs<-peaks_dt[str_detect(peak_called_in,ct)]
  return(ocrs[,cell_type:=ct][,peak_called_in:=str_extract(peak_called_in,paste0(ct,'[0-9,]+'))])
})

peaks_l[[1]]
unique(unlist(strsplit(str_remove(peaks_l[[1]]$peak_called_in,unique(peaks_l[[1]]$cell_type)),',')))
peaks_dt<-rbindlist(lapply(peaks_l, function(x){
  
  groups<-unique(unlist(strsplit(str_remove(x$peak_called_in,unique(x$cell_type)),',')))
  ocrs<-rbindlist(lapply(groups, function(g){
    ocr<-x[str_detect(peak_called_in,g)]
    return(ocr[,group:=g][,-'peak_called_in'])
  }))
  return(ocrs)
}))
peaks_dt
ggplot(peaks_dt)+geom_bar(aes(x=cell_type,fill=group),position = 'dodge')+
  scale_x_discrete(guide = guide_axis(angle = 60))+theme_bw()

peaks_dt[,group_spe:=.N==1,by=.(peak_id,cell_type)]

ggplot(peaks_dt[(group_spe)])+geom_bar(aes(x=cell_type,fill=group),position = 'dodge')+
  scale_x_discrete(guide = guide_axis(angle = 60))+theme_bw()


fwrite(peaks_dt,fp(out,"perDonorGroups_celltype_peaks_groupsinfos.csv.gz"))

#create new feature matrix 
CreateJobForRfile('scripts/04Ji-CountPerCelltype.R',nThreads = 28,maxHours = 48)
RunQsub('scripts/04Ji-CountPerCelltype.R',
        job_name = 'DonorGroupPeakCount',wait_for ='PeakDonorGroup' )

#=> %frag_in_peak before / after
mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/metadata_all_nuclei_celltype_annotated.csv.gz')
mtd[,clinical_status:=ifelse(cognitive_status=='AD',"AD",'NL')]
mtd[,clinical_status:=factor(clinical_status,levels = c('NL','AD'))]
mtd[,APOE4_carrier:=ifelse(apoe4,"yes",'no')]

#mtd$pct_frag_in_peaksCT<- mtd$nCount_peaksCT/mtd$n_tot_fragments
mtd2<-rbind(mtd[,.(cell_id,cell_type,pct_frag_in_peaksCT)][,pct_frag_in_peaks:=pct_frag_in_peaksCT][,peak_calling:='12BestSamples'],mtd[,.(cell_id,cell_type,pct_frag_in_peaksDCT)][,pct_frag_in_peaks:=pct_frag_in_peaksDCT][,peak_calling:='PerSampleGroups'],fill=T)

mtd2<-rbind(mtd2,mtd[,.(cell_id,cell_type,pct_reads_in_peaks)][,pct_frag_in_peaks:=pct_reads_in_peaks][,peak_calling:='bulk_per_donor'],fill=T)
mtd2[,peak_calling:=factor(peak_calling,levels = c('bulk_per_donor','12BestSamples','PerSampleGroups'))]

ggplot(mtd2)+geom_boxplot(aes(x=cell_type,y=pct_frag_in_peaks,fill=peak_calling))+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))

#still less AD/APOE4 in peak? 
ggplot(mtd)+geom_boxplot(aes(x=cell_type,y=pct_frag_in_peaksDCT,fill=disease))+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))

ggplot(mtd)+geom_boxplot(aes(x=cell_type,y=pct_frag_in_peaksDCT,fill=apoe4))+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))

#APOE4 Astrocyte in NL for grant
mtdf<-mtd[clinical_status=='NL'&cell_type=='Astrocyte']


#APOE carrier yes no
ggplot(mtdf)+
  geom_violin(aes(x=individualID,
                  fill=APOE4_carrier,y=pct_frag_in_peaksDCT))+
  scale_y_log10()+scale_x_discrete(limits=unique(mtdf[order(APOE4_carrier)]$individualID))

ggplot(mtdf)+geom_boxplot(aes(x=individualID,fill=APOE4_carrier,y=pct_frag_in_peaksDCT))


#sample lvel - main cell type
mtd[,avg.pct.read.in.peak:=mean(pct_frag_in_peaksDCT),by=c('libraryID')]
mtd[,avg.pct.read.in.peak.ct:=mean(pct_frag_in_peaksDCT),by=c('libraryID','main_cell_type')]

mtsc<-unique(mtd[!(outlier)],by=c('libraryID','main_cell_type'))

mtsc[,transposition.event.outside.peak:=1-avg.pct.read.in.peak.ct]

#AD vs not
ggplot(mtsc)+geom_boxplot(aes(x=main_cell_type,y=transposition.event.outside.peak,fill=clinical_status))+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))



mtsc[,wilcox.test(transposition.event.outside.peak[(disease)],
                  transposition.event.outside.peak[!(disease)])$p.value,by=.(main_cell_type)]

# main_cell_type         V1
# 1:          Oligo 0.22989356
# 2:            Exc 0.02724819
# 3:            Inh 0.08010532
# 4:          Astro 0.28079516
# 5:            Mic 0.65746591
# 6:            OPC 0.06260327
# 7:           VLMC 0.84244971

#APOE4 


ggplot(mtsc,aes(x=APOE4_carrier,y=transposition.event.outside.peak,col=apoe4))+
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(width=0.4,size=0.4)+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))+
  facet_grid(clinical_status~main_cell_type)

mtsc[,wilcox.test(transposition.event.outside.peak[(apoe4)],
                  transposition.event.outside.peak[!(apoe4)])$p.value,by=.(disease,main_cell_type)]
# disease main_cell_type         V1
# 1:    TRUE          Oligo 0.89662106
# 2:    TRUE            Exc 0.48146031
# 3:    TRUE            Inh 0.25960204
# 4:    TRUE          Astro 0.89662106
# 5:    TRUE            Mic 1.00000000
# 6:    TRUE            OPC 0.69643413
# 7:    TRUE           VLMC 0.74172333
# 8:   FALSE            Inh 0.06502848
# 9:   FALSE          Oligo 0.31399012
# 10:   FALSE            Mic 0.17792188
# 11:   FALSE            Exc 0.13587023
# 12:   FALSE            OPC 0.10157674
# 13:   FALSE          Astro 0.01479927 *
# 14:   FALSE           VLMC 0.26782091

#APOE4 Astrocyte in NL for grant
mtscf<-mtsc[clinical_status=='NL'&cell_type=='Astrocyte']

mtscf[,apoe4:=factor(apoe4,levels = c(TRUE,FALSE))]

#n apoe4 
table(mtscf$apoe4)
 # TRUE FALSE 
 #    6    47 
table(mtscf$apoe_genotype)


ggplot(mtscf,aes(x=APOE4_carrier,y=transposition.event.outside.peak,col=apoe4))+
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(width=0.4,size=1)+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))
ggsave('outputs/04-ROSMAP_MIT_ATAC/astrocyte_NL_heterochromatin_figure_for_aging_grant.pdf',width = 5,height = 5)

#barplot +/- std error
# Calculate mean and standard error for each group
mtscf[,mean:=mean(transposition.event.outside.peak),'APOE4_carrier']
mtscf[,se:=sd(transposition.event.outside.peak)/.N,'APOE4_carrier']

# Create the barplot with error bars
ggplot(mtscf, aes(x = APOE4_carrier, y = mean, fill = APOE4_carrier)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Barplot with Standard Error", x = "Group", y = "Mean Value") +
  theme_minimal()

ggplot(mtscf,aes(x=APOE4_carrier,y=transposition.event.outside.peak,col=apoe4))+
  geom_bar(outlier.shape = NA)+ 
  geom_jitter(width=0.4,size=1)+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))

#sample lvel - cell type
# mtd[,avg.pct.read.in.peak:=mean(pct_frag_in_peaksDCT),by=c('libraryID')]
# mtd[,avg.pct.read.in.peak.ct:=mean(pct_frag_in_peaksDCT),by=c('libraryID','cell_type')]
# 
# mtsc<-unique(mtd,by=c('libraryID','cell_type'))
# 
# ggplot(mtsc)+geom_boxplot(aes(x=cell_type,y=avg.pct.read.in.peak.ct,fill=disease))+
#   theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))
# 
# ggplot(mtsc)+geom_boxplot(aes(x=cell_type,y=avg.pct.read.in.peak.ct,fill=apoe4))+
#   theme_bw()+scale_x_discrete(guide = guide_axis(angle = 60))+facet_wrap('disease')
# 
# mtsc[,wilcox.test()]

#stats 
mtd<-fread('outputs/04-ROSMAP_MIT_ATAC/all_final_ATACseq_nuclei_metadata.csv.gz')

ggplot(mtd)+geom_bar(aes(x=libraryID))+facet_wrap('cognitive_status',scales = 'free_x')+theme_bw()
?FeatureMatrix



#BONUS)
#can we recapitulate epigenetic erosion found in the cell paper??

#QUESTIONS####
#Astrocyte cluster analysis
#we have found that Astroyte of APOE4 career have less transposition in peak (open chromatin) region 
#meaning that DNA of APOE4 are more cut in 'heterochromatin' region compared to APOE33 carreer,
#1) is there an APOE4 spe Astrocyte cluster with specific chromatin profile??
#2) can we link them to Astrocyte RNA cluster??
#3) retrotransposon region are more opened??




#10) pseudobulk peak count creation


#11) test a first experiment:
#APOE4vs3 per celltype chrine access
#APOE4 effecton heterochrine disruption
#or APOE TF footprint

