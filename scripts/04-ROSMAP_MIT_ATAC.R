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



#4)cell type level peak calling
#split object per cell type then peak call per groups of donors
#group of donors (to increase number of cells), group according to KNN/SNN celltype spe PCA of celltype spe peak count matrix
#so in any case, need  cellspe peak count matrix first on top20% donors
#then per groups of donors peak calling (group can be different depending of cell type)


#5) pseudobulk peak count creation


#6) test a first experiment

