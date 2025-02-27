#Check pseudobulk ATAC peaks before ftp transfer
out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
out1<-fp(out,'pseudobulk_data')
out2<-fp(out1,'analysis_ready')
dir.create(out2)

mtd<-fread(fp(out,'all_final_ATACseq_nuclei_metadata.csv.gz'))
table(mtd[main_cell_type=='Astro']$individualID)
mtct<-fread('outputs/04-ROSMAP_MIT_ATAC/all_final_ATACseq_nuclei_main_cell_type_level_metadata.csv.gz')
table(unique(mtct,by='libraryID')$sequencingBatch)

mtd_bio<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv')
intersect(mtd_bio[assay=='snMultiome']$specimenID,
          mtct$specimenID) #only snATAC

mtctf<-mtct[,.(libraryID,individualID,sequencingBatch,main_cell_type,n.ct,avg.pct.read.in.peak.ct,
               med.nucleosome_signal.ct,med.n_tot_fragment.ct, med.tss.enrich.ct)] 
mtctf[,.SD[.N>1],by=.(individualID,main_cell_type)]
#need merge info by donor instead of library
mtctf[,n.nuclei:=sum(n.ct),by=.(individualID,main_cell_type)]
mtctf[,avg.pct.read.in.peak.ct:=sum(avg.pct.read.in.peak.ct*n.ct/n.nuclei),
      by=.(individualID,main_cell_type)]
mtctf[,med.nucleosome_signal.ct:=sum(med.nucleosome_signal.ct*n.ct/n.nuclei),
      by=.(individualID,main_cell_type)]
mtctf[,med.n_tot_fragment.ct:=sum(med.n_tot_fragment.ct*n.ct/n.nuclei),
      by=.(individualID,main_cell_type)]
mtctf[,med.tss.enrich.ct:=sum(med.tss.enrich.ct*n.ct/n.nuclei),
      by=.(individualID,main_cell_type)]
mtctf<-unique(mtctf[,-c('libraryID','n.ct')],by=c('individualID','main_cell_type'))

#n peaks n samples
pseudo<-readRDS('outputs/04-ROSMAP_MIT_ATAC/pseudobulk_data/Astro.rds')
pseudo[[]]

dim(pseudo) #531489     90
#QC to do: keep donors only if with >50 nuclei
pseudo$n.nuclei<-mtctf[main_cell_type=='Astro'][colnames(pseudo),on='individualID']$n.nuclei

pseudof<-subset(pseudo,n.nuclei>50)
dim(pseudof)
fwrite(data.table(as.matrix(pseudof@assays$peaksDCT@counts),keep.rownames = 'peak_id'),
       fp(out2,'pseudobulk_peaks_counts_Astro_50nuc.csv.gz'))
#save mtd
fwrite(mtctf[main_cell_type=='Astro'][colnames(pseudof),on='individualID'],
       fp(out2,'metadata_Astro_50nuc.csv'))

#for all

files<-list.files('outputs/04-ROSMAP_MIT_ATAC/pseudobulk_data/',pattern = '.rds',full.names = T)
for(f in files){
  ct<-basename(f)|>str_remove('.rds$')
  print(ct)
  pseudo<-readRDS(f)
  
  #QC to do: keep donors only if with >50 nuclei
  pseudo$n.nuclei<-mtctf[main_cell_type==ct][colnames(pseudo),on='individualID']$n.nuclei
  
  pseudof<-subset(pseudo,n.nuclei>50)
  print(dim(pseudof))
  fwrite(data.table(as.matrix(pseudof@assays$peaksDCT@counts),keep.rownames = 'peak_id'),
         fp(out2,ps('pseudobulk_peaks_counts',ct,'_50nuc.csv.gz')))
         #save mtd
  fwrite(mtctf[main_cell_type==ct][colnames(pseudof),on='individualID'],
         fp(out2,ps('metadata_',ct,'_50nuc.csv')))
}

#n sample per celltype
ggplot(mtctf[n.nuclei>50])+
  geom_bar(aes(x=main_cell_type,
               fill = main_cell_type))+
  theme_bw()+labs(y='# individuals')

#put in ftp

#sanity check
mtctf[main_cell_type=='Astro'],


#samples removed ?
#outlier: outlier.clinical.status(apoe4&apoe2|African)|outlier.cellprop(PC1 distrib>3*IQR)|oligo.donor.outlier|exc.donor.outlier
mts<-fread(fp(out,'all_final_ATACseq_nuclei_sample_level_metadata.csv.gz'))
mts[(donor.outlier)]
unique(mts$individualID)
setdiff(mts$individualID,colnames(pseudo))

mts[individualID=='R4262244']
mts[(outlier.ethnicity)]


#Cells removed?
#outlier: nucleosome_signal<3*IQR|TSS.enrichment<IQR |nFeature_peaksCT>3*IQR