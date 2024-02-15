out<-'outputs/05-ROSMAP_MIT_WGS/'
dir.create(out)

source('../../../utils/r_utils.R')
genotype_file<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/ref-data/ROSMAP/ROSMAP_NIA_geno/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc'

pipelines<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/xqtl-pipeline/pipeline/'

analysis_name='ROSMAP-snRNA' 

#Functions####

Container<-function(name){
  container_address<-'oras://ghcr.io/cumc/*_apptainer:latest'
  return(str_replace(container_address,'\\*',name))
}


#scRNA - WGS samples overlap
fam<-fread(ps(genotype_file,'.fam'))
biosc<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/Metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv')
bio<-fread('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/ROSMAP_biospecimen_metadata.csv')

length(unique(biosc$individualID)) #447
comm_inds<-intersect(bio[specimenID%in%fam$V2][assay=='wholeGenomeSeq']$individualID,
          unique(biosc$individualID)) #388/447


wgs_tokeep<-bio[specimenID%in%fam$V2][assay=='wholeGenomeSeq'][comm_inds,on='individualID']$specimenID
fwrite(data.table(0,wgs_tokeep),fp(out,'scRNA_WGS_common_samples.tsv'),sep='\t',col.names = F)

#1)keep only samples overlap with scRNA samples
genotype_file_match<-fp(out,ps(basename(genotype_file),'_match_snRNA'))

cmd_matchwgs<-paste('plink2' ,
                 '--bfile' ,genotype_file ,
                 '--keep',fp(out,'scRNA_WGS_common_samples.tsv'),
                 '--make-bed', 
                 '--out' ,genotype_file_match,
                 '--threads' ,'8',
                 '--new-id-max-allele-len','1000',
                 '--set-all-var-ids ','chr@:#_\\$r_\\$a')



CreateJobFile(cmd_matchwgs,file = 'scripts/05A-wgs_match_snRNA.qsub',
              loadBashrc = T,proj_name = 'tcwlab')

RunQsub('scripts/05A-wgs_match_snRNA.qsub',job_name = 'scRNAWGS')


#2)PCA analysis to find global structure of the genotype data (i.e. the population ancestry characteristic)
#2a)need to be done on unrelated samples
#use KING to estimate relatedness between individual, remove the related one.

cmds<-list(kinship=paste('sos run', file.path(pipelines,'GWAS_QC.ipynb'),'king',
                         '--cwd',fp(out,'kinship') ,
                         '--genoFile',genotype_file_match ,
                         '--name', analysis_name,
                         '--no-maximize-unrelated',
                         '--mem 100G',
                         '-s force',
                         '--container', Container('bioinfo')
))

CreateJobFile(cmds,file = 'scripts/05B-kinship.qsub',
              micromamba_env = 'pisces-rabbit',nThreads = 28,loadBashrc = T,proj_name = 'tcwlab-adsp'
)
RunQsub('scripts/05B-kinship.qsub',job_name = 'kinship')


#2b) need to be done on High Quality and unrelated vairants: 
#using missingness > 10% ,MAC > 5 , and LD-prunning in preparation for PCA analysis
#run variants QC
CreateJobForRfile('scripts/05C-pca_qc.R',proj_name = 'tcwlab')
RunQsub('scripts/05C-pca_qc.R',job_name = 'pca_qc',wait_for = 'kinship')

#2c) perform the pca on unrelated samples and variants
#run pca
CreateJobForRfile('scripts/05D-pca.R',proj_name = 'tcwlab')
RunQsub('scripts/05D-pca.R',job_name = 'pca',wait_for = 'pca_qc')

#see pca
unrelated_pca_model<-list.files(out,pattern='prune.pca.rds$',full.names = T)

pca<-readRDS(unrelated_pca_model)

pca_dt<-fread(str_replace(unrelated_pca_model,'rds$','txt'))

ggplot(pca_dt)+geom_point(aes(x=PC1,y=PC2))
ps<-lapply(1.5:10*2,function(i)ggplot(pca_dt)+geom_point(aes_string(x=paste0('PC',i),y=paste0('PC',i+1)),size=0.5))
wrap_plots(ps)

#see screeplot
sum(pca$pca_model$pve)
p1<-ggplot(data.frame(PC=1:length(pca$pca_model$pve),PVE=pca$pca_model$pve),
           aes(x=PC,y=PVE))+geom_point()+geom_line()+
  ggtitle("Scree Plot") 
PVE<-pca$pca_model$pve

# Cumulative PVE plot
p2<-ggplot(data.frame(PC=1:length(pca$pca_model$pve),PVE=pca$pca_model$pve),
           aes(x=PC,y=cumsum(PVE)))+geom_point()+geom_line()+
  ggtitle("Cumulative Scree Plot") 

p1+p2
ggsave(fp(out,'scree_plot.png'))




