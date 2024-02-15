
out<-'outputs/05-ROSMAP_MIT_WGS/'
dir.create(out)

source('../../utils/r_utils.R')
genotype_file<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/ref-data/ROSMAP/ROSMAP_NIA_geno/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc'
genotype_file_match<-fp(out,ps(basename(genotype_file),'_match_snRNA.bed'))

pipelines<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/xqtl-pipeline/pipeline/'

analysis_name='ROSMAP-snRNA' 

related_file<-list.files(file.path(out,'kinship'),pattern='\\.related\\.fam$',full.names = T)
having.related<-length(related_file)>0


#Functions####

Container<-function(name){
  container_address<-'oras://ghcr.io/cumc/*_apptainer:latest'
  return(str_replace(container_address,'\\*',name))
}

#ANALYSIS####

prunebed_file<-list.files(file.path(out,'cache'),pattern='plink_qc\\.prune\\.bed$',full.names = T)

cmd_pca<-paste('sos run', file.path(pipelines,'PCA.ipynb'),'flashpca',
               '--cwd',out ,
               '--genoFile',prunebed_file ,
               '--mem 80G',
               '--container', Container('flashpcar')
)

CreateJobFile(cmd_pca,file = 'scripts/05D-pca_unrelated.qsub',
              micromamba_env = 'pisces-rabbit',nThreads = 16,proj_name = 'tcwlab-adsp'
)
jobid<-RunQsub('scripts/05D-pca_unrelated.qsub',job_name = 'pca_unrel',wait_for = 'pca_qc')
WaitQsub('scripts/05D-pca_unrelated.qsub',jobid)

unrelated_pca_model<-list.files(out,pattern='prune.pca.rds$',full.names = T)

#project related individuals
#project back related individuals to the space + detect sample outliers based on PCs
#decide to use the 7 first PCs to compute outliers (the --maha-k parameter), based on the PC scree plot. 
if(having.related){
  
  related_pca_bed<-list.files(file.path(out,'cache'),pattern='extracted.bed$',full.names = T)
  
  # prepared also related ind to pca by keeping same SNPs than unrelated samples
  cmd_proj_rel<-paste('sos run',file.path(pipelines,'PCA.ipynb'), 'project_samples',
                      '--cwd', out,
                      '--genoFile', related_pca_bed,
                      '--pca-model',  unrelated_pca_model, 
                      '--maha-k' ,'2',
                      '--container', Container('flashpcar'))
  
  
  CreateJobFile(cmd_proj_rel,file = 'scripts/05D-proj_related.qsub',
                micromamba_env = 'pisces-rabbit',nThreads = 16,proj_name = 'tcwlab-adsp'
  )
  RunQsub('scripts/05D-proj_related.qsub',job_name = 'proj_related',wait_for = 'pca_unrel')
  
  
}
