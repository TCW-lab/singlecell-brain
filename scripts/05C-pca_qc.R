
out<-'outputs/05-ROSMAP_MIT_WGS/'
dir.create(out)

source('../../utils/r_utils.R')
genotype_file<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/ref-data/ROSMAP/ROSMAP_NIA_geno/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc'
genotype_file_match<-fp(out,ps(basename(genotype_file),'_match_snRNA.bed'))

pipelines<-'/projectnb/tcwlab-adsp/member/adpelle1/projects/fungen-xqtl/xqtl-pipeline/pipeline/'

analysis_name='ROSMAP-snRNA' 

#Functions####

Container<-function(name){
  container_address<-'oras://ghcr.io/cumc/*_apptainer:latest'
  return(str_replace(container_address,'\\*',name))
}

#ANALYSIS####

related_file<-list.files(file.path(out,'kinship'),pattern='\\.related\\.fam$',full.names = T)
having.related<-length(related_file)>0

unrelated_file<-list.files(file.path(out,'kinship'),pattern='unrelated\\.bed$',full.names = T)

if(having.related){
  
  cmd<-paste('sos run', file.path(pipelines,'GWAS_QC.ipynb'),'qc',
             '--cwd',fp(out,'cache') ,
             '--genoFile',unrelated_file,
             '--name', analysis_name ,
             '--mac-filter 5',
             '--mem 80G',
             '-s force',
             '--container', Container('bioinfo')
  )
  
  
}else{
  system(paste('rm',paste(paste0(tools::file_path_sans_ext(unrelated_file),c('.bed','.bim','.fam')),collapse = ' ')))
  
  cmd<-paste('sos run', file.path(pipelines,'GWAS_QC.ipynb'),'qc',
             '--cwd',fp(out,'cache') ,
             '--genoFile',genotype_file_match ,
             '--name', analysis_name ,
             '--mac-filter 5',
             '--mem 80G',
             '-s force',
             '--container', Container('bioinfo')
  )
  
  
}

CreateJobFile(cmd,file = 'scripts/05C-qc_unrelated.qsub',
              micromamba_env = 'pisces-rabbit',nThreads = 16,proj_name = 'tcwlab-adsp'
)
RunQsub('scripts/05C-qc_unrelated.qsub',job_name = 'QCunrel')


if(having.related){
  related_id<-list.files(file.path(out,'kinship'),pattern='\\.related_id$',full.names = T)
  related_file<-str_replace(related_id,'related_id$','related.bed')
  
  # prepared also related ind to pca by keeping same SNPs than unrelated samples
  prune_file<-list.files(file.path(out,'cache'),pattern='.prune.in$',full.names = T)
  
  #create plink file because King step did not create correct file
  
  cmd_prune<-paste('sos run ../xqtl-pipeline/pipeline/GWAS_QC.ipynb qc_no_prune',
                   '--cwd',fp(out,'cache') ,
                   '--genoFile',related_file,
                   '--maf-filter 0',
                   '--geno-filter 0',
                   '--mind-filter 0.1',
                   '--hwe-filter 0',
                   '--keep-variants',prune_file,
                   '--mem 16G',
                   '-s force',
                   '--container', Container('bioinfo'))
  
  CreateJobFile(cmd_list = list('create_related_file'=cmd_relat,
                                'prune_related'=cmd_prune),
                file = 'scripts/05C-prune_related.qsub',
                micromamba_env = 'pisces-rabbit',nThreads = 16
  )
  
  RunQsub('scripts/05C-prune_related.qsub',job_name = 'PruneRelated',wait_for = 'QCunrel')
  
}
