out<-'outputs/01-SEAAD_data/fgsea'
dir.create(out)

source('../../utils/r_utils.R')
library(fgsea)


pathways_dir <- "/projectnb/tcwlab/MSigDB/"
pathways_to_test<-c("GO_all", 'CP_all') # choose pathway gene sets reference you want to test or let like that to test for all
gmt_mtd<-fread(file.path(pathways_dir,'gmt_metadata.csv')) 
pathways_info<-fread(file.path(pathways_dir,'all_CPandGOs_genesets_metadata.csv.gz'))

#load data
res_de<-rbindlist(list(fread('outputs/01-SEAAD_data/pseudobulk_deseq2/res_pseudobulkDESeq2_No.dementia_APOE4_vs_3_astro_DLPFC_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz'),
                   fread('outputs/01-SEAAD_data/pseudobulk_deseq2/res_pseudobulkDESeq2_No.dementia_APOE4_vs_3_astro_MTG_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz')))
res_de
#Prep data
#duplicate gene names?
res_de[,is.dup:=duplicated(gene),by='brain_region']
res_de[,n.dup:=sum(is.dup),by=.(brain_region,gene)]
table(res_de$n.dup)
#    0    1    5    8 
# 6220   14    6    9 
#no need 
# #ok just take the best
res_de[order(brain_region,-abs(stat)),to_rm:=duplicated(gene),by='brain_region']
res_def<-res_de[!(to_rm)]


res_gsea_all<-Reduce(rbind,lapply(unique(res_def$brain_region), function(br){
  
  message(paste('testing enrichement in',br))
  
  res_de<-res_def[brain_region==br]
  #calculate/extract the gene stats
  gene_stats<-setNames(res_de$stat,res_de$gene)
  
  #run fgsea for every pathway source
  res_gsea<-Reduce(rbind,lapply(pathways_to_test, function(p){
    
    message(paste('testing enrichement for',gmt_mtd[name==p]$desc))
    
    pathways<- gmtPathways(file.path(pathways_dir,gmt_mtd[name==p]$gmt))
    
    res<-fgsea(pathways,
               stats=gene_stats,minSize=10,maxSize=2000,scoreType='std',nPermSimple = 10000)
    
    return(res[,source:=p])
    
  }))
  return(res_gsea[,brain_region:=br])
  
  
  
}))

res_gsea_all<-merge(res_gsea_all,pathways_info,by=c('pathway'))[order(brain_region,source,subcat,pval)]
table(res_gsea_all[padj<0.05]$brain_region)
# Astrocyte       MCC 
#        60       482 
fwrite(res_gsea_all,fp(out,"res_fgsea_APOE4vs3_NL_astrocytes.csv.gz"))

res_gsea_all[padj<0.05][str_detect(pathway,'DNA_REPAIR')]
source('../../utils/emmaplot.R')
emmaplot(head(res_gsea_all[padj<0.05][brain_region=='MTG'][order(padj)],100))
ggsave(fp(out,"res_fgsea_APOE4vs3_NL_astrocytes.pdf"),width=10,height=8)
