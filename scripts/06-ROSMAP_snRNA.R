
#check fgsea results
source('../../utils/r_utils.R')
res<-readRDS('../../../acandib/Project/APOE_Ab/Pseudobulk_QC/DESeq2_Results/fgsea_Results.RDS')

res_gsea<-res[["Astrocytes_Pseudobulk"]][["AD"]][["4-3"]][padj<0.05][order(padj)]


#annotate fgsea results
pathways_info<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')
res_gsea<-merge(res_gsea,pathways_info,by=c('pathway'))
table(res_gsea[padj<0.05])

#table
ggplot(res_gsea[padj<0.05])+geom_bar(aes(x=NES>0,fill=subcat))+
  labs(y='number of enriched terms (padj<0.05)')

#toppathways by cell type
source('../../utils/visualisation.R')
res_gsea_top30<-res_gsea[padj<0.05,.SD[head(order(padj),40)],
                             by=.(subcat)]

emmaplot(head(res_gsea_top30[order(pval)],100),label.size=2,cols=c('blue','white','red'))

#only down
emmaplot(head(res_gsea[padj<0.05][order(pval)][NES<0],100),label.size=2,cols=c('blue','white','red'))



#autophagy 
res_gsea[padj<0.05][str_detect(pathway,'AUTO')]
emmaplot(head(res_gsea_top30[cell_type=='Astrocyte'][order(pval)],100),label.size=2,cols=c('blue','white','red'))
#rep to eostradiol..because sex bias ?
table(mtd[ipsc_type=='population'&cell_type=='Astrocyte']$genotype,mtd[ipsc_type=='population'&cell_type=='Astrocyte']$gender)
