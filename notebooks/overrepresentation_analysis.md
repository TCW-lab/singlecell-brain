# Over representation test

Here we will describe how to perform over-representation analysis (ORA). ORA test the stastical significance of an overlap between a list of genes (your gene list of interest, eg the DEGs) and reference gene-sets (e.g. of known biological pathways) using fisher's exact test/ hypergeometric test. With GSEA used by fgsea, ORA is the most common to test use when performing pathway/functional enrichment analysis. The difference with GSEA is are mainly your input data :  ORA use a list of genes of interest (e.g DEGs at certain threshold, genes of a certain module) while GSEA used a ranked list of all your genes tested in an experiment.
Depending of your biological question, one of the others can be more relevant to use, or both can be used also as a validation.

For this example we will used our downladed [MSIGDB pathways/ontologies](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) as reference gene-sets to test functional enrichment in a module of co-expressed genes identified through WGCNA analysis.

## Set up 
We will use the prebuilt function `OR3()` to perform the Over-representation test. this function is available in [alexandre-utils/r_utils.R](https://github.com/TCW-lab/alexandre-utils/blob/main/r_utils.R#L399-L450).
we load the MSIGDB gene sets and filter for pathways/genesets with a too small/high size.



```R
source('../../../utils/r_utils.R')

pathways<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

pathways_infos<-fread('/projectnb/tcwlab/MSigDB/all_CPandGOs_genesets_metadata.csv.gz')

#test pathways <2k et >5
pathwaysf<-pathways[pathway.size>5&pathway.size<2000]
length(unique(pathwaysf$pathway))#12k



```

## Load or query gene list of interest


```R
#Neurons
mods_dt<-fread('../../APOE4_proteomics/outputs/02B-module_conservation/Protein_modules/neurons_modules.csv')

#rm non annotated or unassigned genes
mods_dtf<-mods_dt[!(is.na(gname)|gname=='')&module!='grey'] #grey is for unassigned genes

```

## define the gene background
One important consideration when performing ORA is to know in which 'background' the test is made. Indeed, we need to know what are the genes that have been considered in our study/experiment. For RNA-seq analysis, this is often all the genes that have been tested in your differential expression analysis. Here because we are performing function enrichment in WGCNA modules, we will used as background all genes that have been considered in the WGCNA analysis, ie. all the genes that pass the QC step. Here this correspond to all genes assigned to a module plus those unassigned which are therefore located in the 'grey' module


```R
gene_background=mods_dt$gname
```

## perform the over-representation test for module Red genes
Let's do it for one module, e.gg the module 'red'
The function have 3 required parameter :
- `querys` : vector of your genes of interest. if is a list of several genes list (of several vector in R language), will perform the test for every vector.
- `terms_list`: named list of gene-sets (terms_list) used as reference pathways,
- `background` : vector of gene names background


```R
red_genes=mods_dt[module=='red']$gname

pathways_list=split(pathwaysf$gene,pathwaysf$pathway)

res_red_enr<-OR3(querys=red_genes,terms_list = pathways_list,background =gene_background)

```

we can assess pathways enriched in this module


```R
res_red_enr[padj<0.05]
```

## perform the over-representation test for all modules
Here we will perform perform the test for each modules


```R
modules_list=split(mods_dtf$gname,mods_dtf$module)
res_modules_enr<-OR3(modules_list,
                     terms_list = pathways_list,background =gene_background)


#number of pathways enriched per module
table(res_modules_enr[padj<0.05]$query)


```

## perform the over-representation test separating each category of gene sets
12000 gene set are tested together but are coming from different source (e.g. GO Biological process, KEGG..).  

Just for multiple correction test purpose (i.e `padj`calculation)This is better to perform the OR3 test independantly for each of these sources. To do that we can split first the pathways data.frame per source (`subcat`) and then perform OR3


```R
splitted_pathways<-split(pathwaysf,by='subcat')
res_modules_enr<-rbindlist(lapply(splitted_pathways,function(pathwf)OR3(modules_list,
                                                                        terms_list = split(pathwf$gene,pathwf$pathway),
                                                                        background =gene_background)))

#number of pathways enriched per module
table(res_modules_enr[padj<0.05]$query)
                                  
```


```R
#add subcategory and pathway size info
setnames(pathways_infos,old = 'pathway','term') #for column name compatibility allowing merging the data.frames

res_modules_enr<-merge(res_modules_enr,unique(pathways_infos,by='term'),by='term')[order(query,term,pval)]

```

We can then save these results using e.g. `fwrite` and vizualize it using e.g. `ggplot2` 


```R

```
