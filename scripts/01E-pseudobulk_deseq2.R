out<-'outputs/01-SEAAD_data/pseudobulk_deseq2'
dir.create(out)
source('../../utils/r_utils.R')
source('../../utils/pca_utils.R')
library(DESeq2)
library(edgeR)
library(future)
library(sctransform)

#example for Astrocyte
#DLPFC####
#load data 
pseudo_count<-fread('outputs/01-SEAAD_data/DLPFC/pseudobulk_main_cell_type/Astro.csv.gz')

mat<-as.matrix(data.frame(pseudo_count,row.names = 'gene_id'))

mtd<-fread('outputs/01-SEAAD_data/DLPFC/pseudobulk_main_cell_type/all_final_RNAseq_nuclei_main_cell_type_level_metadata.csv.gz')[main_cell_type=='Astro']


#QC
#remove samples if flagged as outliers (See QC script) 
#here, all donors with clinical or cellular abnormalities are exclude, and if the pseudobulk come from n.cell<50
mtd[(outlier)][,.(Donor.ID,Cognitive.Status,outlier.clinical.status,ethnicity,outlier.ethnicity,cell_prop_dev_pc1,outlier.cellprop,n.cells,outlier.n.cells)]
mtdf<-mtd[!(outlier|outlier.n.cells)]
nrow(mtdf)#71

to_keep<-mtdf$Donor.ID 
mat<-mat[,to_keep]

#genes signicantly express (thr= 1CPM) in less than 10% of sample are removed
isexpr <- rowSums(cpm(mat)>1) >= 0.1 * ncol(mat)
sum(isexpr)#21k
matf <- mat[isexpr,]


#influence of covariates
##first have a look on corelation between our factor of interest (here apoe4) and other covariates
#let's first transform as numerical all covariates which can be
mtdf[,braak.num:=as.numeric(as.factor(Braak))]
mtdf[,.(braak.num,Braak)]
mtdf[,thal_score:=as.numeric(factor(Thal))]
mtdf[,.(Thal,thal_score)]
mtdf[,atherosclerosis:=as.numeric(factor(Atherosclerosis,levels = c('None','Mild','Moderate','Severe')))]

#also levels correctly categorical factor with the good reference
mtdf[,Cognitive.Status:=factor(Cognitive.Status,levels = c('No dementia','Dementia'))] #reference is 'No dementia', so put in first positon,
mtdf[,APOE4.Status:=factor(APOE4.Status,levels = c('N','Y'))]
mtdf[,Sex:=factor(Sex,levels = c('Female','Male'))]

covs_to_check<-list(categorical=c('batch_vendor_name','method',
                                  'Cognitive.Status','Sex','APOE4.Status'),
                    nums=c('avg.pct.mt.ct','n.cells','med.umis.per.cell',
                           'med.genes.per.cell','library_input_ng','pcr_cycles',
                           'PMI','age_at_death_num','Brain.pH','braak.num',
                           'Fresh.Brain.Weight','thal_score','atherosclerosis'))

factor_of_int<-'APOE4.Status'

res_cor_fctr<-rbindlist(lapply(setdiff(covs_to_check$nums,factor_of_int), function(f){
  mod<-lm(unlist(mtdf[,lapply(.SD, as.numeric),.SDcols=factor_of_int])~unlist(mtdf[,..f]))
  summstats<-summary(mod)
  data.table(factor=f,
             p=summstats$coefficients[2,4],
             beta=summstats$coefficients[2,1],
             R2=summstats$adj.r.squared)
}))
res_cor_fctr[p<0.05] #asso with braak.num, thal score, which is expected


#we then perform a PCA to assess influence of each covariates in the main transcriptome variance
#based on that, we could choose for which technical covariates or putative confounding bioligical factor to correct for
#scale the data using a variance stabilizing transformation (vst)to UMI count data using a regularized Negative Binomial regression model.
#plan(strategy = "multicore", workers = 4)

matf_norm_scaled <- sctransform::vst(matf)$y #y is pearson residual, so already normalized data


pca<-RunPca(matf_norm_scaled,scale = F,center = F)



for(i in seq.int(from = 1,to = length(unlist(covs_to_check)),by=4)){
  covs_to_plot<-head(unlist(covs_to_check)[i:length(unlist(covs_to_check))],4)
  PcaPlot(pca,mtdf,
          group.by =covs_to_plot ,
          sample_col = 'Donor.ID',ncol=2)
  
  ggsave(fp(out,ps('pca_astrocyte_pseudobulk_covs_',paste(covs_to_plot,collapse = '_'),'.png')),width = 7,height = 6)
  
}

res_cor_pcs<-CorrelCovarPCs(pca,mtdf,
                        sample_col = 'Donor.ID',
               vars_num = covs_to_check$nums,
               vars_cat = covs_to_check$categorical)

plotPvalsHeatMap(res_cor_pcs,p.thr = 0.1,p_col='padj',
                 labels_col =paste0(paste0('PC',1:10),"(",round(pctPC(pca,rngPCs = 1:10)*100,0),"%)") )



#avg.mt highly assoc with PC1 but very visible in the PcaPlot(), more visible if we logtransform data?  
mtdf_pcs<-merge(data.table(pca$x,keep.rownames = 'Donor.ID'),mtdf[,.SD,.SDcols=unique(colnames(mtdf))])
ggplot(mtdf_pcs)+geom_point(aes(x=PC1,y=PC2,col=log(avg.pct.mt.ct)))+theme_bw()#yes

## choose covariates to include
#choose to correct for putative technical differences : avg.pct.mt.ct, PMI, library input ng, pcr_cycles, n.cells, med.umis.cell
#choose to correct for putative biological confounder : sex, atherosclerosis
design= ~avg.pct.mt.ct+PMI+library_input_ng+pcr_cycles+n.cells+med.umis.per.cell+atherosclerosis+Sex+APOE4.Status

#scale numerical covariates to improve GLM convergence
covs_to_scale<-c('PMI','avg.pct.mt.ct','library_input_ng','pcr_cycles','n.cells','med.umis.per.cell')
mtdf_scaled<-copy(mtdf)
mtdf_scaled[,(covs_to_scale):=lapply(.SD,scale),.SDcols=covs_to_scale]

#run DESEQ2
dds <- DESeqDataSetFromMatrix(matf, 
                              colData = data.frame(mtdf_scaled,row.names="Donor.ID")[colnames(matf),], 
                              design = design)
dds <- DESeq(dds)

#extract contrast of interest
mod_mat <- model.matrix(design(dds), colData(dds))

apoe4 <- colMeans(mod_mat[dds$APOE4.Status== "Y", ])
apoe3 <- colMeans(mod_mat[dds$APOE4.Status == "N", ])

res <- results(dds,contrast = apoe4-apoe3,alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

#add some annotation to the results
res[,cell_type:="Astro"][,brain_region:='DLPFC']
res[padj<0.25][order(padj)]

res[,design:=paste0('~',as.character(design)[2])]

fwrite(res,fp(out,"res_pseudobulkDESeq2_APOE4_vs_3_astro_DLPFC_Cov_sex_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz"))
res<-fread(fp(out,"res_pseudobulkDESeq2_APOE4_vs_3_astro_DLPFC_Cov_sex_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz"))

#test correcting for braak stage
design= ~avg.pct.mt.ct+PMI+library_input_ng+pcr_cycles+n.cells+med.umis.per.cell+braak.num+atherosclerosis+Sex+APOE4.Status

dds <- DESeqDataSetFromMatrix(matf, 
                              colData = data.frame(mtdf_scaled,row.names="Donor.ID")[colnames(matf),], 
                              design = design)
dds <- DESeq(dds)

#extract contrast of interest
mod_mat <- model.matrix(design(dds), colData(dds))

apoe4 <- colMeans(mod_mat[dds$APOE4.Status== "Y", ])
apoe3 <- colMeans(mod_mat[dds$APOE4.Status == "N", ])

res_braak<- results(dds,contrast = apoe4-apoe3,alpha = 0.05)

res_braak<-data.table(as.data.frame(res_braak),keep.rownames="gene")

#add some annotation to the results
res_braak[,cell_type:="Astro"][,brain_region:='DLPFC']
res_braak[padj<0.25][order(padj)]

res_braak[,design:=paste0('~',as.character(design)[2])]

#compa with previous model
res_merge<-merge(res,res_braak,by=c('gene','cell_type','brain_region'))
ggplot(res_merge)+geom_point(aes(x=stat.x,y=stat.y)) 

#very similar results, but some are better when modeling for braak stage
#we want effect of APOE4 at 'early' level, so braak stage which denote of late neuropathological markers, make sens

fwrite(res_braak,fp(out,"res_pseudobulkDESeq2_Dementia_APOE4_vs_3_astro_DLPFC_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz"))

#do within control, dementia people (wihtout braak correction this time)
design= ~avg.pct.mt.ct+PMI+library_input_ng+pcr_cycles+n.cells+med.umis.per.cell+atherosclerosis+Sex+APOE4.Status

for(d in c('Dementia','No dementia')){
  mtdf_scaledf<-mtdf_scaled[Cognitive.Status==d]
  matff<-matf[,mtdf_scaledf$Donor.ID]
  dds <- DESeqDataSetFromMatrix(matff, 
                                colData = data.frame(mtdf_scaledf,row.names="Donor.ID")[colnames(matff),], 
                                design = design)
  dds <- DESeq(dds)
  
  #extract contrast of interest
  mod_mat <- model.matrix(design(dds), colData(dds))
  
  apoe4 <- colMeans(mod_mat[dds$APOE4.Status== "Y", ])
  apoe3 <- colMeans(mod_mat[dds$APOE4.Status == "N", ])
  
  res<- results(dds,contrast = apoe4-apoe3,alpha = 0.05)
  
  res<-data.table(as.data.frame(res),keep.rownames="gene")
  
  #add some annotation to the results
  res[,cell_type:="Astro"][,brain_region:='DLPFC'][,Cognitive.Status:=d]
  res[padj<0.25][order(padj)]
  
  res[,design:=paste0('~',as.character(design)[2])]
  

  fwrite(res,fp(out,paste0("res_pseudobulkDESeq2_",make.names(d),"_APOE4_vs_3_astro_DLPFC_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz")))
  
}

#for MTG####
#load data 
pseudo_count<-fread('outputs/01-SEAAD_data/MTG/pseudobulk_main_cell_type/Astro.csv.gz')

mat<-as.matrix(data.frame(pseudo_count,row.names = 'gene_id'))

mtd<-fread('outputs/01-SEAAD_data/MTG/pseudobulk_main_cell_type/all_final_RNAseq_nuclei_main_cell_type_level_metadata.csv.gz')[main_cell_type=='Astro']

#QC
#remove samples if flagged as outliers (See QC script) 
#here, all donors with clinical or cellular abnormalities are exclude, and if the pseudobulk come from n.cell<50
mtd[(outlier)][,.(Donor.ID,Cognitive.Status,outlier.clinical.status,ethnicity,outlier.ethnicity,cell_prop_dev_pc1,outlier.cellprop,n.cells,outlier.n.cells)]
mtdf<-mtd[!(outlier|outlier.n.cells)]
nrow(mtdf)#77

to_keep<-mtdf$Donor.ID 
mat<-mat[,to_keep]

#Filter lowly expressed genes
nrow(mat)
#genes signicantly express (thr= 1CPM) in less than 10% of sample are removed
isexpr <- rowSums(cpm(mat)>1) >= 0.1 * ncol(mat)
sum(isexpr)#21k
matf <- mat[isexpr,]


#influence of covariates
##first have a look on corelation between our factor of interest (here apoe4) and other covariates
#let's first transform as numerical all covariates which can be
mtdf[,braak.num:=as.numeric(as.factor(Braak))]
mtdf[,.(braak.num,Braak)]
mtdf[,thal_score:=as.numeric(factor(Thal))]
mtdf[,.(Thal,thal_score)]
mtdf[,atherosclerosis:=as.numeric(factor(Atherosclerosis,levels = c('None','Mild','Moderate','Severe')))]

#also levels correctly categorical factor with the good reference
mtdf[,Cognitive.Status:=factor(Cognitive.Status,levels = c('No dementia','Dementia'))] #reference is 'No dementia', so put in first positon,
mtdf[,APOE4.Status:=factor(APOE4.Status,levels = c('N','Y'))]
mtdf[,Sex:=factor(Sex,levels = c('Female','Male'))]

covs_to_check<-list(categorical=c('batch_vendor_name','method',
                                  'Cognitive.Status','Sex','APOE4.Status'),
                    nums=c('avg.pct.mt.ct','n.cells','med.umis.per.cell',
                           'med.genes.per.cell','library_input_ng','pcr_cycles',
                           'PMI','age_at_death_num','Brain.pH','braak.num',
                           'Fresh.Brain.Weight','thal_score','atherosclerosis'))

factor_of_int<-'APOE4.Status'

res_cor_fctr<-rbindlist(lapply(setdiff(covs_to_check$nums,factor_of_int), function(f){
  mod<-lm(unlist(mtdf[,lapply(.SD, as.numeric),.SDcols=factor_of_int])~unlist(mtdf[,..f]))
  summstats<-summary(mod)
  data.table(factor=f,
             p=summstats$coefficients[2,4],
             beta=summstats$coefficients[2,1],
             R2=summstats$adj.r.squared)
}))
res_cor_fctr[p<0.05] #asso with braak.num, thal score, which is expected


#we then perform a PCA to assess influence of each covariates in the main transcriptome variance
#based on that, we could choose for which technical covariates or putative confounding bioligical factor to correct for
#scale the data using a variance stabilizing transformation (vst)to UMI count data using a regularized Negative Binomial regression model.
#plan(strategy = "multicore", workers = 4)

matf_norm_scaled <- sctransform::vst(matf)$y #y is pearson residual, so already normalized data


pca<-RunPca(matf_norm_scaled,scale = F,center = F)



for(i in seq.int(from = 1,to = length(unlist(covs_to_check)),by=4)){
  covs_to_plot<-head(unlist(covs_to_check)[i:length(unlist(covs_to_check))],4)
  PcaPlot(pca,mtdf,
          group.by =covs_to_plot ,
          sample_col = 'Donor.ID',ncol=2)
  
  ggsave(fp(out,ps('pca_astrocyte_MTG_pseudobulk_covs_',paste(covs_to_plot,collapse = '_'),'.png')),width = 7,height = 6)
  
}

res_cor_pcs<-CorrelCovarPCs(pca,mtdf,
                            sample_col = 'Donor.ID',
                            vars_num = covs_to_check$nums,
                            vars_cat = covs_to_check$categorical)

plotPvalsHeatMap(res_cor_pcs,p.thr = 0.1,p_col='padj',
                 labels_col =paste0(paste0('PC',1:10),"(",round(pctPC(pca,rngPCs = 1:10)*100,0),"%)") )

#interestingly, this time Cognitive status is associated with main source of variance (PC1)


## choose covariates to include
#choose to correct for putative technical differences : avg.pct.mt.ct, PMI, library input ng, pcr_cycles, n.cells, med.umis.cell
#choose to correct for putative biological confounder : sex, atherosclerosis, braak
design= ~avg.pct.mt.ct+PMI+library_input_ng+pcr_cycles+n.cells+med.umis.per.cell+braak.num+atherosclerosis+Sex+APOE4.Status

#scale numerical covariates to improve GLM convergence
covs_to_scale<-c('PMI','avg.pct.mt.ct','library_input_ng','pcr_cycles','n.cells','med.umis.per.cell')
mtdf_scaled<-copy(mtdf)
mtdf_scaled[,(covs_to_scale):=lapply(.SD,scale),.SDcols=covs_to_scale]

#run DESEQ2
dds <- DESeqDataSetFromMatrix(matf, 
                              colData = data.frame(mtdf_scaled,row.names="Donor.ID")[colnames(matf),], 
                              design = design)
dds <- DESeq(dds)

#extract contrast of interest
mod_mat <- model.matrix(design(dds), colData(dds))

apoe4 <- colMeans(mod_mat[dds$APOE4.Status== "Y", ])
apoe3 <- colMeans(mod_mat[dds$APOE4.Status == "N", ])

res <- results(dds,contrast = apoe4-apoe3,alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

#add some annotation to the results
res[,cell_type:="Astro"][,brain_region:='MTG']
res[padj<0.25][order(padj)]

res[,design:=paste0('~',as.character(design)[2])]

fwrite(res,fp(out,"res_pseudobulkDESeq2_APOE4_vs_3_astro_MTG_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz"))


#do within control or dementia people (wihtout braak correction this time)
design= ~avg.pct.mt.ct+PMI+library_input_ng+pcr_cycles+n.cells+med.umis.per.cell+atherosclerosis+Sex+APOE4.Status

for(d in c('Dementia','No dementia')){
  mtdf_scaledf<-mtdf_scaled[Cognitive.Status==d]
  matff<-matf[,mtdf_scaledf$Donor.ID]
  dds <- DESeqDataSetFromMatrix(matff, 
                                colData = data.frame(mtdf_scaledf,row.names="Donor.ID")[colnames(matff),], 
                                design = design)
  dds <- DESeq(dds)
  
  #extract contrast of interest
  mod_mat <- model.matrix(design(dds), colData(dds))
  
  apoe4 <- colMeans(mod_mat[dds$APOE4.Status== "Y", ])
  apoe3 <- colMeans(mod_mat[dds$APOE4.Status == "N", ])
  
  res<- results(dds,contrast = apoe4-apoe3,alpha = 0.05)
  
  res<-data.table(as.data.frame(res),keep.rownames="gene")
  
  #add some annotation to the results
  res[,cell_type:="Astro"][,brain_region:='MTG'][,Cognitive.Status:=d]
  res[padj<0.25][order(padj)]
  
  res[,design:=paste0('~',as.character(design)[2])]
  
  
  fwrite(res,fp(out,paste0("res_pseudobulkDESeq2_",make.names(d),"_APOE4_vs_3_astro_MTG_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz")))
  
}

#compare DLPFC and MTG
#full data
res_merge<-merge(fread(fp(out,"res_pseudobulkDESeq2_APOE4_vs_3_astro_DLPFC_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz")),
                       fread(fp(out,"res_pseudobulkDESeq2_APOE4_vs_3_astro_MTG_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz")),
                 by=c('gene','cell_type'),suffixes = c(".DLPFC", ".MTG"))


ggplot(res_merge)+geom_point(aes(x=stat.DLPFC,y=stat.MTG),size=0.2) +theme_bw() 

res_merge[,summary(lm(stat.DLPFC~stat.MTG))]
#correl but r2 low
#R2= 0.07339 
#p < 2e-16

#by cognitive status
res_merge<-merge(rbind(fread(fp(out,"res_pseudobulkDESeq2_Dementia_APOE4_vs_3_astro_DLPFC_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz")),
                       fread(fp(out,"res_pseudobulkDESeq2_No.dementia_APOE4_vs_3_astro_DLPFC_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz"))),
                 rbind(fread(fp(out,"res_pseudobulkDESeq2_Dementia_APOE4_vs_3_astro_MTG_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz")),
                       fread(fp(out,"res_pseudobulkDESeq2_No.dementia_APOE4_vs_3_astro_MTG_Cov_sex_braak_atherosclerosis_PMI_ncells_nUmis_LibInput_PCRcycles_avg.mt.csv.gz"))),
                 by=c('gene','cell_type','Cognitive.Status'),suffixes = c(".DLPFC", ".MTG"))


ggplot(res_merge)+geom_point(aes(x=stat.DLPFC,y=stat.MTG),size=0.2) +theme_bw() +facet_wrap('Cognitive.Status')

res_merge[,summary(lm(stat.DLPFC~stat.MTG))$adj.r.squared,by='Cognitive.Status']
# Cognitive.Status        V1
# 1:         Dementia 0.3609409
# 2:      No dementia 0.2284995

#bigger R2 when cognitive status separate, and specifically in Dementia people
#regarding than in MTG only, dementia have a bif effect on transcriptome, 
#we expect than astro response to dementia is visible in MTG because affected in a later stage of the D., 
#so MTG brain environment look more DLPFC as this stage
