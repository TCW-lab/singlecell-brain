source('../../utils/r_utils.R')
out<-'outputs/01-SEAAD_data/'
dir.create(out)
#QC SEA-AD####
#DLPFC
#0) annotate main celltype
mtd<-fread('outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz')
mtd
mtd[,cell_type:=ifelse(str_detect(Class,'Glut'),paste0('Exc_',Subclass),
                       ifelse(str_detect(Class,'GABA'),paste0('Inh_',Subclass),Subclass))]
unique(mtd$cell_type)
mtd[,main_cell_type:=str_extract(cell_type,'Oligo|Exc|Inh|Astro|Mic|Endo|VLMC|OPC')]
mtd[,main_cell_type:=factor(main_cell_type,levels = c('Exc','Inh','Oligo','Astro','OPC','Mic','VLMC','Endo'))]
mtd[,PMI:=as.numeric(PMI)]

ggplot(mtd)+
  geom_bar(aes(x=main_cell_type,fill=main_cell_type),position='dodge')+theme_bw()

#1) Flag outliers donors

#convert column name to minuscule and without space
var_names<-colnames(mtd)
setnames(mtd,var_names,make.names(var_names))

#Clinical Traits####
mts<-unique(mtd,by=c('Donor.ID'))
#Age, PMI, ethnicity
#flag missing info
table(mts$Age.at.Death)
table(mts$Cognitive.Status)
ggplot(mts)+geom_bar(aes(x=Cognitive.Status,fill=Cognitive.Status))+theme_bw( )+scale_fill_manual(values=c(scales::hue_pal()(2),'grey'))
#missing cognitive status
mts[,missing.cognitive.status:=Cognitive.Status=='Reference']


#age
mts[!(missing.cognitive.status),age_decade:=paste0(str_extract(Age.at.Death,'^[1-9]'),'0')]
mts[!(missing.cognitive.status),age_at_death_num:=as.numeric(ifelse(Age.at.Death=='90+',91,Age.at.Death))]

table(mts[!(missing.cognitive.status)]$age_decade,mts[!(missing.cognitive.status)]$Cognitive.Status)
  #    Dementia No dementia
  # 60        2           0
  # 70        3           3
  # 80       12          15
  # 90       22          23

ggplot(mts[!(missing.cognitive.status)])+geom_boxplot(aes(x=Cognitive.Status,y=age_at_death_num,fill=Cognitive.Status))+theme_bw()
ggplot(mts[!(missing.cognitive.status)])+geom_boxplot(aes(y=age_at_death_num)) #ok

mts[,outlier.clinical.status:=missing.cognitive.status]


#PMI
ggplot(mts[!(missing.cognitive.status)])+geom_boxplot(aes(x=Cognitive.Status,y=PMI,fill=Cognitive.Status))+theme_bw() #ok

mts[,outlier.clinical.status:=missing.cognitive.status]

#ethnicity
table(mts[!(missing.cognitive.status)]$Race..choice.White.)
  # Checked Unchecked 
  #      77         3 
races_cols<-colnames(mts)[str_detect(colnames(mts),'^Race')]
races<-str_remove(str_remove(races_cols,'Race..choice.'),'\\.$')
mts[!(missing.cognitive.status),ethnicity:=paste(races[sapply(.SD,function(x)x=='Checked')],collapse = '_'),by='Donor.ID',.SDcols=races_cols]

table(mts$ethnicity)
# Asian                                      White White_American.Indian..Alaska.Native_Other                                White_Other 
# 3                                         74                                          1                                          2 

ggplot(mts[!(missing.cognitive.status)])+
  geom_bar(aes(x=ethnicity,fill=Cognitive.Status))+theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(str_replace_all(x,'\\.',' '),'_','/'), width = 10))

mts[,outlier.ethnicity:=ethnicity=='Asian']

mts[,outlier.clinical.status:=missing.cognitive.status|outlier.ethnicity]

mts[(outlier.clinical.status)] #6

#add this information to the main metadata
mtd<-merge(mtd,mts[,.SD,.SDcols=c('Donor.ID',setdiff(colnames(mts),colnames(mtd)))])


#Cellular/Molecular Traits QC####
#based on %MT
ggplot(mtd[!(outlier.clinical.status)])+geom_violin(aes(x=Donor.ID,y=Fraction.mitochondrial.UMIs))

mtd[,avg.pct.mt:=mean(`Fraction.mitochondrial.UMIs`),by=c('Donor.ID')]
mtd[,avg.pct.mt.ct:=mean(`Fraction.mitochondrial.UMIs`),by=c('Donor.ID','main_cell_type')]

mtsc<-unique(mtd,by=c('Donor.ID','main_cell_type'))

ggplot(mtsc[!(outlier.clinical.status)])+
  geom_boxplot(aes(x=main_cell_type,y=avg.pct.mt.ct,fill=Cognitive.Status))+theme_bw()
#no removal because associated to Dementia


#based on Cellular Distribution
ggplot(mtd[!(outlier.clinical.status)])+
  geom_bar(aes(x=main_cell_type,fill=main_cell_type),position='dodge')+
  facet_wrap(~`Donor.ID`)+theme_bw()

mtd[,n.cells:=.N,by=.(`Donor.ID`)]
mtd[,pct.ct:=.N/n.cells,by=.(`Donor.ID`,main_cell_type)]
mtd[,n.ct:=.N,by=.(`Donor.ID`,main_cell_type)]
mtd[,n.ct.all:=.N,by=.(main_cell_type)]
mtsc<-unique(mtd,by=c('Donor.ID','main_cell_type'))

ggplot(mtsc[!(outlier.clinical.status)])+
  geom_boxplot(aes(x=main_cell_type,y=pct.ct,fill=`Cognitive.Status`),position='dodge')+
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

iqr_skew_mat<-dcast(mtsc,Donor.ID~main_cell_type,value.var = 'pct.from.IQR')
iqr_skew_pca<-RunPca(t(data.frame(iqr_skew_mat,row.names = 'Donor.ID')),scale = T)
iqr_skew_pca$x

mts<-merge(mts,unique(mtsc[,.(Donor.ID,avg.pct.mt)]))

iqr_skew_pca_dt<-PcaPlot(iqr_skew_pca,mts,group.by ='Cognitive.Status',sample_col = 'Donor.ID',return_pcs_mtd = T)
iqr_skew_pca_dt[PC1<(-10)]
ggplot(iqr_skew_pca_dt)+geom_boxplot(aes(x=Cognitive.Status,y=PC1))

#outlier if sample > 3*IQR of PC1
iqr_skew_pca_dt[,outlier.cellprop:=PC1%in%boxplot.stats(PC1,coef=3)$out]
iqr_skew_pca_dt[(outlier.cellprop)]

PcaPlot(iqr_skew_pca,mts,group.by ='Cognitive.Status',
        sample_col = 'Donor.ID',label = iqr_skew_pca_dt$outlier.cellprop)


# iqr_skew_pca_dt[,outlier.cellprop.pc2:=PC2%in%boxplot.stats(PC2,coef=3)$out]
# iqr_skew_pca_dt[(outlier.cellprop.pc2)]
# PcaPlot(iqr_skew_pca,mts,group.by ='Cognitive.Status',
#         sample_col = 'Donor.ID',label = iqr_skew_pca_dt$outlier.cellprop.pc2)

# PcaPlot(iqr_skew_pca,mts,group.by ='avg.pct.mt',
#                             sample_col = 'Donor.ID',label = iqr_skew_pca_dt$outlier.cellprop.pc2,return_pcs_mtd = T)

ggplot(iqr_skew_pca_dt)+geom_point(aes(x=PC1,y=avg.pct.mt,col=outlier.cellprop))+theme_bw()

iqr_skew_pca_dt[,cell_prop_dev_pc1:=PC1]

iqr_skew_pca_dt[,outlier.cellular.status:=outlier.cellprop]
iqr_skew_pca_dt[(outlier.cellular.status)]$Donor.ID

mtsc<-merge(mtsc,iqr_skew_pca_dt[,.(Donor.ID,outlier.cellprop,cell_prop_dev_pc1,outlier.cellular.status)])


mtd<-merge(mtd,iqr_skew_pca_dt[,.(Donor.ID,outlier.cellprop,cell_prop_dev_pc1,outlier.cellular.status)])


mtd[,outlier:=outlier.cellular.status|outlier.clinical.status]



mts<-unique(mtd,by=c('Donor.ID'))
tot.outliers<-mts[(outlier)]$Donor.ID
length(tot.outliers) #12/84

#save 
fwrite(mtd,'outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_metadata.csv.gz')

mts<-RemoveUselessColumns(mts,key_cols=c('Donor.ID'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.|Class|Subclass|Supertype|cell_type')

fwrite(mts,'outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_sample_level_metadata.csv.gz')


mtsc<-RemoveUselessColumns(mtsc,key_cols=c('Donor.ID','main_cell_type'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.|Subclass|Supertype')

fwrite(mtsc,'outputs/01-SEAAD_data/DLPFC/all_final_RNAseq_nuclei_main_cell_type_level_metadata.csv.gz')

#MTG####
#1) Flag outliers donors
mtd<-fread('outputs/01-SEAAD_data/MTG/all_final_RNAseq_nuclei_metadata.csv.gz')
mtd
mtd[,cell_type:=ifelse(str_detect(Class,'Glut'),paste0('Exc_',Subclass),
                       ifelse(str_detect(Class,'GABA'),paste0('Inh_',Subclass),Subclass))]
unique(mtd$cell_type)
mtd[,main_cell_type:=str_extract(cell_type,'Oligo|Exc|Inh|Astro|Mic|Endo|VLMC|OPC')]
mtd[,main_cell_type:=factor(main_cell_type,levels = c('Exc','Inh','Oligo','Astro','OPC','Mic','VLMC','Endo'))]
mtd[,PMI:=as.numeric(PMI)]


#convert column name to minuscule and without space
var_names<-colnames(mtd)
setnames(mtd,var_names,make.names(var_names))

#Clinical Traits
mts<-unique(mtd,by=c('Donor.ID'))
#Age, PMI, ethnicity
#flag missing info
table(mts$Cognitive.Status)
ggplot(mts)+geom_bar(aes(x=Cognitive.Status,fill=Cognitive.Status))+theme_bw( )+scale_fill_manual(values=c(scales::hue_pal()(2),'grey'))
#missing cognitive status
mts[,missing.cognitive.status:=Cognitive.Status=='Reference']


#age
mts[!(missing.cognitive.status),age_decade:=factor(paste0(str_extract(Age.at.Death,'^[1-9]'),'0'))]
mts[!(missing.cognitive.status),age_at_death_num:=as.numeric(ifelse(Age.at.Death=='90+',91,Age.at.Death))]

table(mts[!(missing.cognitive.status)]$age_decade,mts[!(missing.cognitive.status)]$Cognitive.Status)
# Dementia No dementia
# 60        3           0
# 70        3           3
# 80       13          16
# 90       23          23

ggplot(mts[!(missing.cognitive.status)])+geom_boxplot(aes(x=Cognitive.Status,y=age_at_death_num,fill=Cognitive.Status))+theme_bw()
ggplot(mts[!(missing.cognitive.status)])+geom_boxplot(aes(y=age_at_death_num)) #ok



#PMI
ggplot(mts[!(missing.cognitive.status)])+geom_boxplot(aes(x=Cognitive.Status,y=PMI,fill=Cognitive.Status))+theme_bw() #ok


#ethnicity
table(mts[!(missing.cognitive.status)]$Race..choice.White.)
# Checked Unchecked 
#      77         3 
races_cols<-colnames(mts)[str_detect(colnames(mts),'^Race')]
races<-str_remove(str_remove(races_cols,'Race..choice.'),'\\.$')
mts[!(missing.cognitive.status),ethnicity:=paste(races[sapply(.SD,function(x)x=='Checked')],collapse = '_'),by='Donor.ID',.SDcols=races_cols]
table(mts$ethnicity)

# Asian                                      White 
# 3                                         78 
# White_American.Indian..Alaska.Native_Other                                White_Other 
# 1                                          2 
  
ggplot(mts[!(missing.cognitive.status)])+geom_bar(aes(x=ethnicity,fill=Cognitive.Status))+theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(str_replace_all(x,'\\.',' '),'_','/'), width = 10))

mts[,outlier.ethnicity:=ethnicity=='Asian']

mts[,outlier.clinical.status:=missing.cognitive.status|outlier.ethnicity]

mts[(outlier.clinical.status)] #8

mtd<-merge(mtd,mts[,.SD,.SDcols=c('Donor.ID',setdiff(colnames(mts),colnames(mtd)))])
# table(mts[!(outlier.clinical.status)]$Cognitive.Status)
# Dementia No dementia 
# 40          41 

#Cellular/Molecular Traits
#based on %MT
ggplot(mtd[!(outlier.clinical.status)])+geom_violin(aes(x=Donor.ID,y=Fraction.mitochondrial.UMIs))

mtd[,avg.pct.mt:=mean(`Fraction.mitochondrial.UMIs`),by=c('Donor.ID')]
mtd[,avg.pct.mt.ct:=mean(`Fraction.mitochondrial.UMIs`),by=c('Donor.ID','main_cell_type')]

mtsc<-unique(mtd,by=c('Donor.ID','main_cell_type'))

ggplot(mtsc[!(outlier.clinical.status)])+
  geom_boxplot(aes(x=main_cell_type,y=avg.pct.mt.ct,fill=Cognitive.Status))+theme_bw()

#no removal because associated to Dementia

#based on Cellular Distribution
ggplot(mtd[!(outlier.clinical.status)])+
  geom_bar(aes(x=main_cell_type,fill=main_cell_type),position='dodge')+
  facet_wrap(~`Donor.ID`)+theme_bw()

mtd[,n.cells:=.N,by=.(`Donor.ID`)]
mtd[,pct.ct:=.N/n.cells,by=.(`Donor.ID`,main_cell_type)]
mtd[,n.ct:=.N,by=.(`Donor.ID`,main_cell_type)]
mtd[,n.ct.all:=.N,by=.(main_cell_type)]
mtsc<-unique(mtd,by=c('Donor.ID','main_cell_type'))

ggplot(mtsc[!(outlier.clinical.status)])+
  geom_boxplot(aes(x=main_cell_type,y=pct.ct,fill=`Cognitive.Status`),position='dodge')+
  theme_bw()


#PCA of the % from IQR
source('../../utils/pca_utils.R')

mtsc[,pct.ct.q25:=quantile(pct.ct,0.25),by='main_cell_type']
mtsc[,pct.ct.q75:=quantile(pct.ct,0.75),by='main_cell_type']

mtsc[,pct.ct.IQR:=pct.ct.q75-pct.ct.q25]
mtsc[pct.ct<pct.ct.q25,pct.from.IQR:=(pct.ct.q25-pct.ct)/pct.ct.IQR]
mtsc[pct.ct>pct.ct.q75,pct.from.IQR:=(pct.ct-pct.ct.q75)/pct.ct.IQR]
mtsc[is.na(pct.from.IQR),pct.from.IQR:=0]
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.003323  0.327975  0.362668 13.114971 

iqr_skew_mat<-dcast(mtsc,Donor.ID~main_cell_type,value.var = 'pct.from.IQR')
iqr_skew_pca<-RunPca(t(data.frame(iqr_skew_mat,row.names = 'Donor.ID')),scale = T)
iqr_skew_pca$x

mts<-merge(mts,unique(mtsc[,.(Donor.ID,avg.pct.mt)]))

iqr_skew_pca_dt<-PcaPlot(iqr_skew_pca,mts,group.by ='Cognitive.Status',sample_col = 'Donor.ID',return_pcs_mtd = T)
iqr_skew_pca_dt[PC1<(-10)]
ggplot(iqr_skew_pca_dt)+geom_boxplot(aes(x=Cognitive.Status,y=PC1))

iqr_skew_pca_dt[,outlier.cellprop:=PC1%in%boxplot.stats(PC1,coef=3)$out]
iqr_skew_pca_dt[(outlier.cellprop)]

PcaPlot(iqr_skew_pca,mts,group.by ='Cognitive.Status',
        sample_col = 'Donor.ID',label = iqr_skew_pca_dt$outlier.cellprop)


# iqr_skew_pca_dt[,outlier.cellprop.pc2:=PC2%in%boxplot.stats(PC2,coef=3)$out]
# iqr_skew_pca_dt[(outlier.cellprop.pc2)]
# PcaPlot(iqr_skew_pca,mts,group.by ='Cognitive.Status',
#         sample_col = 'Donor.ID',label = iqr_skew_pca_dt$outlier.cellprop.pc2)

# PcaPlot(iqr_skew_pca,mts,group.by ='avg.pct.mt',
#                             sample_col = 'Donor.ID',label = iqr_skew_pca_dt$outlier.cellprop.pc2,return_pcs_mtd = T)

ggplot(iqr_skew_pca_dt)+geom_point(aes(x=PC1,y=avg.pct.mt,col=outlier.cellprop))+theme_bw()



iqr_skew_pca_dt[,cell_prop_dev_pc1:=PC1]

iqr_skew_pca_dt[,outlier.cellular.status:=outlier.cellprop]
iqr_skew_pca_dt[(outlier.cellular.status)]$Donor.ID

mtsc<-merge(mtsc,iqr_skew_pca_dt[,.(Donor.ID,outlier.cellprop,cell_prop_dev_pc1,outlier.cellular.status)])


mtd<-merge(mtd,iqr_skew_pca_dt[,.(Donor.ID,outlier.cellprop,cell_prop_dev_pc1,outlier.cellular.status)])


mtd[,outlier:=outlier.cellular.status|outlier.clinical.status]



mts<-unique(mtd,by=c('Donor.ID'))
tot.outliers<-mts[(outlier)]$Donor.ID
length(tot.outliers) #12/84

#save 
fwrite(mtd,'outputs/01-SEAAD_data/MTG/all_final_RNAseq_nuclei_metadata.csv.gz')

mts<-RemoveUselessColumns(mts,key_cols=c('Donor.ID'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.|Class|Subclass|Supertype|cell_type')
fwrite(mts,'outputs/01-SEAAD_data/MTG/all_final_RNAseq_nuclei_sample_level_metadata.csv.gz')


mtsc<-RemoveUselessColumns(mtsc,key_cols=c('Donor.ID','main_cell_type'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.|Subclass|Supertype')
fwrite(mtsc,'outputs/01-SEAAD_data/MTG/all_final_RNAseq_nuclei_main_cell_type_level_metadata.csv.gz')


#Single Cell Object and Pseudobulk data Creations####
#generate the per cell type object. 
unique(mtd$Supertype)
#run 01Di-
CreateJobForRfile('scripts/01Di-create_seurat_celltype.R',nThreads = 28)
RunQsub('scripts/01Di-create_seurat_celltype.R',job_name = 'SEActSeur')

#create pseudobulk matrix
#by cell_type
#for 1
library(Seurat)
astro<-readRDS('outputs/01-SEAAD_data/DLPFC/Astrocyte.rds')
pseudo_mat<-AggregateExpression(astro,assays = 'RNA',slot = 'count',group.by = 'Donor.ID',
                    return.seurat = FALSE)
fwrite(data.table(pseudo_mat$RNA,keep.rownames = 'gene_id'),'outputs/01-SEAAD_data/DLPFC/Astrocyte_pseudobulk.csv.gz')


#for all
#run 01Dii
CreateJobForRfile('scripts/01Dii-create_pseudobulk_cell_type.R',nThreads = 28)
RunQsub('scripts/01Dii-create_pseudobulk_cell_type.R',wait_for = 'SEActSeur',job_name ='SEAPseudo' )

#merge the sets, create metadata by cell type, aggregate by main cell type, create metdata by main cell type

#1) merge the sets
for(region in c('DLPFC','MTG')){
  print(region)
  out1<-fp(out,region)
  mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))
  
  mtd[,original_cell_type_name:=str_remove(str_replace_all(cell_type,' ','_'),'^Inh_|Exc_')]
  mtd[,original_cell_type_name:=str_replace_all(original_cell_type_name,'/','-')]
  fwrite(mtd,fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))
  

  
  mtsc<-unique(mtd,by=c('Donor.ID','cell_type'))
  
  pseudo_files<-list.files(out1,pattern = '\\_pseudobulk\\.csv\\.gz',full.names = T)
  pseudo_files_dt<-data.table(file=pseudo_files,
                              original_cell_type_name=str_remove(basename(pseudo_files),'_pseudobulk\\.csv\\.gz'))
  pseudo_files_dt[,original_cell_type_name:=str_remove(original_cell_type_name,'\\_set[0-9]+')]
  
  
  pseudo_files_dt<-merge(pseudo_files_dt,unique(mtsc[,.(original_cell_type_name,cell_type,main_cell_type)]))
  pseudo_files_dt[,n_file:=.N,by='cell_type']
  
  for(ct in unique(pseudo_files_dt[n_file>1]$original_cell_type_name)){
    pseudo_list<-lapply(pseudo_files_dt[original_cell_type_name==ct]$file, function(f)fread(f))
    
    #lacking samples column
    samples<-Reduce(union,lapply(pseudo_list,colnames))
    pseudo_list<-lapply(pseudo_list,function(x){
      samples_lacking<-setdiff(samples,colnames(x))
      x[,(samples_lacking):=0]
    })
    
    #transform as matrix
    pseudo_list<-lapply(pseudo_list, function(x)as.matrix(data.frame(x,row.names = 'gene_id')))
    
    #aggregate count by sample / featute
    features<-rownames(pseudo_list[[1]])
    samples<-colnames(pseudo_list[[1]])
    
    pseudo_merge<-Reduce(`+`,lapply(pseudo_list,function(x)x[features,samples]))
    print(head(pseudo_merge[,1:10]))
    fwrite(data.table(pseudo_merge,keep.rownames = 'gene_id'),fp(out1,paste0(ct,'pseudobulk.csv.gz')))
    
  }
  
  #remove sets
  system(paste('rm',paste(pseudo_files_dt[n_file>1]$file,collapse = ' ')))
  
}
  
#2) create mtd for each pseudobulk 
for(region in c('DLPFC','MTG')){
  print(region)
  out1<-fp(out,region)
  mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))
  mtd[,tot.cells.donor:=.N,by=.(`Donor.ID`)]
  mtd[,n.cells:=.N,by=.(`Donor.ID`,cell_type)]
  mtd[,prop.cells:=n.cells/tot.cells.donor,by=.(`Donor.ID`,cell_type)]
  
  mtd[,med.umis.per.cell:=median(Number.of.UMIs,na.rm = T),by=.(`Donor.ID`,cell_type)]
  mtd[,med.genes.per.cell:=median(Genes.detected,na.rm = T),by=.(`Donor.ID`,cell_type)]
  mtd[,avg.pct.mt.per.cell:=mean(`Fraction.mitochondrial.UMIs`),by=c('Donor.ID','cell_type')]
  
  mtsc<-unique(mtd,by=c('Donor.ID','cell_type'))
  
  #flag donors with not enough cells
  mtsc[,pass.threshold.n.cells:=n.cells>50]
  mtsc[,outlier.n.cells:=!pass.threshold.n.cells]
  
  mtscf<-RemoveUselessColumns(mtsc,key_cols=c('Donor.ID','cell_type'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.')
  
  fwrite(mtscf,fp(out1,'all_final_RNAseq_nuclei_cell_type_level_metadata.csv.gz'))
}

#3) merge by main_cell_type
for(region in c('DLPFC','MTG')){
  print(region)
  out1<-fp(out,region)
  mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))
  
  out2<-fp(out1,'pseudobulk_main_cell_type')
  dir.create(out2)
  pseudo_files<-list.files(out1,pattern = '\\_pseudobulk\\.csv\\.gz',full.names = T)
  to_merge<-data.table(file=pseudo_files,
                       original_cell_type_name=str_remove(basename(pseudo_files),'_pseudobulk\\.csv\\.gz'))
  
  
  to_merge<-merge(to_merge,unique(mtsc[,.(original_cell_type_name,cell_type,main_cell_type)]))
  
  for(ct in unique(to_merge$main_cell_type)){
    message(ct)
    pseudo_list<-lapply(to_merge[main_cell_type==ct]$file, function(f)fread(f))
    
    #lacking samples column
    samples<-Reduce(union,lapply(pseudo_list,colnames))
    pseudo_list<-lapply(pseudo_list,function(x){
      samples_lacking<-setdiff(samples,colnames(x))
      x[,(samples_lacking):=0]
    })
    
    #transform as matrix
    pseudo_list<-lapply(pseudo_list, function(x)as.matrix(data.frame(x,row.names = 'gene_id')))
    
    #aggregate count by sample / featute
    features<-rownames(pseudo_list[[1]])
    samples<-colnames(pseudo_list[[1]])
    
    pseudo_merge<-Reduce(`+`,lapply(pseudo_list,function(x)x[features,samples]))
    print(head(pseudo_merge[,1:10]))
    fwrite(data.table(pseudo_merge,keep.rownames = 'gene_id'),fp(out2,paste0(ct,'.csv.gz')))
    
    
  }
}
  
#4) create mtd for each pseudobulk main cell type 
for(region in c('DLPFC','MTG')){
  print(region)
  out1<-fp(out,region)
  mtd<-fread(fp(out1,'all_final_RNAseq_nuclei_metadata.csv.gz'))
  out2<-fp(out1,'pseudobulk_main_cell_type')
  dir.create(out2)
  
  
  mtd[,tot.cells.donor:=.N,by=.(`Donor.ID`)]
  mtd[,n.cells:=.N,by=.(`Donor.ID`,main_cell_type)]
  
  mtd[,prop.cells:=n.cells/tot.cells.donor,by=.(`Donor.ID`,main_cell_type)]
  
  mtd[,med.umis.per.cell:=median(Number.of.UMIs,na.rm = T),by=.(`Donor.ID`,main_cell_type)]
  mtd[,med.genes.per.cell:=median(Genes.detected,na.rm = T),by=.(`Donor.ID`,main_cell_type)]
  
  mtd[,avg.pct.mt.per.cell:=mean(`Fraction.mitochondrial.UMIs`),by=c('Donor.ID','main_cell_type')]
  
  mtsc<-unique(mtd,by=c('Donor.ID','main_cell_type'))
  #flag donors with not enough cells
  mtsc[,pass.threshold.n.cells:=n.cells>50,by=c('Donor.ID','main_cell_type')]
  mtsc[,outlier.n.cells:=!pass.threshold.n.cells,by=c('Donor.ID','main_cell_type')]
  
  mtscf<-RemoveUselessColumns(mtsc,key_cols=c('Donor.ID','main_cell_type'),pattern_to_exclude = 'ATAC|Multiome|Doublet|Number.of|Genes.')
  
  fwrite(mtscf,fp(out2,'all_final_RNAseq_nuclei_main_cell_type_level_metadata.csv.gz'))
  
  for(ct in unique(mtscf$main_cell_type)){
    
    fwrite(mtscf[main_cell_type==ct],fp(out2,paste0(ct,'_metadata.csv.gz')))
    
  }
  
}



#Cell Level QC####
# BiocManager::install("flexmix")
# BiocManager::install("miQC")
library(miQC)
#remotes::install_github('satijalab/seurat-wrappers@community-vignette')
library(Seurat)

library(SeuratWrappers)
library(flexmix)


#nFeature, nRNA : 3*IQR 
#%MT <0.05
#for 1####

library(Seurat)
astro<-readRDS('outputs/01-')
VlnPlot(astro,features =  c('Fraction.mitochondrial.UMIs','Genes.detected','Number.of.UMIs'),combine = F)
ggsave(fp(out2,ps(ct,'qc_cell_metrics.png')),width = 8,height = 6)

boxres<-boxplot.stats(astro$Number.of.UMIs,coef = 3)
astro$UMIs.outlier<-colnames(astro)%in%names(boxres$out)
VlnPlot(astro,features =  'Number.of.UMIs',split.by = 'UMIs.outlier')
ggsave(fp(out2,ps(ct,'qc_nUMIs.png')),width = 4,height = 6)

VlnPlot(astro,features =  'Number.of.UMIs',split.by = 'UMIs.outlier',group.by = 'Donor.ID',fill.by = 'outlier.cellprop')
ggsave(fp(out2,ps(ct,'qc_nUMIs_per_donor.png')),width = 10,height = 6)

boxres<-boxplot.stats(astro$Genes.detected,coef = 3)
astro$Genes.outlier<-colnames(astro)%in%names(boxres$out)
VlnPlot(astro,features =  'Genes.detected',split.by = 'Genes.outlier')
ggsave(fp(out2,ps(ct,'qc_nGenes.png')),width = 4,height = 6)

VlnPlot(astro,features =  'Genes.detected',split.by = 'Genes.outlier',group.by = 'Donor.ID',fill.by = 'outlier.cellprop')
ggsave(fp(out2,ps(ct,'qc_nGenes_per_donor.png')),width = 10,height = 6)

#percent.mt
#try use MIQC
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.
#?filterCells
astro_sce<-as.SingleCellExperiment(astro)
astro_sce@colData$Genes.detected
model <- mixtureModel(astro_sce,
                      subsets_mito_percent = 'Fraction.mitochondrial.UMIs',detected = 'Genes.detected')

model <- flexmix(Fraction.mitochondrial.UMIs ~ Genes.detected, data = astro@meta.data,
                 k = 2)

intercept1 <- parameters(model, component = 1)[1]
intercept2 <- parameters(model, component = 2)[1]
if (intercept1 > intercept2) {
  compromised_dist <- 1
  intact_dist <- 2
}else {
  intact_dist <- 1
  compromised_dist <- 2
}

astro$Mito.outlier <- post[, compromised_dist] > 0.99
VlnPlot(astro,features =  'Fraction.mitochondrial.UMIs',split.by = 'Mito.outlier')
FeatureScatter(astro,'Fraction.mitochondrial.UMIs','Genes.detected',group.by ='Mito.outlier' )
#==> do not use MIQC for scNuc https://github.com/TCW-lab/SingleCell_APOE44/issues/2

boxres<-boxplot.stats(astro$Fraction.mitochondrial.UMIs,coef = 3)


astro$Mito.high<-colnames(astro)%in%names(boxres$out)
VlnPlot(astro,features =  'Fraction.mitochondrial.UMIs',split.by = 'Mito.high')
ggsave(fp(out2,ps(ct,'qc_Mito.png')),width = 4,height = 6)

VlnPlot(astro,features =  'Fraction.mitochondrial.UMIs',split.by = 'Mito.high',group.by = 'Donor.ID',fill.by = 'outlier.cellprop')
ggsave(fp(out2,ps(ct,'qc_Mito_per_donor.png')),width = 10,height = 6)

#For scNuc, keep hardthreshold of 5% mito
astro$Mito.outlier<-astro$Fraction.mitochondrial.UMIs>0.05

#Final outlier is:
astro$donor.outlier<-astro$outlier #conserved donor outlier metadata compute in previous step
astro$cell.outlier<-astro$Genes.outlier|astro$UMIs.outlier|astro$Mito.outlier

astro$outlier<-astro$donor.outlier|astro$cell.outlier

message(round(sum(astro$donor.outlier)/nrow(astro),digits = 1),'% cells flagged as donors outliers')
message(round(sum(astro$cell.outlier)/nrow(astro),digits = 1),'% cells flagged as cells outliers' )

message('In total ',round(sum(astro$outlier)/nrow(astro),digits = 1),'% cells flagged as outliers')

#save the full object
saveRDS(astro,file)

#create the QCed good quqlity ~10Kcell reduced object 
celltype_qc<-astro[,!astro$outlier]
celltype_qc 

#get the top10% donors
good_qual_donors<-unique(astro$Donor.ID[!astro$outlier&astro$avg.pct.mt<=quantile(astro$avg.pct.mt,0.10)&astro$cell_prop_dev_pc1<=quantile(astro$cell_prop_dev_pc1,0.90)])
length(good_qual_donors)#12
celltype_qc_small<-subset(astro,Donor.ID%in%good_qual_donors)
saveRDS(celltype_qc_small,fp(out2,ps(ct,'_qc_small.rds')))


#for all####
#run 01Diii
CreateJobForRfile('scripts/01Diii-cell_level_QC.R',nThreads = 28)
RunQsub('scripts/01Diii-cell_level_QC.R',job_name ='SEACellQC' )


#+ Create HighQual small rds object by cell type
#run 01Div
CreateJobForRfile('scripts/01Div-create_small_QCed_object.R',nThreads = 28)
RunQsub('scripts/01Div-create_small_QCed_object.R',job_name ='SEACellSmall' )

#Next TODO: how to use/presentation of this dataset
#+ how to perform pseudobulk Analysis
#+ in a notebook
#+ 


