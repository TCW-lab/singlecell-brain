#create the QCed good quqlity ~10Kcell reduced object 

source('../../utils/r_utils.R')
library(Seurat)
out<-'outputs/01-SEAAD_data/'

# library(SeuratWrappers)
# library(flexmix)


for(region in c('DLPFC','MTG')){
  message('Within ', region,' data...')
  out1<-fp(out,region)
  out2<-fp(out1,'QC_small')
  dir.create(out2)
  rds_files<-list.files(out1,pattern = '.rds$',full.names = TRUE)
  
  for(file in rds_files){
    ct<-str_remove(basename(file),pattern = '.rds$')
    message('cells QC of', ct)
    celltype<-readRDS(file)
    
    celltype$outlier<-celltype$donor.outlier|celltype$cell.outlier
    saveRDS(celltype,file)
    #create the QCed good quqlity ~10Kcell reduced object 
    celltype_qc<-celltype[,!celltype$outlier]
    celltype_qc 
    
    #get the top10% donors
    good_qual_donors<-unique(celltype_qc$Donor.ID[!celltype_qc$outlier&celltype_qc$avg.pct.mt<=quantile(celltype_qc$avg.pct.mt,0.10)&celltype_qc$cell_prop_dev_pc1<=quantile(celltype_qc$cell_prop_dev_pc1,0.90)])
    length(good_qual_donors)#12
    celltype_qc_small<-subset(celltype_qc,Donor.ID%in%good_qual_donors)
    saveRDS(celltype_qc_small,fp(out2,ps(ct,'_qc_small.rds')))
  }
  
}
