#MultiBrainRegion
out<-'outputs/07-ROSMAP_MIT_MultiRegion'
dir.create(out)
source('../../utils/r_utils.R')
library(Seurat)


angular=readRDS('~/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/MultiRegion/Angular_gyrus.rds')

length(table(angular$projid))
