
source('../../utils/r_utils.R')
out<-'outputs/04-ROSMAP_MIT_ATAC/'
fread('../../../../ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/D19-12513.fragments.tsv.gz')
fragment_files<-list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/',pattern = 'fragments\\.tsv\\.gz$',full.names = TRUE)

n_fragments<-rbindlist(lapply(fragment_files, function(file){
  libID=str_remove(basename(file),'.fragments.tsv.gz')
  fragments<-fread(file,select = 4,col.names = 'cell_id')
  return(fragments[,.N,by='cell_id'][,libraryID:=libID][N>500])
}))
fwrite(n_fragments,fp(out,'n_fragment_per_barcode_alllibs.csv.gz'))

table(n_fragments$libraryID)

