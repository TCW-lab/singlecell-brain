
out<-'outputs/04-ROSMAP_MIT_ATAC'
dir.create(out)

source('../../utils/r_utils.R')
library(Seurat)
library(Signac)


#produce the peak count matrix from fragments.tsv file

#1) call peaks by donors, then combine peaks using GenomicRanges::reduce function, and create the read count matrix from these combined called
fragment_files<-list.files('/projectnb/tcwlab/ShareSpace/MIT_ROSMAP_Multiomics/snATAC-seq/',pattern = 'fragments\\.tsv\\.gz$',full.names = TRUE)
donor_ids<-str_remove(basename(fragment_files),'.fragments.tsv.gz')
output_peaks_name<-'initial_peak_calling_by_indivs'

peaks_list<-lapply(fragment_files,function(file){
  donor_id<-str_remove(basename(file),'.fragments.tsv.gz')
  message('Call Peaks for', donor_id)
  gr<-CallPeaks(file,macs2.path = '/projectnb/tcwlab/LabMember/adpelle1/micromamba/envs/macs2/bin/macs2',
                name=donor_id)
  gr$donor_id<-donor_id
  return(gr)
})


#2) combine peaks and reduce, maintaining ident information
peaks.combined <- Reduce(f = c, x = peaks_list)
peaks <- reduce(x = peaks.combined, with.revmap = TRUE)

peaks$peak_called_in <-sapply(1:length(peaks), function(i){
  datasets <-  peaks.combined$donor_id[peaks$revmap[[i]]]
  return(paste(unique(x = datasets), collapse = ","))
})
peaks
peaks$revmap <- NULL


saveRDS(peaks,fp(out,ps(output_peaks_name,"_GenomeRange.rds")))
peaks_dt<-data.table(as.data.frame(peaks))
peaks_dt[,peak_id:=paste(seqnames,start,end,sep="-")]
peaks_dt[,chr:=seqnames]
fwrite(peaks_dt,fp(out,ps(output_peaks_name,".csv.gz")))

#get also the tidy format of these peaks call (keep information of from which donors, each peak have been called)

peaks_dt_inds<-Reduce(rbind,lapply(donor_ids, function(d){
  peaks_dt_sub<-peaks_dt[sapply(strsplit(peak_called_in,','), function(ids)d%in%ids)]
  
  return(peaks_dt_sub[,donor_id:=d][,-'peak_called_in'])
  
}))


peaks_dt_inds[,n.donors:=.N,by="peak_id"]
table(peaks_dt_inds[n.donors==1]$donor_id)
# D19-12513 D19-12514 
#      6543     42380 

fwrite(peaks_dt_inds,fp(out,ps(output_peaks_name,"_tidy.csv.gz")))


#3)  count the fragments falling within each peak in each cell (produce the read-count matrix)
#need the tabix index for each fragment file
# #for 1 
# 
# system('tabix -p bed D19-12514.fragments.tsv.gz')
# fragments_list <- lapply(fragment_files[1],function(file)CreateFragmentObject(file,)) 
# 
# peaks_mat<-FeatureMatrix(fragments_list,
#                          features =peaks,
#                          process_n = 2000)
# dim(peaks_mat)#89647 202463
# 
# as.matrix(peaks_mat[1:10,1:10])
# #cell Calling
# 
# reads_counts<-data.table(cell_id=colnames(peaks_mat),count=colSums(peaks_mat))
# reads_counts[,cell_rank:=rank(-count)]
# ggplot(reads_counts)+
#   geom_line(aes(x=cell_rank,y=count))+scale_x_log10()+scale_y_log10()+
#   geom_hline(yintercept = 500)
# 
# cells<-reads_counts[count>500]$cell_id #5483
# 
# # save the raw and filtered matrix in rds 
# saveRDS(peaks_mat,fp(out1,donor_ids,))
# 
out1<-fp(out,'peak_count_matrices')
dir.create(out1)

message('Peak cell matrix creation')
for(file in fragment_files){
  donor_id<-str_remove(basename(file),'.fragments.tsv.gz')
  message('for ',donor_id)
  
  out2<-fp(out1,donor_id)
  #coutnt fragment falling in peaks
  fragments <- CreateFragmentObject(file)
  
  peaks_mat<-FeatureMatrix(fragments,
                           features =peaks,
                           process_n = 2000)
  dim(peaks_mat)#89647 202463
  
  as.matrix(peaks_mat[1:10,1:10])
  
  # Call Cell if more than 500 reads in peaks
  
  reads_counts<-data.table(cell_id=colnames(peaks_mat),count=colSums(peaks_mat))
  reads_counts[,cell_rank:=rank(-count)]
  ggplot(reads_counts)+
    geom_line(aes(x=cell_rank,y=count))+scale_x_log10()+scale_y_log10()+
    geom_hline(yintercept = 500)
  ggsave(fp(out2,'knee_plot_read_falling_peaks.png'),width = 6,height = 6)
  
  cells<-reads_counts[count>500]$cell_id #5483
  message(length(cells),' cells called')
  # save the raw and filtered matrix in rds 
  saveRDS(peaks_mat,fp(out2,'raw_peak_count_matrix.rds'))
  
  peaks_mat<-peaks_mat[,cells]
  saveRDS(peaks_mat,fp(out2,'raw_peak_count_matrix.rds'))
  
  
  #4) create signac/seurat object with associated metadata

  chrom_assay <- CreateChromatinAssay(
    counts = peaks_mat,
    sep = c("-", "-"),
    fragments = file,
    min.cells = 10,
    min.features = 200
  )
  
  brain <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
    # meta.data = metadata
  )
  brain$donor_id<-donor_id
  
  #save the final object
  saveRDS(brain,fp(out2,"signac_object.rds"))
  
  
}

