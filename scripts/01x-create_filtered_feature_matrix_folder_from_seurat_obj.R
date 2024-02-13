
#save matrix and metadata for exercise

library(Matrix)
out1<-fp(out,'Young_coding_exercise')
dir.create(out1)
# save sparse matrix
out2<-fp(out1,'SEA_AD34k_filtered_feature_bc_matrix')
dir.create(out2)

writeMM(obj = brain20kf@assays$RNA@counts, file=fp(out2,"matrix.mtx"))
#system: gzip matrix.mtx 

# save genes and cells names
fwrite(list(rownames(brain20kf@assays$RNA@counts)), file = fp(out2,"features.tsv.gz"))
fwrite(x = list(colnames(brain20kf@assays$RNA@counts)), file = fp(out2,"barcodes.tsv.gz"))


#metadata without cell type info
head(brain20kf[[]])
fwrite(data.table(brain20kf@meta.data,keep.rownames = 'cell_id')[,.(cell_id,donor_id,tissue,disease,sex,self_reported_ethnicity,`Age at death`,`Years of education`,`APOE4 status`)],
       fp(out1,'SEA_AD34k_metadata.csv.gz'))

#test
mat<-Read10X(out2,gene.column = 1)
dim(mat)
head(mat[,1:10])
colnames(mat)
mtd<-fread(fp(out1,'SEA_AD34k_metadata.csv.gz'))
mtd     
#OK
