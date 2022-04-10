
# prepare input files for matrix-eqtl

# -- load libraries
library(dplyr)

# -- load command line args

# -- load PDUI matrix (imputed and normalied)
pdui_mat <- read.table("./Matrix_eQTL/pdui_mat.imputed_qnorm.txt",header=T,sep="\t")

# -- load covariates
covariate_mat <- read.table("./Matrix_eQTL/pdui.peer.covariates.txt",header=T,sep="\t",stringsAsFactors=F)

# -- reorder columns of PDUI matrix
id_order <- colnames(pdui_mat)[-1]
covariate_mat.reorder <- covariate_mat %>% dplyr::select("id",all_of(id_order))
rm(covariate_mat)
# -- write covariates and phenotype matrix to file
write.table(pdui_mat,file="./Matrix_eQTL/Phenotype_matrix.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(covariate_mat.reorder,file="./Matrix_eQTL/Covariate_matrix.txt",quote=F,sep="\t",row.names=F,col.names=T)


# -- load genotype matrix
gt_mat <- read.table("./Matrix_eQTL/genotype_matrix.bed",header=T,sep="\t")

gt_mat.reorder <- gt_mat %>% dplyr::select("id",all_of(id_order))
rm(gt_mat,id_order)
write.table(gt_mat.reorder,file="./Matrix_eQTL/Genotype_matrix.txt",quote=F,sep="\t",row.names=F,col.names=T)

print("Done!")
