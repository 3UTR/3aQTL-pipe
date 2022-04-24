#!/opt/app/languages/R-3.6.3/bin/Rscript
# 2022-01-08
library(optparse)

# -- global variable
option_list <- list(
	make_option(c("-p", "--pheno_data"),type = "character", default = "Dapars2_res.all_chromosomes.txt", action = "store", help = "Proivde the merged ouput from DaPars2, Dapars2_res.all_chromosomes.txt in default"),
	make_option(c("-g","--geno_pca"),type = "character", default = "./Matrix_eQTL/genotype_pca.eigenvec", action = "store", help = "the eigenvector of genotpye pca analysis, ./Matrix_eQTL/genotype_pca.eigenvec in default"),
	make_option(c("-c","--known_covs"),type = "character", default = "NA", action = "store", help = "input a text file contains known covariates if available, NA in defaut"),
	make_option(c("-n","--top_N_pca"),type = "integer", default = "5", action = "store", help = "specify the top N PCA on gt would be used, defaut value = 5")
		    )
opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))
apa_res_file <- opt$pheno_data
gtPCA_file <- opt$geno_pca
known_cov_file <- opt$known_covs
topN_pca <- opt$top_N_pca

cat('Arguments:','\n',
    '--pheno_data',apa_res_file,'\n',
    '--geno_pca',gtPCA_file,'\n',
    '--known_covs',known_cov_file,'\n',
    '--top_N_pca',topN_pca,'\n')

cat('Current directory:')
getwd()

library(dplyr)
library(peer)
library(impute)
# --------------- prepare covariates -----------------
# load genotype pca
gt_pca <- read.table(gtPCA_file,header=F,sep=" ",stringsAsFactors=F)
N <- as.integer(topN_pca) + 1
gt_pca$V1 <- NULL;gt_pca <- gt_pca[,1:N]
names(gt_pca) <- c("subject_id",paste0("PC_",1:(N-1)))
rm(N)
# add known covariates into topN gt_pca if available
if(known_cov_file!="NA"){
	known_cov <- read.table(known_cov_file,header=T,sep="\t",stringsAsFactors=F)
	dim(known_cov)
	N <- dim(known_cov)[2]
	for(i in 2:N){
		if(class(known_cov[,i])=="character"){
			known_cov[,i] <- as.factor(known_cov[,i])
			known_cov[,i] <- as.numeric(known_cov[,i])
		}
	}
	col_names <- names(known_cov)[2:N]
	names(known_cov) <- c("subject_id",col_names)
	gt_pca <- merge(gt_pca,known_cov,by="subject_id")
}


# convert covariates' data.frame to matrix 
rownames(gt_pca) <- gt_pca$subject_id;gt_pca <- as.matrix(gt_pca[,-1])
cat('Dimension of gt_pca:',dim(gt_pca),'\n')
rm(known_cov,known_cov_file)

cat("Start phenotype matrix\n","Open APA results file:",apa_res_file,"\n")
# --------------- prepare phenotype matrix ------------------
pdui_mat <- read.table(apa_res_file, stringsAsFactors=FALSE, header=TRUE,sep="\t",check.names=FALSE)
pdui_mat <- pdui_mat[,-c(2,3,4)]

pdui_mat.sel <- pdui_mat %>% dplyr::select(all_of(rownames(gt_pca)))
pdui_mat.sel <- as.matrix(pdui_mat.sel)
rownames(pdui_mat.sel) <- pdui_mat[,1]


#remove genes with more than 50% entries missing and individuals with more than 80% missing data
pdui_mat.sel <- pdui_mat.sel[, colMeans(is.na(pdui_mat.sel)) <= 0.8];pdui_mat.sel <-  pdui_mat.sel[rowMeans(is.na(pdui_mat.sel)) < 0.5,]
class(pdui_mat.sel) <- 'numeric'


# run peer to estimate confounders
cat("Start covariate analysis by peer...")
#save.image(file="run_peer_impute.RData")
model <- PEER()
covs_se <- gt_pca

PEER_setCovariates(model, covs_se)
dim(PEER_getCovariates(model))
#impute missing values in PDUI matrix
mat.ds <- pdui_mat.sel
mat_impute <- impute.knn(mat.ds)
#quantile normalization
df_w <- as.data.frame(mat_impute$data)
for(gene in 1:nrow(df_w)){
	mat = df_w[gene,]
	mat = apply(mat,1,rank,ties.method = "average")
	mat = qnorm(mat / (ncol(df_w)+1))
	df_w[gene,] = mat
}

pdui_mat <- cbind(rownames(mat_impute$data),df_w)
y <- colnames(pdui_mat)[-1]
colnames(pdui_mat) <- c("Gene",y)
id_order <- colnames(pdui_mat)[-1]
cat("Output phenotype matrix:\n")
write.table(pdui_mat,file="./Matrix_eQTL/Phenotype_matrix.txt",row.names=F,col.names=T,quote=F,sep="\t")
rm(y)

PEER_setPhenoMean(model, t(as.matrix(mat_impute$data)))

dim(PEER_getPhenoMean(model))

# set number of peer factors
## N < 150, use 15  PEERs, 150<=N<250, use 30 PEERs, N >=250 use 35 PEERs
if (ncol(mat.ds) < 150) {
	numcov <- 15
} else if (ncol(mat.ds) < 250) {
	numcov <- 30
} else if (ncol(mat.ds) >= 250) {
	numcov <- 35
}

PEER_setNk(model, numcov)
PEER_getNk(model)

PEER_update(model)

# diag
pdf('peer.diag.pdf', width=6, height=8)
PEER_plotModel(model)
dev.off()


factors = t(PEER_getX(model))
weights = PEER_getW(model)
precision = PEER_getAlpha(model)

residuals = t(PEER_getResiduals(model))
rownames(residuals) <- rownames(mat.ds)
colnames(residuals) <- colnames(mat.ds)

rownames(factors) <- c(colnames(gt_pca), paste0("PEER_",1:numcov))
colnames(factors) <- colnames(mat.ds)

residuals.ds <- residuals

	#png(paste0(loop.pop[i], '.expr.peer.clust.png'), width=8, height=8, res=150, units='in')
	#heatmap.2(as.matrix(residuals.ds), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
        #  trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(TISSUE, '\nexpr clustering'))
	#dev.off()

gz1 <- "pdui.peer.residuals.txt"
write.table(cbind(rownames(residuals), residuals), file=gz1, row.names=FALSE, col.names=c("id",colnames(residuals)), quote=FALSE, sep='\t')
rm(model,mat.ds,mat_impute,weights,precision,residuals,gz1)

# --------------------- prepare genotype matrix
cat("Load genotype_matrix.bed:\n")
gt_mat <- read.table("./Matrix_eQTL/genotype_matrix.bed",header=T,sep="\t",check.names=FALSE)
dim(gt_mat)
gt_mat.reorder <- gt_mat %>% dplyr::select("id",all_of(id_order))

cat("Output genotype:\n")
write.table(gt_mat.reorder,file="./Matrix_eQTL/Genotype_matrix.txt",quote=F,sep="\t",row.names=F,col.names=T)

rm(gt_mat)

# -------------------- prepare covariates matrix
factors.df <- cbind(rownames(factors),factors)
colnames(factors.df) <- c("id",colnames(factors))
factors.df <- as.data.frame(factors.df)
factors.reorder <- factors.df %>% dplyr::select("id",all_of(id_order))
covariate_file <- "./Matrix_eQTL/Covariate_matrix.txt"
write.table(factors.reorder, file=covariate_file, row.names=FALSE,quote=FALSE, sep='\t',col.names=T)
rm(factors,factors.df,id_order)
