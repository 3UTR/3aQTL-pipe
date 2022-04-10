#!/opt/app/languages/R-3.6.3/bin/Rscript
# 2022-01-08

# print usage
usage <- function() {
  cat(
    'usage: Rscript 3UTR_impute_peer.R known.cov.txt topN
')
}

library(dplyr)
library(peer)
library(impute)

args <- commandArgs(T)
known_cov_file <- args[1]
topN_pca <- args[2]
cat('Arguments:','\n','topN_pca:',topN_pca,'\n',
    'known_cov_file:',known_cov_file,'\n')

cat('Current directory:')
getwd()

# --------------- load known covariates -----------------## --------------------------------
# genotype pca
gtPCA_file <- "./Matrix_eQTL/genotype_pca.eigenvec"
gt_pca <- read.table(gtPCA_file,header=F,sep=" ",stringsAsFactors=F)
N <- as.integer(topN_pca) + 1
gt_pca$V1 <- NULL;gt_pca <- gt_pca[,1:N]
names(gt_pca) <- c("subject_id",paste0("PC_",1:(N-1)))
rm(N)
# process known covariates
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


cat('Dimension of gt_pca:',dim(gt_pca),'\n')
rm(known_cov,known_cov_file)

# --------------- load PDUI matrix
pdui_mat <- read.table("./Dapars2_res.all_chromosomes.txt", stringsAsFactors=FALSE, header=TRUE,sep="\t")
pdui_mat <- pdui_mat[,-c(2,3,4)]

pdui_mat.sel <- pdui_mat %>% dplyr::select(all_of(gt_pca$subject_id))
pdui_mat.sel <- as.matrix(pdui_mat.sel)
rownames(pdui_mat.sel) <- pdui_mat[,1]


#remove genes with more than 50% entries missing and individuals with more than 80% missing data
pdui_mat.sel <- pdui_mat.sel[, colMeans(is.na(pdui_mat.sel)) <= 0.8];pdui_mat.sel <-  pdui_mat.sel[rowMeans(is.na(pdui_mat.sel)) < 0.5,]
class(pdui_mat.sel) <- 'numeric'

# convert covariates' data.frame to matrix 
rownames(gt_pca) <- gt_pca$subject_id;gt_pca <- as.matrix(gt_pca[,-1])

# ----------------- run peer
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

x <- cbind(rownames(mat_impute$data),df_w)
y <- colnames(x)[-1]
colnames(x) <- c("Gene",y)

write.table(x,file="./Matrix_eQTL/pdui_mat.imputed_qnorm.txt",row.names=F,col.names=T,quote=F,sep="\t")
rm(x,y)
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

covariate_file <- "./Matrix_eQTL/pdui.peer.covariates.txt"
write.table(cbind(rownames(factors), factors), file=covariate_file, row.names=FALSE, col.names=c("id",colnames(factors)), quote=FALSE, sep='\t')

gz1 <- "pdui.peer.residuals.txt"
write.table(cbind(rownames(residuals), residuals), file=gz1, row.names=FALSE, col.names=c("id",colnames(residuals)), quote=FALSE, sep='\t')
rm(model,mat.ds,mat_impute,factors,weights,precision,residuals,covariate_file,gz1)
