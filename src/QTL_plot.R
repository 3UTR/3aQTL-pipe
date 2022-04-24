library(optparse)

option_list <- list(
	make_option(c("-s","--snp"),type="character",default="NA",action="store",help="specify a SNP"),
	make_option(c("-g","--gene"),type="character",default="NA",action="store",help="specify a gene"),
	make_option(c("-G","--genotype"),type="character",default="./Matrix_eQTL/Genotype_matrix.txt",action="store",help="specify the genotype matrix used in 3'aQTL mapping, default is ./Matrix_eQTL/Genotype_matrix.txt"),
	make_option(c("-P","--phenotype"),type="character",default="./Matrix_eQTL/Phenotype_matrix.txt",action="store",help="specify the phenotype matrix used in 3'aQTL mapping, default is ./Matrix_eQTL/Phenotype_matrix.txt")
		    )

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

#load Genotype matrix and Phenotype matrix
gt <- read.table(opt$genotype,header=T,sep="\t", check.names=FALSE)
pt <- read.table(opt$phenotype,header=T,sep="\t", check.names=FALSE)

rownames(gt) <- gt[,1]
rownames(pt) <- pt[,1]
gt <- gt[,-1]
pt <- pt[,-1]

snp <- as.character(opt$snp)
gene <- as.character(opt$gene)
geneName <- strsplit(gene,split="|",fixed=T)[[1]][2]


e1 = as.numeric(pt[which(rownames(pt)==gene),])
s1 = as.numeric(gt[which(rownames(gt)==snp),])

lm1 = lm(e1 ~ s1)
pdf(paste(snp, geneName,"pdf", sep="."))
boxplot(e1 ~ s1, lwd = 2, xaxt="n",xlab="Genotype",ylab="Normalized PDUI",main=paste(snp,gene,sep=" || "))
axis(1,at=c(1:3),labels=c("REF","HET","ALT"))
stripchart(e1 ~ s1, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = c(rgb(102,194,165,max=255),rgb(252,141,98,max=255),rgb(141,160,203,max=255)))
dev.off()
