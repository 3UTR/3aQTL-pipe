args <- commandArgs(T)


#load Genotype matrix and Phenotype matrix
gt <- read.table("./Matrix_eQTL/Genotype_matrix.txt",header=T,sep="\t")
pt <- read.table("./Matrix_eQTL/Phenotype_matrix.txt",header=T,sep="\t")

rownames(gt) <- gt[,1]
rownames(pt) <- pt[,1]
gt <- gt[,-1]
pt <- pt[,-1]

snp <- as.character(args[1])
gene <- as.character(args[2])
geneName <- strsplit(gene,split="|",fixed=T)[[1]][2]


e1 = as.numeric(pt[which(rownames(pt)==gene),])
s1 = as.numeric(gt[which(rownames(gt)==snp),])

lm1 = lm(e1 ~ s1)
pdf(paste(snp, geneName,"pdf", sep="."))
boxplot(e1 ~ s1, lwd = 2, xaxt="n",xlab="Genotype",ylab="Normalized PDUI",main=paste(snp,gene,sep=" || "))
axis(1,at=c(1:3),labels=c("REF","HET","ALT"))
stripchart(e1 ~ s1, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = c(rgb(102,194,165,max=255),rgb(252,141,98,max=255),rgb(141,160,203,max=255)))
dev.off()
