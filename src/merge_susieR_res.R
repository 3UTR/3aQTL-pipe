args <- commandArgs(trailingOnly=TRUE)
basedir <- args[1] # the basic path of susieR analysis
aGene_file <- args[2]
setwd(basedir)
cat('Options:\n','basedir:',basedir,'\naGene_file:',aGene_file,'\n')
gene_list <- read.table(paste0("./input/",aGene_file),header=F,sep="\t")
gene_list <- gene_list[,1]
independent_snp_count <- c()

susie_df <- data.frame(locus_id=c(),variant_id=c(),pip=c(),cs=c(),cs_size=c(),cs_purity=c())
for(idx in 1:length(gene_list)){
	file_name <- paste0("./output/",gene_list[idx],"/3aQTL.SuSiE.txt")
	cat(file_name,"\n")
	if (file.exists(file_name)){
		df <- read.table(file_name,header=T,sep=" ")
		independent_snp_count[idx] <- dim(df)[1]

		if (dim(df)[1]>0){
			df$locus_id <- gene_list[idx]
			susie_df <- rbind(susie_df,df)
		}
	}else{
		independent_snp_count[idx] <- NA
	}
}

summary_susie <- data.frame(Gene=gene_list,Count=independent_snp_count)

write.table(susie_df,file="susieR_res.all_genes.txt",quote=F,row.names=F,sep="\t")
write.table(summary_susie,file="susieR_res.stat.txt",quote=F,row.names=F,sep="\t")
