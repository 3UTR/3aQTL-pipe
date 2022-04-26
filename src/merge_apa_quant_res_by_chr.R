# merge Dapars2 output by chromosome
library(optparse)
# -- global variable
option_list <- list(
		    make_option(c("-d","--dir_prefix"), type = "character", default = "Dapars2_out",
				action = "store", help = "Specify the directory prefix of DaPars2 output"),
		    make_option(c("-f", "--file_prefix"), type = "character", default = "Dapars2",
				action = "store", help = "Specify the file prefix of DaPars2 output"),
		    make_option(c("-s", "--sample_list"), type = "character", default = "sample_list.txt",
				action = "store", help = "A file contains the sample list"),
		    make_option(c("-c", "--chr_list"), type = "character", default = "chrList.txt",
				action = "store", help = "A file contains the chromosome list"),
		    make_option(c("-o", "--output"), type = "character", default = "Dapars2_res.all_chromosomes.txt",
				action = "store", help = "Specify the output file name")
		    )

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))
dir_pre <- opt$dir_prefix
file_pre <- opt$file_prefix
bamList <- opt$sample_list
chromList <- opt$chr_list
outFile <- opt$output
cat("dir_pre:",dir_pre,"\nfile_pre:",file_pre,"\nchromList:",chromList,"\n")
# -- functions
load_dapars2_res <- function(chromosome,new_header){
	input_file <- paste0(dir_pre,"_",chromosome,"/",file_pre,"_result_temp.",chromosome,".txt")
	dap_res <- read.table(input_file,header=T, sep="\t")
	names(dap_res) <- new_header

	return(dap_res)
}


# -- main

# load samples
dat <- read.table(bamList,header=F)
sample_list <- as.character(dat$V1)
col_names <- c("Gene","fit_value","Predicted_Proximal_APA","Loci",sample_list)
chrs_list <- read.table(chromList,header=F)
chrs_vec <- as.character(chrs_list$V1)
rm(chrs_list)
if(substr(chrs_vec[1],1,3)!="chr"){
	chrs_vec <- paste0("chr",chrs_vec)
}

chrs_vec
res.df <- data.frame()

for(chr in chrs_vec){
	temp.df <- load_dapars2_res(chr,col_names)
	print(paste(chr,dim(temp.df)[1],sep=":"))
	res.df <- rbind(res.df,temp.df)
}

dim(res.df)
write.table(res.df,file=outFile,quote=F,sep="\t",row.names=F)


