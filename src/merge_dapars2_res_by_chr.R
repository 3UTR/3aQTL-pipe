#!/opt/app/languages/R-4.1.0/bin/Rscript

# merge Dapars2 output by chromosome
# -- global variable
args <- commandArgs(T)
dir <- args[1] # Prefix of Dapars2 output directory,default is "Dapars2_out"
prefix <- args[2]
bamList <- args[3] # bam file list file
chromList <- args[4] # chromosome list file: e.g. chrList.txt
# -- functions

load_dapars2_res <- function(chromosome,new_header){
	input_file <- paste0(dir,"_",chromosome,"/",prefix,"_result_temp.",chromosome,".txt")
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
chrs_vec <- as.character(chrs_list)
rm(chrs_list)
if(! "chr" %in% chrs_vec){
	chrs_vec <- paste0("chr",chrs_vec)
}

res.df <- data.frame()

for(chr in chrs_vec){
	temp.df <- load_dapars2_res(chr,col_names)
	print(paste(chr,dim(temp.df)[1],sep=":"))
	res.df <- rbind(res.df,temp.df)
}

dim(res.df)
write.table(res.df,file="Dapars2_res.all_chromosomes.txt",quote=F,sep="\t",row.names=F)


