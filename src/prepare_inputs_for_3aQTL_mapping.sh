#!/bin/bash

# -- usage function
script_name=$0
function usage(){
	echo "#=============================="
	echo "Default usage:"
	echo "bash $script_name -g <vcf_list.txt> -p <Dapars2_res.all_chromosomes.txt> -c <known_covariates.txt> -s <sample_list.txt> -m <0.05> -n <5>"
	echo "Options:"
	echo "        -g  text file,input a text file contains VCF file(s), default=vcf_list.txt"
	echo "        -p  text file,input the merged DaPars2 results, default=Dapars2_res.all_chromosomes.txt"
	echo "        -c  text file,input the known covariates like age and gender, default=NA"
	echo "        -s  text file,input a text file contains the list of samples, default=sample_list.txt"
	echo "        -m  float,minor allele frequency threshold for selecting common genetic variants, default=0.05"
	echo "        -n  integer,the top N genotype PCA components to be used as covariates, default=5"
	echo "        -h  print the help information"
	exit 1

}
# define global variables from command parameters
currDir=`pwd`
VCFLIST="vcf_list.txt"
APA_RES="Dapars2_res.all_chromosomes.txt"
KNOWN_COV="NA"
SAMPLES="sample_list.txt"
MAF="0.05"
TOP_N="5"

while getopts g:p:c:s:m:n:h opt
do
	case $opt in
		g)
			if [ ! -f "$OPTARG" ];then echo "File $OPTARG not found";exit 0;fi
			VCFLIST="$OPTARG"
		;;
		p)
			if [ ! -f "$OPTARG" ];then echo "File $OPTARG not found";exit 0;fi
			APA_RES="$OPTARG"
		;;
		c)
			if [ ! -f "$OPTARG" ];then echo "File $OPTARG not found";exit 0;fi
			KNOWN_COV="$OPTARG"
		;;
		s)
			if [ ! -f "$OPTARG" ];then echo "File $OPTARG not found";exit 0;fi
			SAMPLES="$OPTARG"
		;;
		m)
			MAF="$OPTARG"
		;;
		n)
			TOP_N="$OPTARG"
		;;
		h)
			echo "Help message:"
			usage
		;;
		:)
			echo "The option -$OPTARG requires an argument."
			exit 1
		;;
		?)
			echo "Invalid option: $OPTARG"
			usage
			exit 2
		;;
	esac

done

# -- Basic settings
if [ ! -d "${currDir}/tmp" ]
then
	mkdir -p ${currDir}/tmp
	echo "Create a directory called tmp/"
fi

if [ ! -d "${currDir}/Matrix_eQTL" ]
then
	mkdir -p ${currDir}/Matrix_eQTL
	echo "Create a directory called Matrix_eQTL/"
fi

# -- Main function --
function main(){
	date
	echo "Start..."
	echo "Convert VCF file into 012 format..."
	generate_gt_matrix $VCFLIST $MAF $SAMPLES
	echo "PCA analysis on genotype by PLINK1.9 ..."
	PCA_on_genotype $VCFLIST $SAMPLES $MAF
	echo "Curate phenotype matrix, genotype matrix, and covariate matrix for Matrix-eQTL"
	curate_pheno_geno_covariates $APA_RES $KNOWN_COV $TOP_N
	echo "extracting 3UTR location and SNP location:"
	snp_and_3utr_location $APA_RES
	echo "Done!"
	date
}

function snp_and_3utr_location(){
	dapars2_res=$1
	python ./src/extract_SNP_location.py --genotype_bed ./Matrix_eQTL/genotype_matrix.bed --output ${currDir}/Matrix_eQTL/snp_location.txt &
	wait
	python ./src/extract_3UTR_location.py --dapars_res $dapars2_res --output ${currDir}/3UTR_location.txt &
	wait
	echo "Two location files: snp_location.txt and 3UTR_location.txt are generated!"
}
function curate_pheno_geno_covariates(){
	dapars2_res=$1
	known_cov=$2
	topN=$3
	Rscript ./src/curate_pheno_geno_covariates.R -p $dapars2_res -c $known_cov -n $topN
}
# PCA analysis on genotype by PLINK 1.9
function PCA_on_genotype(){
	vcfList=$1
	keep_inds=$2
	maf=$3

	if [ ! -f "$vcfList" ]
	then
		echo "$vcfList not exits!"
		exit
	fi

	if [ ! -f "$keep_inds" ]
	then
		echo "File $keep_inds not found!"
		exit
	fi


	for vcf in `cat $vcfList`
	do
		if [ ! -f "$vcf" ]
		then
			echo "$vcf not exists!"
			exit
		else
			filename=`basename $vcf .gz`
			vcf_body=${filename%.*}

			plink --vcf $vcf --const-fid --out ${currDir}/tmp/${vcf_body}.plink &
			wait

			echo "${vcf_body}.plink.bed ${vcf_body}.plink.bim ${vcf_body}.plink.fam" >> ${currDir}/tmp/merge_list.txt
		fi
	done
	sleep 10

	cd ${currDir}/tmp
	N=`cat merge_list.txt|wc -l`

	if [ $N -gt 1 ]
	then
		tmp=`cat merge_list.txt|head -n 1|awk '{print $1}'`
		firstVCF=${tmp%.*}
		cat merge_list.txt|tail -n+2 > tmp.txt
		cat tmp.txt > merge_list.txt
		rm tmp.txt

		echo "Merging genotype of multi-chromosomes into one file: merged_plink.* ..."
		plink --bfile $firstVCF --merge-list merge_list.txt --out merged_plink --allow-extra-chr &
		wait

	else
		cd $currDir
		vcf=`cat $vcfList`
		filename=`basename $vcf .gz`
		vcf_body=${filename%.*}
		mv ${currDir}/tmp/${vcf_body}.plink.bed ${currDir}/tmp/merged_plink.bed
		mv ${currDir}/tmp/${vcf_body}.plink.bim ${currDir}/tmp/merged_plink.bim
		mv ${currDir}/tmp/${vcf_body}.plink.fam ${currDir}/tmp/merged_plink.fam

	fi

	# extract selected samples in plink file
	cd ${currDir}
	cat $keep_inds | awk '{print "0",$1}' > ${currDir}/tmp/keep.list
	cd ${currDir}/tmp
	plink --bfile merged_plink --keep keep.list --geno 0.02 --hwe 0.000001 --maf $maf --make-bed --out merged_plink_QC &
	wait
	sleep 20

	# pca analysis
	plink --bfile merged_plink_QC --indep-pairwise 50 5 0.2 --out merged_plink_QC &
	wait
	plink --bfile merged_plink_QC --extract merged_plink_QC.prune.in --pca 30 --out genotype_pca &
	wait

	# mv pca results to input file for furhter analysis
	if [ -f "genotype_pca.eigenvec" ]
	then
		cp genotype_pca.eigenvec ${currDir}/Matrix_eQTL
		echo "move genotype_pca.eigenvec to ${currDir}/Matrix_eQTL"

		cd $currDir
	else
		echo "genotype_pca.eigenvec not found!"
		exit
	fi

}

# recode genotype into 012 format from VCF
function generate_gt_matrix(){
	vcfList=$1
	maf=$2
	keep_inds=$3
	
	cat $keep_inds |cut -f1 > ${currDir}/tmp/keep_inds.txt

	for vcf in `cat $vcfList`
	do
		if [ ! -f "$vcf" ]
		then
			echo "$vcf not exists!"
			exit
		else
			last_suffix=${vcf##*.}
			if [ $last_suffix = "gz" ]
			then
				filename=`basename $vcf`
				tmp=${filename%.*}
				vcf_body=${tmp%.*}
				echo "Extract genotype from $vcf ..."
				vcftools --gzvcf $vcf --out ${currDir}/tmp/${vcf_body}.gt_filtering --remove-filtered-all --keep $keep_inds --maf $maf --max-missing-count 10 --extract-FORMAT-info GT &
				
				echo "Extract allele frequence from $vcf ..."
				vcftools --gzvcf $vcf --out ${currDir}/tmp/${vcf_body}.gt_filtering --remove-filtered-all --keep $keep_inds  --maf $maf --max-missing-count 10 --freq &
				
				wait
				echo "$vcf Done!"
			elif [ $last_suffix = "vcf" ]
			then
				filename=`basename $vcf`
				vcf_body=${filename%.*}
				echo "Extract genotype from $vcf ..."
				vcftools --vcf $vcf --out ${currDir}/tmp/${vcf_body}.gt_filtering --remove-filtered-all --keep $keep_inds --maf $maf --max-missing-count 10 --extract-FORMAT-info GT &
				
				echo "Extract allele frequence from $vcf ..."
				vcftools --vcf $vcf --out ${currDir}/tmp/${vcf_body}.gt_filtering --remove-filtered-all --keep $keep_inds --maf $maf --max-missing-count 10 --freq &
				wait
				echo "$vcf Done!"
			else
				echo "Unrecognized file format."
				exit

			fi

		fi
	done &
	wait

	for frq in `ls ${currDir}/tmp/*.frq`
	do
		filename=`basename $frq`
		sample=${filename%.*}
		python ./src/recode_with_012.py --frq $frq --GT ${currDir}/tmp/${sample}.GT.FORMAT --output ${currDir}/tmp/${sample}.GT.bed &
		wait
	done &
	wait

	cat ${currDir}/tmp/*.GT.bed |head -n 1 > ${currDir}/Matrix_eQTL/genotype_matrix.bed &
	wait
	for bed in `ls ${currDir}/tmp/*.GT.bed`
	do
		cat $bed |tail -n+2 >> ${currDir}/Matrix_eQTL/genotype_matrix.bed &
		wait
	done

	echo "Genotype matrix has been generated: ./Matrix_eQTL/genotype_matrix.bed"
	rm ${currDir}/tmp/*

}

# -- main
main
