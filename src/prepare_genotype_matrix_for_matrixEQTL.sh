#!/bin/bash


# define global variables from command parameters
currDir=`pwd`
VCFLIST=$1
MAF=$2
SAMPLES=$3
REF=$4
sourceDir="./src"

# -- Basic settings
if [ ! -d "${currDir}/tmp" ]
then
	mkdir -p ${currDir}/tmp
fi

if [ ! -d "${currDir}/Matrix_eQTL" ]
then
	mkdir -p ${currDir}/Matrix_eQTL
fi

# -- Main function --
function main(){
	date
	echo "Start..."
	echo "Convert VCF file into 012 format..."
	generate_gt_matrix $VCFLIST $MAF $SAMPLES $REF
	echo "PCA analysis on genotype by PLINK1.9 ..."
	PCA_on_genotype $VCFLIST $SAMPLES $MAF
	echo "Done!"
	date
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
	ref_build=$4
	
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
	done	

	sleep 20
	for frq in `ls ${currDir}/tmp/*.frq`
	do
		filename=`basename $frq`
		sample=${filename%.*}
		if [ ! -f "${currDir}/tmp/${sample}.GT.FORMAT" ]
		then
			echo "No GT.FORMAT file for sample: $sample"
			exit
		else
			python ${sourceDir}/recode_with_012.py --frq $frq --GT ${currDir}/tmp/${sample}.GT.FORMAT --Reference $ref_build --output ${currDir}/tmp/${sample}.GT.bed &
			wait
		fi
	done &
	wait
	sleep 10
	cat ${currDir}/tmp/*.GT.bed |head -n 1 > ${currDir}/Matrix_eQTL/genotype_matrix.bed &
	for bed in `ls ${currDir}/tmp/*.GT.bed`
	do
		cat $bed |tail -n+2 >> ${currDir}/Matrix_eQTL/genotype_matrix.bed
	done

	echo "Genotype matrix has been generated: ./Matrix_eQTL/genotype_matrix.bed"
	rm ${currDir}/tmp/*

}

# -- main
main
