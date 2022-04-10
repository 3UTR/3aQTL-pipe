#!/bin/bash
# Prepare input data for fine-mapping by susieR
# @Xudong Zou, zouxd@szbl.ac.cn
# 2022-03-30

# -- Usage function
script_name=$0
function usage(){
	echo "#=============================="
	echo "Default usage:"
	echo "#=============================="
	echo "bash $script_name -w 1000000 -q 0.05"
	echo "Options:"
	echo "        -w  integer,setting the window size around aGenes for fine-mapping"
	echo "        -q  float, specify the maximum of FDR to filtering significant 3'aQTL association,default 0.05"
	echo "        -h print the help information"
	exit 1
}

# define global variables from command parameters
currDir=`pwd`
window="1000000"
FDR="0.05"
sourceDir="./src"
while getopts :w:q:h opt
do
	case $opt in 
		w)
			window="$OPTARG"
		;;
		q)
			FDR="$OPTARG"
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
if [ ! -d "${currDir}/FineMapping/input" ]
then
	mkdir -p ${currDir}/FineMapping/input
fi

if [ ! -d "${currDir}/output/FineMapping/output" ]
then
	mkdir -p ${currDir}/FineMapping/output
fi

if [ ! -d "${currDir}/Matrix_eQTL" ]

then
	echo "No Matrix-eQTL output found!"
	exit
fi

# -- Main function --
function main(){
	echo "Running $script_name with the following parameters:"
	echo "*************************************************"
	echo "-w: $window"
	echo "-q: $FDR"
	echo "*************************************************"
	date
	echo "Prepare input for susieR"
	prepare_input $window $FDR
	echo "Done!"
	date
}


# -- Other functions --
function prepare_input(){
	w=$1
	fdr=$2
	if [ ! -f "${currDir}/Matrix_eQTL/3UTR_location.txt" ]
	then
		echo "File ${currDir}/Matrix_eQTL/3UTR_location.txt not found!"
		exit
	fi

	if [ ! -f "${currDir}/Matrix_eQTL/Genotype_matrix.txt" ]
	then
		echo "File ${currDir}/Matrix_eQTL/Genotype_matrix.txt not found!"
		exit
	fi

	echo "Prepare unique aGenes set..."
	python ${sourceDir}/prepare_susieR_uniqGene_location.py --utr_loc_file ${currDir}/Matrix_eQTL/3UTR_location.txt \
		--aQTL_map ${currDir}/Matrix_eQTL/Cis_3aQTL_all_control_gene_exprs.txt \
		--extend_size $w \
		--Max_FDR $fdr \
		--outdir ${currDir}/FineMapping/input \
		--output aGenes.loc_${w}.txt &
	wait

	echo "Prepare SNP files in bed format..."
	python ${sourceDir}/genotype_2_bed.py --genotype ${currDir}/Matrix_eQTL/Genotype_matrix.txt \
		--out_bed ${currDir}/FineMapping/input/Genotype_matrix.bed \
		--out_header ${currDir}/FineMapping/input/Header.txt &
	wait

	if [ -f "${currDir}/FineMapping/input/Genotype_matrix.bed" ]
	then
		sort -k1,1 -k2,2n ${currDir}/FineMapping/input/Genotype_matrix.bed > tmp.bed &
		wait
		mv tmp.bed ${currDir}/FineMapping/input/Genotype_matrix.bed &
		wait
	else
		echo "File ${currDir}/FineMapping/input/Genotype_matrix.bed not found!"
		exit
	fi
#  create a unique workspace for each gene in "aGenes.loc_${w}.txt", and generate a "expr.phen" for each gene
	if [ ! -f "${currDir}/Matrix_eQTL/Phenotype_matrix.txt" ]
	then
		echo "File ${currDir}/Matrix_eQTL/Phenotype_matrix.txt not found!"
		exit
	fi
	echo "Make a directory and generate a expr.phen for each gene"
	for gene in `cat ${currDir}/FineMapping/input/aGenes.loc_${w}.txt| cut -f1`
	do
		mkdir -p ${currDir}/FineMapping/output/$gene
		cat ${currDir}/Matrix_eQTL/Phenotype_matrix.txt | awk -v aGENE=$gene -F"\t" 'BEGIN{OFS="\t"} {if(NR==1){for (i=2;i<NF;++i) SAMPLES[i]=$i} if ($1==aGENE){ for(i=2;i<NF;++i) print SAMPLES[i],SAMPLES[i],$i}}' > ${currDir}/FineMapping/output/$gene/expr.phen &
		wait
	done

	echo "Select SNPs around a window of ${w}bp of the gene and generate 3aQTL.vcf in the gene's directory"
	while read line
	do
		gene=`echo $line|awk '{print $1}'`
		loc=`echo $line|awk '{print $2}'`
		cd ${currDir}/FineMapping/output/$gene
		CHR=${loc%:*}
		COORD=${loc#*:}
		S=${COORD%-*}
		E=${COORD#*-}
		echo -e "$CHR\t$S\t$E" > gene_loc.bed
		cat ${currDir}/FineMapping/input/Header.txt > 3aQTL.vcf
		bedtools intersect -a ${currDir}/FineMapping/input/Genotype_matrix.bed -b gene_loc.bed -wa |cut -f4- >> 3aQTL.vcf &
		wait
		rm gene_loc.bed
	done < ${currDir}/FineMapping/input/aGenes.loc_${w}.txt

	cd ${currDir}
}

# - main
main
