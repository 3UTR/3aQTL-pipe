#!/bin/bash
# This is the sub-pipe of 3'aQTL-pipe, here in this script we will perform fine-mapping of 3'aQTL detected by Matrix-eQTL with SuSieR
# @Xudong Zou, zouxd@szbl.ac.cn
# 2022-03-30

# -- Usage function
script_name=$0
function usage(){
	echo "#=============================="
	echo "Default usage:"
	echo "#=============================="
	echo "bash $script_name -w 1000000 -p 0.1 -L 10 -V 0.2 -t 8"
	echo "Options:"
	echo "        -w  integer,setting the window size around aGenes for fine-mapping"
	echo "        -p  float, specify the minimum PIP for filtering fine mapped 3'aQTLs"
	echo "        -L  integer, specify the L value in susieR, default 10"
	echo "        -V  float, specify the variance used in susieR, default 0.2"
	echo "        -t  integer, setting threads to run susieR in parallel, default 8"
	echo "        -h  print the help information"
	exit 1
}

# define global variables from command parameters
currDir=`pwd`
sourceDir="./src"
PIP="0.1"
Variance="0.2"
L="10"
Threads="8"
window=""

while getopts :w:p:L:V:t:h opt
do
	case $opt in
		w)
			window="$OPTARG"
		;;
		p)
			PIP="$OPTARG"
		;;
		L)
			L="$OPTARG"
		;;
		V)
			Variance="$OPTARG"
		;;
		t)
			Threads="$OPTARG"
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
	echo "${currDir}/FineMapping/input not found!"
	exit
fi

if [ ! -d "${currDir}/FineMapping/output" ]
then
	echo "${currDir}/FineMapping/output not found!"
	exit
fi

# -- Main function --
function main(){
	echo "Running $script_name with the following parameters:"
	echo "*************************************************"
	echo "-w: $window"
	echo "-p: $PIP"
	echo "-L: $L"
	echo "-V: $Variance"
	echo "-t: $Threads"
	echo "*************************************************"
	echo "Start 3'aQTL fine mapping ..."
	date
	echo "Run fine-mapping analysis by susieR"
	run_fine_mapping $window $PIP $L $Variance $Threads
	echo "Done!"
	date
}


# -- Other functions --
function run_fine_mapping(){
	w=$1
	min_PIP=$2
	L=$3
	Var=$4
	threads=$5
	if [ ! -f "${currDir}/FineMapping/input/aGenes.loc_${w}.txt" ]
	then
		echo "File ${currDir}/FineMapping/input/aGenes.loc_${w}.txt not found!"
		exit
	fi

	if [ $threads -eq 1 ]
	then
		for gene in `cat ${currDir}/FineMapping/input/aGenes.loc_${w}.txt|cut -f1`
		do
			echo "Analyzing $gene"
			if [ -d "${currDir}/FineMapping/output/$gene" ]
			then
				GeneDir=${currDir}/FineMapping/output/$gene
				if [ -f "${GeneDir}/3aQTL.vcf" -a -f "${GeneDir}/expr.phen" ]
				then
					Rscript ${sourceDir}/SuSiE_Finemapping.R ${GeneDir} $L $Var $min_PIP &
					wait
				else
					echo "${gene}:File 3aQTL.vcf and expr.phen not found!"
					continue
				fi
			else
				echo "${gene} not exits!"
				continue
			fi
		done
		cd $currDir
	else
		if [ ! -d "${currDir}/FineMapping/output/tmp" ]
		then
			mkdir -p ${currDir}/FineMapping/output/tmp
		fi
		split -l $threads -d ${currDir}/FineMapping/input/aGenes.loc_${w}.txt ${currDir}/FineMapping/output/tmp/finemap_task_ &
		wait
		for task in `ls ${currDir}/FineMapping/output/tmp/finemap_task_*`
		do
			for gene in `cat $task |cut -f1`
			do
				echo "Analyzing $gene"
				if [ -d "${currDir}/FineMapping/output/$gene" ]
				then
					GeneDir=${currDir}/FineMapping/output/$gene
					if [ -f "${GeneDir}/3aQTL.vcf" -a -f "${GeneDir}/expr.phen" ]
					then
						Rscript ${sourceDir}/SuSiE_Finemapping.R ${GeneDir} $L $Var $min_PIP &
					else
						echo "${gene}:File 3aQTL.vcf and expr.phen not found!"
						continue
					fi
				else
					echo "${gene} not exits!"
					continue
				fi
			done
			wait
		done &
		wait
		cd $currDir
	fi

}

# - main
main
