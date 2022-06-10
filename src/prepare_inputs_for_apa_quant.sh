#!/bin/bash
# This script takes bam files, a gene annotation (bed), a ID mapping file (between refseq and gene symbol),
# and generate bedgraph files, read depth, 3'UTR reference region
# Two input files represent the gene annotation file (in bed format) and the ID mapping between Refseq transcript ID and gene symbol are required for executing this script 
# @Xudong Zou, zouxd@szbl.ac.cn
# 2022-04-21

# -- Usage function
script_name=$0
function usage(){
	echo "#=============================="
	echo "Default usage:"
	echo "#=============================="
	echo "bash $script_name -s <sample_list> -g <geneAnnotation> -r <refseqID_to_GeneName> -t <N_threads> -c <Coverage_cutoff>"
	echo "Options:"
	echo "        -s  text file,input a text file contains all samples (column 1) and corresponding bam files (column 2)"
	echo "        -g  text file,provide a RefSeq gene annotation file extracted from UCSC"
	echo "        -r  text file,provide a file list the ID mapping between RefSeq transcript and gene name"
	echo "        -t  integer,specify the number of threads used to parallelly running Dapars2,default=8"
	echo "        -c  integer,define the threshold of reads coverage of the alternative APA site,default=15"
	echo "        -o  file name,specify a name for the configure file which will be used by DaPars2"
	echo "        -h  print the help information"
	exit 1
}

# define global variables from command parameters
currDir=`pwd`
BAMlist="sample_list.txt"
GeneAnno=""
RefIDmap=""
P=8
Cutoff_Cov=15
Config="Dapars2_running_configure.txt"

while getopts :s:g:r:t:c:o:h opt
do
	case $opt in 
		s)
			if [ ! -f "$OPTARG" ];then echo "File not exists"; exit;fi
			BAMlist="$OPTARG"
		;;
		g)
			if [ ! -f "$OPTARG" ];then echo "File $OPTARG not found!";exit;fi
			GeneAnno="$OPTARG"
		;;
		r)
			if [ ! -f "$OPTARG" ];then echo "File $OPTARG not found!";exit;fi
			RefIDmap="$OPTARG"
		;;
		t)
			P="$OPTARG"
		;;
		c)
			Cutoff_Cov="$OPTARG"
		;;
		o)
			Config="$OPTARG"
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


# -- Main function --
function main(){
	echo "Running $script_name with the following parameters:"
	echo "*************************************************"
	echo "-s: $BAMlist"
	echo "-g: $GeneAnno"
	echo "-r: $RefIDmap"
	echo "-t $P"
	echo "-c: $Cutoff_Cov"
	echo "-o: $Config"
	echo "*************************************************"
	echo "Start preparing inputs for Dapars2 ..."
	date
	echo "Convert bam to bedgraph format ..."
	run_bam2bedgraph $BAMlist $P
	echo "Counting total aligned reads by samtools flagstat ..."
	run_samtools_flagstat $BAMlist $P
	echo "Generating the 3' UTR reference ..."
	generate_3utr_reference $GeneAnno $RefIDmap
	echo "Generating the wigFile_and_readDepth.txt ..."
	generate_wigFileList_with_readDepth $BAMlist
	echo "Generating a configure file for Dapars2 ..."
	generate_configure_for_dapars2 $Cutoff_Cov $P
	wait
	echo "Done!"
	date
}


# -- Other functions --
function generate_configure_for_dapars2(){
	N_cov=$1
	N_threads=$2
	python ./src/generate_configure_for_Dapars2.py --annotation_3utr ${currDir}/refseq_3utr_annotation.bed \
		--wigFile_depth ${currDir}/wigFile_and_readDepth.txt \
		--coverage_threshold $N_cov \
		--threads $N_threads \
		--out_config_name ${currDir}/Dapars2_running_configure.txt &
	wait
	echo "Generate file Dapars2.allSamples_joint.configure.txt in ${currDir}/input/"
}
function generate_wigFileList_with_readDepth(){
	bamList=$1
	if [ ! -f "$bamList" ]
	then
		echo "File $bamList not found!"
		exit
	fi
	python ./src/extract_read_depth.py --sample_list $bamList --path_wig ${currDir}/wig --output ${currDir}/wigFile_and_readDepth.txt &
	wait
	echo "Generate file wigFile_and_readDepth.txt"

}
function generate_3utr_reference(){
	gene_anno=$1
	refID2Symbol=$2
	python ./src/DaPars_Extract_Anno.py -b ${gene_anno} -s $refID2Symbol -o ${currDir}/refseq_3utr_annotation.bed &
	wait
	echo "Generate refseq_3utr_annotation.bed"
}

function run_bam2bedgraph(){
	bamList=$1
	N_jobs=$2
	if [ ! -d "${currDir}/tmp" ]
	then
		mkdir -p ${currDir}/tmp
	fi

	if [ ! -d "${currDir}/wig" ]
	then
		mkdir -p ${currDir}/wig
	fi

	if [ ! -f "$bamList" ]
	then
		echo "File $bamList not found!"
		exit
	fi

	N_samples=`cat $bamList|wc -l`
	echo "$N_samples bam files waiting for processing."
	if [ $N_jobs -lt $N_samples ]
	then
		split -l $N_jobs $bamList -d ${currDir}/tmp/task_ &
		wait
	else
		echo -e "Number of parallel jobs exceeds the total number of tasks.\n$N_samples threads will be used!"
		cat $bamList > $currDir/tmp/task_00 &
		wait
	fi

	i_task=1
	for task in `ls ${currDir}/tmp/task_*`
	do
		echo "Start subtask ${i_task}..."
		while read line
		do
			sample=`echo $line | awk '{print $1}'`
			bam=`echo $line | awk '{print $2}'`
			echo $sample
			bedtools genomecov -ibam ${bam} -bga -split -trackline > ${currDir}/wig/${sample}.wig &
		done < $task
		wait
		echo "Subtask $i_task finished!"
		(( i_task += 1 ))
	done &
	wait

	echo "$N_samples bam files processed"
	rm ${currDir}/tmp/*
	rmdir ${currDir}/tmp
	date
}

function run_samtools_flagstat(){
	bamList=$1
	N_jobs=$2
	if [ ! -d "${currDir}/tmp" ]
	then
		mkdir -p ${currDir}/tmp
	fi

	if [ ! -f "$bamList" ]
	then
		echo "File $bamList not found!"
		exit
	fi
	N_samples=`cat $bamList|wc -l`
	echo "$N_samples bam files waiting for processing."
	if [ $N_jobs -lt $N_samples ]
	then
		split -l $N_jobs $bamList -d ${currDir}/tmp/task_ &
		wait
	else
		echo -e "Number of parallel jobs exceeds the total number of tasks.\nN_samples will be used!"
		cat $bamList > ${currDir}/tmp/task_00 &
		wait
	fi

        i_task=1
        for task in `ls ${currDir}/tmp/task_*`
        do
                echo "Start subtask ${i_task}..."
                while read line
                do
			sample=`echo $line | awk '{print $1}'`
			bam=`echo $line | awk '{print $2}'`
                        echo $sample
                        samtools flagstat -@ 2 ${bam} > ${currDir}/tmp/${sample}.flagstat &
                done < $task
                wait
                echo "Subtask $i_task finished!"
                (( i_task += 1 ))
        done
}


# -- run main
main
