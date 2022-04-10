#!/bin/bash

# use samtools flagstat to obtain statistics of a bam file
# use extract_read_depth.py to extract total number of mapped reads
# a text file "sample_list.txt" contains all samples and corresponding bam files are required
# usage: bash extract_sequencing_depth.sh sample_list.txt [/path/to/wig] 8 [output_file_name]
echo "Start ..."
date
sampleList=$1
wigDir=$2
Threads=$3
outFile="wigFile_and_readDepth.txt"
currDir=`pwd`

if [ ! -z "$4" ]
then
	outFile=$4
fi

if [ ! -f "${currDir}/$sampleList" ]
then
	echo "File $sampleList not found!"
	exit
fi

if [ ! -d "$wigDir" ]
then
	echo "$wigDir not exists!"
	exit
fi

N=`cat $sampleList | wc -l`
if [ ! -d "${currDir}/tmp" ]
then
	mkdir -p ${currDir}/tmp
fi

if [ $Threads -lt $N ]
then
	split -l $Threads $sampleList -d ${currDir}/tmp/task_ &
	wait
else
	cp $sampleList ${currDir}/tmp/task_00
fi

task_i=1
for task in `ls ${currDir}/tmp/task_*`
do
	while read line
	do
		sample=`echo $line | awk '{print $1}'`
		bam=`echo $line | awk '{print $2}'`
		echo $sample
		samtools flagstat -@ 2 $bam > ${currDir}/tmp/${sample}.flagstat &
	done < $task
	wait
	echo "Task_$task_i finished!"
	(( task_i += 1 ))
done &
wait

echo "samtools flagstat finished!"
echo "Start extracting total aligned reads in each sample ..."
python ${currDir}/src/extract_read_depth.py --sample_list $sampleList --path_flagstat ./tmp --path_wig $wigDir --output wigFile_and_readDepth.txt &
wait
echo "Done!"
rm ${currDir}/tmp/*
rmdir ${currDir}/tmp
date
