#!/bin/bash
# use bedtools genomecov to convert bam files into bedgraph (*.wig)
# before running this script, make a text file (e.g. sample_list.txt) containing all samples (column 1) and their bam files (column 2)
# usage: bash bam2wig.sh sample_list.txt 8
echo "Start ..."
date
sampleList=$1
Threads=$2
currDir=`pwd`

if [ ! -f "$sampleList" ]
then
	echo "File $sampleList not found!"
	exit
fi

# extract the count of samples in sampleList
N=`cat $sampleList|wc -l`
if [ ! -d "${currDir}/tmp" ]
then
	mkdir -p ${currDir}/tmp
fi

if [ ! -d "${currDir}/wig" ]
then
	mkdir -p ${currDir}/wig
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
		sample=`echo $line |awk '{print $1}'`
		bam=`echo $line |awk '{print $2}'`
		echo $sample
		bedtools genomecov -ibam $bam -bga -split -trackline > ${currDir}/wig/${sample}.wig &
	done < $task
	wait
	echo "Task_$task_i finished!"
	(( task_i += 1 ))
done &
wait
echo "Done!"

rm ${currDir}/tmp/*
rmdir ${currDir}/tmp
date
