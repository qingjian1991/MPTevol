#!/bin/sh
# initial setting 

#this is for 192.168.120.51
index=/data1/database/human/hg19/bwaIndex/genome.fa

wig50=/data1/database/human/hg19/sequenza/hg19_genome_gc50.wig.gz

# WES data
captureRegions=/data1/qingjian/Soft/bedfile/Exon_haploxbed_hg19.bed

if test $# != 3
then
	echo "sh run_freec_paired_WXS.sh [normal.bam] [tumor.bam] [sampleIDs]"
	exit 1
fi

#test normal and tumor bam names.
if test ! -e $1
then
	normal="${1}_dedup_bqsr.bam"
	tumor="${2}_dedup_bqsr.bam"
else
	normal=$1
	tumor=$2
fi

samplename=$3

echo "sampleid: $out"
echo "normal: $normal"
echo "tumor: $tumor"


#sequenza-utils https://sequenza-utils.readthedocs.io/en/latest/guide.html

#1)Generate GC reference file in 50bp windows.
#need run only-one time. 
#sequenza-utils gc_wiggle --fasta $index -w 50 -o hg19_genome_gc50.wig.gz

#2)Normal and tumor pileup files

echo "##################################################################"
echo "######   Step1: pileup file.       ##################             "
echo "##################################################################"

date
if test ! -e ${samplename}.tumor.pileup.gz
then
	samtools mpileup -l $captureRegions -q 20 $tumor | gzip >${samplename}.tumor.pileup.gz & 
fi

if test ! -e ${samplename}.normal.pileup.gz
then
	samtools mpileup -l $captureRegions -q 20 $normal | gzip >${samplename}.normal.pileup.gz &
fi
date

wait

echo "##################################################################"
echo "######   Step2: seqz       ##################             "
echo "##################################################################"

date

if test ! -e ${samplename}.seqz.gz
then
	sequenza-utils bam2seqz \
	    --normal ${samplename}.normal.pileup.gz \
	    --tumor ${samplename}.tumor.pileup.gz \
	    --fasta $index \
	    -gc ${wig50} --output ${samplename}.seqz.gz --pileup
fi

date

#2.2) Binning seqz, reduce memory
if test ! -e ${samplename}.bin50.seqz.gz
then
	sequenza-utils seqz_binning --seqz ${samplename}.seqz.gz --window 50 \
	    -o ${samplename}.bin50.seqz.gz
fi

echo "##################################################################"
echo "######   Step3: R visulaization       ##################             "
echo "##################################################################"

Rscript sequenza.R --seqz.file=${samplename}.bin50.seqz.gz  --sample.id=${samplename} 

echo "Done"
