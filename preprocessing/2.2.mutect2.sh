#!/bin/sh 
#############Pipe line for WGS/WES/targeting sequencing data analysis################
#Version 1.1     								                                    #
#2019-10-25								                                            #
#Authored by Qingjian Chen  from Sun Yat-sen University Cancer center			    #
#email: chenqj@sysucc.org.cn							                            #

# mutect2(gatk version 4.1.4.0 )

################# Mutect2: Call somatic SNVs and indels via local assembly of haplotypes#####
#

###########################Input file options+###############################
#reference
index=/data1/database/human/hg19/bwaIndex/genome.fa
#population mutations
af_only_gnomad=/data1/database/human/hg19/forGATK/af-only-gnomad.raw.sites.hg19.vcf.gz

#PON
pon_sample=/data1/database/human/TCGA.PON/hg19/gatk4_mutect2_4136_pon.hg19.PASS.vcf.gz


#load gatk, using gatk 4.1.4.0
#module avil python/2.7.16-gcc-5.4.0
#gatk=/CLS/sysu_rj_1/data/chenqj/Soft/gatk-4.1.4.0/gatk

#create folder
if test ! -e stderr; then mkdir stderr;echo "mkdir stderr"; fi
if test ! -e tmp;then mkdir tmp; echo "mkdir tmp"; fi

#########################Path configuration#####################################

if test $# != 4
then
	echo "Help information:

	gatk mutect2 for somatic mutations calling;
	sh 2.2.mutect2.sh [{normal}_dedup_bqsr.bam] [{tumor}_dedup_bqsr.bam] [sampleIDs] [threads:1,2,5,6,12,20]

	"
	exit 1
fi

#partition=$4  #divide the genomes into x regions.
partition=$4
echo "Split regions into $partition parts."

#chr Region
#cut the whole genomes into [10,18 or 24] regions using ***gatk SplitIntervals*****.
intervalsfolder=/data1/qingjian/Soft/bedfile/Exon_haploxbed_hg19_${partition}/
prefix=$(ls $intervalsfolder)

if test ! -e $intervalsfolder
then
    echo "no folder: $prefix $intervalsfolder"
    exit 1
else
	echo "cutting genomes into $partition parts"
fi

#test normal and tumor bam names.
if test ! -e $1
then
	normalbam="${1}_dedup_bqsr.bam"
	tumorbam="${2}_dedup_bqsr.bam"
else
	normalbam=$1
	tumorbam=$2
fi

samplename=$3
normalname=${normalbam%%_dedup_bqsr.bam}

if test ! -e tmp/${samplename};then mkdir tmp/${samplename}; fi
if test ! -e stderr/${samplename};then mkdir stderr/${samplename}; fi

#rm tmp/${samplename}/*
#rm stderr/${samplename}/*

#Steps 1. cut the genome into 24 parts and run mutect2 seperately.

echo "Start for $samplename--$normalname (tumor-normal) calling"

echo "STEP 1: calling mutect2 into different regions"


date

for i in $prefix
do
	id=$(echo $i|sed s/-scattered.interval_list//)

    gatk --java-options "-Xmx4g" \
     Mutect2 \
	-R $index \
	-I $tumorbam \
	-I $normalbam \
	-normal $normalname \
	-O tmp/${samplename}/${samplename}.${id}.mutect2.vcf.gz \
	-L "${intervalsfolder}${i}" \
	--native-pair-hmm-threads 2 \
	--germline-resource $af_only_gnomad \
	--panel-of-normals $pon_sample \
	2>stderr/${samplename}/${samplename}.mutect2.${id}.stderr &

done 

wait # waiting for calling variants.

date

echo "STEP 2: concatenate vcfs into one files"

filenum=$(ls tmp/${samplename}/${samplename}.*.mutect2.vcf.gz.tbi|wc -l)

if test $filenum = $partition
then
	ls tmp/${samplename}/${samplename}.*.mutect2.vcf.gz > ${samplename}.vcfs.list
	gatk MergeVcfs \
	 -I ${samplename}.vcfs.list \
	 -O ${samplename}.mutect2.vcf.gz \
	 2>stderr/${samplename}.MergeVcfs.mutect2.stderr
else
	echo "the sub-vcf is imcomplete, please check the file. $filenum"
	exit 1
fi

echo "STEP 3: get stats from separete tasks"

if test $partition > 1 
then
	ls tmp/${samplename}/${samplename}.*.mutect2.vcf.gz.stats| awk '{print "-stats",$1}' > ${samplename}.stats.list

	gatk --java-options "-Xmx2g" MergeMutectStats \
	 --arguments_file ${samplename}.stats.list \
	 -O ${samplename}.mutect2.vcf.gz.stats \
	 2>stderr/${samplename}.MergeMutectStats.mutect2.stderr

fi

echo "Start for ${samplename} filtering variants"

echo "STEP 4: FilterMutectCalls"

gatk FilterMutectCalls \
   -R $index \
   -V ${samplename}.mutect2.vcf.gz \
   -O ${samplename}.filtered.vcf.gz \
   2>stderr/${samplename}.FilterMutectCalls.stderr

echo "STEP 5: SelectVariants"

gatk SelectVariants \
-R $index \
-V  ${samplename}.filtered.vcf.gz \
-select "vc.isNotFiltered()" \
-O ${samplename}.filtered.PASS.vcf.gz \
2>stderr/${samplename}.SelectVariants.stderr


echo "2.2 -Mutect2: done $samplename"
