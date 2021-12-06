#!/bin/sh 
#############Pipe line for WGS/WES/targeting sequencing data analysis################
#Version 1.1     								    #
#2019-10-25								    #
#Authored by Qingjian Chen  from Sun Yat-sen University Cancer center			    #
#email: chenqj@sysucc.org.cn							    #
#####Brief Introdution ####
#This data analysis pipeline were writing for processing raw sequencing data from   #
#WGS/WXS experiment. It contains several preprocessing steps and summary steps for  #
#evaluting data quality. Before using this shell script in you own system, you must #
#make sure the following softwares and tools were properly installed and configured,#
#
#see: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
#which are:									    #
#bwa,
#samtools,
#gatk4                                      #
####################################################################################
#####################© 2016 Qi Zhao All Rights Reserved#############################
####################################################################################

###########################Input file options+###############################
#reference
index=/data1/database/human/hg19/bwaIndex/genome.fa

# bed file for analysis copynumber 
#bedfile=/CLS/sysu_rj_1/database/hg19/exome_seq_bed/bed_exomeseq_V6/Exon_haploxbed

#gatk knownsitefile
knowfile1=/data1/database/human/hg19/forGATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
knowfile2=/data1/database/human/hg19/forGATK/1000G_phase1.indels.hg19.vcf
knowfile3=/data1/database/human/hg19/forGATK/dbsnp_138.hg19.vcf


#Target sequencing Interval file(bed format with 6 columns)
interval=/data1/database/human/hg19/exome_seq_bed/haplox/Exon_haploxbed.bed

#########################Path configuration#####################################

if test ! -e stderr
then
    mkdir stderr
    echo "mkdir stderr"
fi

pwd=`pwd`
TEMPDIR=${pwd}/tmp
#create a temp dir avoiding inadequate disk space in default tmp folder of yor linux system
if [ ! -d "./tmp" ]; then
mkdir $TEMPDIR
fi


############Engage mapping step##########################

Start="Mapping"

while getopts 'mb' opt
do
  case $opt in
    m) Start="MD" 
    echo "starting from MarkDuplicates"
    shift
      ;;
    b) Start="BQSR" 
    echo "starting from BQSR"
    shift
      ;;
    *)
    echo "Help information:
    
    1) when Input is the fq:        sh 1.mapping.all.sh samplename reads1 reads2 nthread

    2)when Input is the bam:
    2.1 bam is for MarkDuplicates:  sh 1.mapping.all.sh -m {samplename}_sorted.bam nthread
    2.2 bam is for BQSR:            sh 1.mapping.all.sh -b {samplename}_sort_dedup.bam nthread
    "
    exit 1
      ;;


  esac
done

if test $# = 2
then
	if test $Start = "BQSR" -o $Start = "MD"
	then
		echo "Input is bam"
		#get pure sample name information
	    bamfile=$1
		#thread num; default is 24
		nthread=$2
    else
		echo "Please set MD(-m) or BQSR(-b) when inputted file is bam"
		exit 1
    fi
elif [[ $# = 4 ]]; 
then
	#statements
	echo "Input is fq"
	#input fq/fq.gz file for paired end sequencing data 
	samplename=$1
	#thread num
	name[1]=$2
	name[2]=$3
	nthread=$4
	#reads1 fastq
	file1=${name[1]}
	#reads2 fastq
	file2=${name[2]}
else
	echo "Please checking input formate, The paramets number is : $#"
	exit 1
fi

#setting sample names.
case $Start in
      MD) samplename=${bamfile%%_sorted.bam};;
      BQSR) samplename=${bamfile%%_sort_dedup.bam};;  
esac

echo $Start
echo $samplename
echo $nthread

filename=$samplename


if test $Start = "Mapping"

then

	echo "${samplename} is mapping now"
	echo "${file1}"
	echo "${file2}"
	############################################step 1: mapping
	#
	
	if test ! -e ${samplename}_sorted.bam
	then
		date

		echo "bwa mem"
		bwa mem -t $nthread -M -R '@RG\tID:noID\tPL:ILLUMINA\tLB:noLB\tSM:'${samplename}'' $index $file1 $file2 \
		1>${filename}.sam 2>stderr/${filename}.bwa.stderr
	    
	    ###############Filter,convert & sorting, summary of mapping data####################
		date

		#bwa mem -t $nthread -M -R '@RG\tID:noID\tPL:ILLUMINA\tLB:noLB\tSM:'${samplename}'' $index $file1 $file2 \
		#2>stderr/${filename}.bwa.stderr | sambamba  view -S -f bam -t $nthread /dev/stdin 1>${filename}.bam

		samfile=${filename}.sam
		sample=$filename
		echo "sam to bam"
		sambamba  view -S -f bam -t $nthread $samfile 1>${filename}.bam  2>stderr/${filename}.samtoolsView.stderr

		sambamba sort  -t $nthread  --tmpdir=tmp -o ${filename}_sorted.bam ${filename}.bam
		date
		rm ${filename}.sam

	else
		echo "file exists, ${samplename}_sorted.bam"
	fi	

fi

############################################step 2: MarkDuplicatesSpark

# MarkDuplicatesSpark processing can replace both the MarkDuplicates and SortSam steps of the Best Practices single sample pipeline. After flagging duplicate sets, the tool automatically coordinate-sorts the records.

#provide reads-sorted bam rather than coordinated-sorted bam.
# 
#The tool is optimized to run on queryname-grouped alignments (that is, all reads with the same queryname are together in the input file). If provided coordinate-sorted alignments, the tool will spend additional time first queryname sorting the reads internally. 


if test $Start = "Mapping" -o $Start = "MD"
then
	if test ! -e ${samplename}_sort_dedup.bam
    then
		echo "MarkDuplicates"
		gatk MarkDuplicates \
		-I ${filename}_sorted.bam \
		-O ${filename}_sort_dedup.bam \
		-M ${filename}_marked_dup_metrics.txt \
		2>stderr/${samplename}.GATK_MarkDuplicates.stderr
        samtools index ${filename}_sort_dedup.bam 
	else
		echo "file exists, ${samplename}_sort_dedup.bam"
    fi
fi

############################################step 3: MarkDuplicatesSpark

echo "Base quality score recalibration"
       
if test ! -e ${samplename}_dedup_bqsr.grp
then
 	gatk  BaseRecalibrator \
	   -R $index \
	   -I ${samplename}_sort_dedup.bam \
	   -L $interval \
	   --known-sites $knowfile3 \
	   --known-sites $knowfile1 \
	   --known-sites $knowfile2 \
	   -O ${samplename}_dedup_bqsr.grp \
	   --verbosity INFO \
	   2>stderr/${samplename}.GATK_BaseRecalibrator.stderr
else
       echo "file exists, ${samplename}_dedup_bqsr.grp"
fi
       
date


echo "ApplyBQSR"

if test ! -e ${samplename}_dedup_bqsr.bam
then
     gatk  ApplyBQSR \
       -R $index \
       -I ${samplename}_sort_dedup.bam \
       -L $interval \
       -bqsr ${samplename}_dedup_bqsr.grp \
       -O ${samplename}_dedup_bqsr.bam \
       --tmp-dir tmp \
       2>stderr/${samplename}.GATK_ApplyBQSR.stderr
    date
else
	echo "file exists, ${samplename}_dedup_bqsr.bam"
fi

# rm the inter-median files.
if test -e ${samplename}_dedup_bqsr.bam
then
	rm ${samplename}_sorted.bam* ${samplename}_sort_dedup.bam* ${filename}.bam
fi

echo "Done: 1-mapping"
