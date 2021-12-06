#!/bin/sh

#this is for 192.168.120.51
varpath=/data1/soft/VarScan.v2.3.9.jar
index=/data1/database/human/hg19/bwaIndex/genome.fa

#The following softwares are required.
#1) bedtools
#2）bam-readcount

if test $# != 3 & test $# != 4
then
	echo "
sh var_somaticSNP.sh [normal.bam] [tumor.bam] [sampleIDs] {bamfolder}
or
sh var_somaticSNP.sh [normal.bam] [tumor.bam] [sampleIDs]
	"
	exit 1
fi

if test ! -e stderr
then
    mkdir stderr
    echo "mkdir stderr"
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

out=$3

if test $# == 4
then
     folder=$4
else
     folder="."
fi

echo "sampleid: $out"
echo "normal: $normal"
echo "tumor: $tumor"

# remove zero coverage reads
echo "mpileup"
samtools mpileup -f $index -q 1 -B $folder/$normal $folder/$tumor | awk -F "\t" '$4 > 0 && $7 > 0' 1> ${out}.mpileup 2>stderr/${out}.mpileup

############################Running somatic mutations.+###############################

echo "somatic mutations"
#--strand-filter - If set to 1, removes variants with >90% strand bias [0]
time java -jar $varpath somatic ${out}.mpileup $out --strand-filter 1 --mpileup 1 --output-vcf 1 --min-var-freq 0.01  2>stderr/${out}.varscan.somatic

#filter
# --min-tumor-freq - Minimum variant allele frequency in tumor [0.10]
java -jar $varpath processSomatic ${out}.snp.vcf --min-tumor-freq 0.05  2>stderr/${out}.varscan.processSomatic.snp
java -jar $varpath processSomatic ${out}.indel.vcf --min-tumor-freq 0.1  2>stderr/${out}.varscan.processSomatic.indel

#Additional filter: identify and remove somatic mutations that are likely false positives due to alignment problems near indels
#--min-strands2 Minimum # of strands on which variant observed (1 or 2) [1]
#--min-var-freq  Minimum variant allele frequency threshold [0.20]
#--min-reads2    Minimum supporting reads for a variant [2]
#--min-coverage  Minimum read depth [10]
#--p-value	Default p-value threshold for calling variants [5e-02]

java -jar $varpath somaticFilter ${out}.snp.Somatic.hc.vcf \
--output-file  ${out}.snp.Somatic.hc.filter.vcf \
--indel-file ${out}.indel.Somatic.hc.vcf \
--min-strands2 2 \
--min-var-freq 0.05 \
--min-reads2 4 \
--min-coverage  15 \
--p-value 0.05 \
2>stderr/${out}.varscan.somaticFilter

#merge somatic snp and indels.
egrep -v "#" ${out}.indel.Somatic.hc.vcf | cat ${out}.snp.Somatic.hc.filter.vcf - > ${out}.snp_indel.Somatic.hc.vcf 


###########################Running Germline+###############################
 
#echo "germline mutations"

#merge somatic snp and indels.
egrep -v "#" ${out}.indel.Germline.hc.vcf | cat ${out}.snp.Germline.hc.vcf - > ${out}.snp_indel.Germline.hc.vcf 

## 3.1: Prepare a BED file from the high-confidence germline mutation VCF file to be used with bam-readcount:
egrep -hv "^#" ${out}.snp_indel.Germline.hc.vcf | awk 'OFS="\t" {print $1, $2-1, $2+1}' | sort -k1,1 -k2,2n | bedtools merge -i - >${out}.germline.hc.bed

## 3.2: Run bam-readcount:
time bam-readcount -f $index -q 20 -b 25 -l ${out}.germline.hc.bed -w 1 $normal >  ${out}.germline.hc.bamRC

# 3.3: Run fpfilter:
# --min-var-count		Minimum number of variant-supporting reads [4]
#--min-var-freq		Minimum variant allele frequency [0.05]
java -jar  $varpath fpfilter  ${out}.snp_indel.Germline.hc.vcf  ${out}.germline.hc.bamRC --output-file  ${out}.snp_indel.Germline.hc.fpfilterPassed.vcf --filtered-file  ${out}.snp_indel.Germline.hc.fpfilterFailed.vcf --min-var-count 5  --min-var-freq 0.20

###########################Running copynumber+###############################

#echo "copynumber"

# Run VarScan copynumber on normal and tumor mpileup output
#time java -jar $varpath copynumber ${out}.mpileup $out --mpileup 1 --output-vcf 1 2>stderr/${out}.varscan.copynumber

#Run VarScan copyCaller to adjust for GC content and make preliminary calls.
#java -jar $varpath copyCaller ${out}.copynumber --output-file ${out}.copynumber.called --output-homdel-file ${out}.copynumber.called.homdel 2>stderr/${out}.varscan.copyCaller

rm ${out}.mpileup 

echo "Done"
