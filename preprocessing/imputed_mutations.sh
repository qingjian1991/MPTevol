#!/bin/sh 
#############################
#Version 1.1     								                                    #
#2019-10-25								                                            #
#Authored by Qingjian Chen  from Sun Yat-sen University Cancer center			    #
#email: chenqj@sysucc.org.cn							                            #

################# #####
#

###########################Input file options+###############################
#reference
index=/data1/database/human/hg19/ucsc/hg19.ucsc.genome.fa
annovarpath=/data1/soft/Annovar/Annovar/
baseDir=/data1/soft/ExomePipe

if test ! -e stderr; then mkdir stderr;echo "mkdir stderr"; fi
if test ! -e tmp;then mkdir tmp; echo "mkdir tmp"; fi

#########################Path configuration#####################################

if test $# != 3 && test $# != 4
then
	echo "Help information:

<1> get allele frequency across samples;
	sh imputed_mutations.sh [paramater.txt] [samplename] [maf] {removed.sample.pattern}

	<1.1> Example paramter.text, including: SampleGroups Bamfile Tumor_ID.
	:
		-normal  Breast.Normal.recal.bam  Breast_normal
		-I  Breast_1.recal.bam  Breast_1
		-I  Breast_2.recal.bam  Breast_2
    
    the SampleGroups \"-normal\" indicates the normal samples. 

	<1.2> samplename: Patients ID.

	<1.3> maf: The annotated maf file.

	<1.4> removed.sample.pattern: samples to removed analysis.
	
    
<2> Details: From the annotated maf, we get the mutations within patients acorss different sequencing regions. Then we impute these mutations by gatk GetPileupSummaries. In the end, we get the imputed mutations in maf formate.

	"
	exit 1
fi

#Prepare the params across samples.

#samples=$(ls *bam| awk -F "-" '{print $1}'|sort|uniq)
#for ii in $samples
#do
#	echo $ii
#	ls ${ii}*.bam|sed "s/_dedup_bqsr.bam//g"| awk '{ OFS="\t"; print "-I",$1"_dedup_bqsr.bam",$1}'| perl -alne 'if($F[2] !~/N$/ ){print $_}else{print "-normal\t$F[1]\t$F[2]"  } ' >${ii}.param
#done

echo "STEP 1: GetPileupSummaries"

parafile=$1
samplename=$2
maf=$3

if test $#==4; then
	removed=$4
	params=$(cat $parafile|grep -v $removed )
else
	removed=""
	params=$(cat $parafile)
fi


echo "parafile: $parafile.  samplename: $samplename"

#vcfs=($(cat $parafile| awk '{print $2}'))
bams=($(echo $params| xargs -n 3 |awk '{print $2}'))
names=($(echo $params| xargs -n 3 |awk '{print $3}'))
num=$(expr ${#bams[@]} - 1 )
#get the number of normal. Note: The shell array is start from 0.
id_normal=$(echo $params| xargs -n 3| grep -n  '\-normal'|awk -F ":"  '{print $1}')
id_normal=$(expr $id_normal - 1)

if test -z $id_normal
then
	echo "There is no normal sample, using -normal to indicate normal samples".
fi

if test ! -e stderr/${samplename};then mkdir stderr/${samplename}; fi
if test ! -e ${samplename};then mkdir ${samplename}; echo "mkdir ${samplename}"; fi


###################################################################################################
### To get the high confidence mutations, we set the cutoff: AF >=0.05 & AD.1(ALternative AD) >=3 & DP >= 10
###
###################################################################################################

echo "------------<1>-------------------- Prepare input files."

#parafile=499556.param
#samplename="499556"
#maf=/data1/qingjian/Rproject/QiuMZ/isma/gastric/gastric.Mutect_Merged.maf
if test -z $removed
then
	cat $maf| awk -va=$samplename '$18 ~ a { OFS="\t"; print "chr"$2,$3,".",$7,$9,".","PASS","AF=0.1"}'|perl -alne 'if($F[3]=~/[ATCG]{1}/ && $F[4]=~/[ATCG]{1}/ ){ print $_ }'  |cat header - |bcftools view --types snps - -O z>${samplename}/${samplename}.merge.gnomAD.vcf.gz
else
	cat $maf| awk -va=$samplename -vb=$removed '$18 ~ a && $18 !~ b { OFS="\t"; print "chr"$2,$3,".",$7,$9,".","PASS","AF=0.1"}'|perl -alne 'if($F[3]=~/[ATCG]{1}/ && $F[4]=~/[ATCG]{1}/ ){ print $_ }'  |cat header - |bcftools view --types snps - -O z>${samplename}/${samplename}.merge.gnomAD.vcf.gz
fi

#bgzip ${samplename}/${samplename}.merge.gnomAD.vcf
tabix -p vcf ${samplename}/${samplename}.merge.gnomAD.vcf.gz

echo "------------<2>-------------------- run GetPileupSummaries"

date

for i in `seq 0 $num`
do
	echo "GetPileupSummaries: ${names[$i]}"
	gatk GetPileupSummaries \
	  -mmq 20	\
      -I ${bams[$i]} \
      -V ${samplename}/${samplename}.merge.gnomAD.vcf.gz \
      -L ${samplename}/${samplename}.merge.gnomAD.vcf.gz \
      -O ${samplename}/${names[$i]}.pileups.table 2>stderr/${names[$i]}.GetPileupSummaries.stderr & 
done

wait

date

#notice normal
normalname=${names[$id_normal]}
sed -n '3,$p' ${samplename}/${normalname}.pileups.table >${samplename}/${normalname}.pileups.noheader.table
echo "The normal sample is ${normalname}"

cmd=""

for i in `seq 0 $num`
do
	
if test $i -ne $id_normal
then
	echo ${names[$i]}
	sed -n '3,$p' ${samplename}/${names[$i]}.pileups.table >${samplename}/${names[$i]}.pileups.noheader.table
	less -S ${samplename}/${samplename}.merge.gnomAD.vcf.gz|grep -v "#"|join2filesByallMultKeyVOUTVTab.pl - ${samplename}/${names[$i]}.pileups.noheader.table 0,1 0,1|join2filesByallMultKeyVOUTVTab.pl - ${samplename}/${normalname}.pileups.noheader.table 0,1 0,1 |awk -va=${names[$i]} '{OFS="\t"; print $1,$2,$2,$4,$5,$11,$12, $17,$18, a }'| awk '$6!="-" && $8!="-" ' >${samplename}/${names[$i]}.pileups
	cmd="$cmd ${samplename}/${names[$i]}.pileups "
fi

done

#combied with maf files.
cat $cmd >${samplename}/${samplename}.pileups 

Rscript ${baseDir}/bin/getAlleleFreq_from_mafs.R --maf=$maf --Pileup=${samplename}/${samplename}.pileups  --samplename=${samplename}/${samplename}

echo "Done. ${samplename}"
