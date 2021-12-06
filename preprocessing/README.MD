## Main
This fold contains scripts including mapping, mutation calling, imputed mutations, copy number alterations(CNAs) calling. 

### 1. mapping
mapping the sequencing reads to the human references (hg19)

**scripts**: 1.mapping.all.sh

**Usage**: 

    sh 1.mapping.all.sh [samplename] [reads1] [reads2] [nthread]

### 2. calling somatic mutation from Mutect2
Calling mutations from the bam by Mutect2.

**scripts**: 2.2.mutect2.sh

**Usage**: 
    
    sh 2.2.mutect2.sh [normal.bam] [tumor.bam] [sampleIDs] [Threads]


### 3. imputed mutations.
we impute the mutations that were present in a single patients from different regions. These mutations were stored in a maf file. The imputation were achieved by gatk *GetPileupSummaries*. 

**scripts**: imputed_mutations.sh

**Usage**: 

    sh imputed_mutations.sh [paramater.txt] [samplename] [maf] [removed.sample.pattern]

    <1> The paramter.text includes three columns: SampleGroups, Bamfile, and Tumor_ID. The SampleGroups **-normal** indicates the normal samples.
    
    The example paramter.text is as following:
                -normal  Breast.Normal.recal.bam  Breast_normal
                -I  Breast_1.recal.bam  Breast_1
                -I  Breast_2.recal.bam  Breast_2

    <2> samplename: Patients ID.

    <3> maf: The annotated maf file.

    <4> removed.sample.pattern: samples patterns that were removed for further analysis.



### 4. calling CNAs from sequenza
Calling CNAs from sequenza. Two supporting files, including **sequenza.R** and **extract_Cellularity_Ploidy_sequenza.R** were also provided.

**scripts**: run_sequenza_paired_WSX.sh

**Usage**: 
    
    sh run_sequenza_paired_WSX.sh [normal.bam] [tumor.bam] [sampleIDs]
    
    

### 5. calling somatic and germline mutations from varscan2

**scripts**: varscan_somaticSNP.sh

**Usage**: 
    
    sh varscan_somaticSNP.sh [normal.bam] [tumor.bam] [sampleIDs]
