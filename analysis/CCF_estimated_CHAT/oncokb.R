#See oncokb sites in the tumors.

target = rbind( 
  read.table("oncokb/oncokb_biomarker_drug_associations_breast.tsv", header = T, sep = "\t") %>%
  mutate(Type = "Breast"),
  
  read.table("oncokb/oncokb_biomarker_drug_associations_lung.tsv", header = T, sep = "\t") %>%
    mutate(Type = "Lung"),
  
  read.table("oncokb/oncokb_biomarker_drug_associations_coad.tsv", header = T, sep = "\t") %>%
    mutate(Type = "Coad")
  
)

##################################################################################################


# Maf object with CCF information
data.type = "split1"
maf <- readMaf(mafFile = sprintf( "MesKit/meskit.%s.mutation.txt", data.type),
               ccfFile = sprintf("MesKit/meskit.%s.CCF.txt", data.type),
               clinicalFile  = sprintf("MesKit/meskit.%s.clinical.txt", data.type),
               refBuild = "hg19") 



onckgenes = list()

onckgenes$breast = target %>% filter(Type == "Breast") %>% pull(Gene) %>% unique()
onckgenes$lung = target %>% filter(Type == "Lung") %>% pull(Gene) %>% unique()
onckgenes$coad = target %>% filter(Type == "Coad") %>% pull(Gene) %>% unique()

onckmutations = list()

onckmutations$breast =  maf$Breast@data %>%
  filter(Hugo_Symbol  %in% onckgenes$breast)

onckmutations$lung = maf$Lung@data %>%
  filter(Hugo_Symbol  %in% onckgenes$lung)

onckmutations$coad = maf$Met1@data %>%
  filter(Hugo_Symbol  %in% onckgenes$coad)



#comment
## (1): PIK3CA_p.H1047R :public clonal mutations, the aviable site.
## (2): FGFR3  Missense_Mutation in Breast_3, VAF_adj = 0.310
onckmutations$breast


#comment
#ARID1A: public clonal mutations, truncating mutations

onckmutations$lung


#comment
#CDK12: frameshit mutations in Coad_1, VAF_adj = 0.5
onckmutations$coad


### other types 

### MSI

### ERBB2 amplications.
### ERBB2 is located in chr17, whereas the CNVs in Breast were in chr1q,chr10,chr16

###ERBB2 without amplications.
#(base) [qingjian@xu2 panCancer]$ less -S all_thresholded.by_genes.txt|grep -w ERBB2
#ERBB2	ENSG00000141736.9	17q12	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0    0

