#Analysis MRS data using MesKit.

library(MesKit)
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
#-----------Prepare mutation data ---------------

#How to do:
#For the CCF data, the CCF data from CHAT are filtered with Low-quality data, so we limit the maf data with in the filtered CCF data.

#sampABlist is in CCF_adjust.R

sampABlist= readRDS("data/sampABlist.rds")

#The adjust CCF is stored in sampABlist$entire

#Tumor_Sample_Barcode: 
#Tumor_ID: tumor site, including primary, met, liver-M and so on.
#Patient_ID: Patients
#Tumor_Sample_Label: names for label.

#type is merge or split or split.all

#merge: all samples are labeled from one patients.
#split: split all sites according to their sites.
#split.1: combined Coad and related metastasis.

data.type = "split" 

if(data.type == "split1"){
  
  #cmobined all CRC relative tumor into one.
  clinical.info = data.frame(
    Tumor_Sample_Barcode = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                              paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                              paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7)) )) %>%
    mutate(
      Tumor_ID = mapply(function(x) x[1], str_split(Tumor_Sample_Barcode, "_") ),
      Patient_ID = c( rep("Breast", 5), rep("Lung", 5), rep("Met1", 17), "Uterus", "Met1", rep("Uterus", 4)  ),
      Tumor_Sample_Label = Tumor_Sample_Barcode
    )
  
}else if(data.type == "split"){
  #default, each samples from the same organ.
  clinical.info = data.frame(
    Tumor_Sample_Barcode = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                              paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                              paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7)) )) %>%
    mutate(
      Tumor_ID = mapply(function(x) x[1], str_split(Tumor_Sample_Barcode, "_") ),
      Patient_ID =Tumor_ID,
      Tumor_Sample_Label = Tumor_Sample_Barcode
    )
  
}else if(data.type == "merge"){
  
  #combined: we combined all samples into one patient. We will see the different between different sites.
  
  clinical.info = data.frame(
    Tumor_Sample_Barcode = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                              paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                              paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7)) )) %>%
    mutate(
      Tumor_ID = c( rep("Breast", 5), rep("Lung", 5), rep("Coad", 5), 
                    rep("OveryLM", 5), rep("OveryRM", 6),
                    "UterusM", "UterusM", "UterusM",rep("Uterus", 4)  ),
      Patient_ID = "P1",
      Tumor_Sample_Label = Tumor_Sample_Barcode
    )
  
}else if(data.type = "split2"){
  #split UterusM into 2 groups.
  
  #cmobined all CRC relative tumor into one.
  clinical.info = data.frame(
    Tumor_Sample_Barcode = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                              paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                              paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7)) )) %>%
    mutate(
      Tumor_ID = mapply(function(x) x[1], str_split(Tumor_Sample_Barcode, "_") ),
      Patient_ID = c( rep("Breast", 5), 
                      rep("Lung", 5), 
                      rep("Coad", 5),
                      rep("OveryLM", 5), 
                      rep("OveryRM", 6),
                      "UterusM1","UterusM","UterusM1",
                      rep("UterusM", 4) ),
      Tumor_Sample_Label = Tumor_Sample_Barcode
    )
  
}


#sampABlist$entire %>% select(starts_with("Lung_1") ) %>% head()

#get all adjust maf.
ccf.data = sampABlist$entire %>%
  mutate(mutation_id = str_c(chr , from , ref , alt, sep = ":")) %>%
  dplyr::select(mutation_id, ends_with("ccf")  & !starts_with("merge") ) %>%
  pivot_longer(
    cols =  c(ends_with("ccf") ),
    names_to = "Tumor_Sample_Barcode",
    names_pattern = "(.*)",
    values_to = "CCF"
  ) %>%
  mutate(Tumor_Sample_Barcode = str_remove(Tumor_Sample_Barcode, "ccf"))


ccfsd.data = sampABlist$entire %>%
  mutate(mutation_id = str_c(chr , from , ref , alt, sep = ":")) %>%
  dplyr::select(mutation_id, ends_with("ccfSD")  & !starts_with("merge") ) %>%
  pivot_longer(
    cols =  c(ends_with("ccfSD") ),
    names_to = "Tumor_Sample_Barcode",
    names_pattern = "(.*)",
    values_to = "CCF_Std"
  ) %>%
  mutate(Tumor_Sample_Barcode = str_remove(Tumor_Sample_Barcode, "ccfSD"))

ccf.data = inner_join( ccf.data, ccfsd.data ) %>%
  filter(CCF >0) %>%
  left_join(
    clinical.info %>% dplyr::select(Tumor_Sample_Barcode, Patient_ID)
  )
 
ccf.data = ccf.data %>%
  cbind(
    str_split(ccf.data$mutation_id, ":", simplify = T)  %>% as.data.frame() %>% setNames(c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2"))
        ) %>%
  mutate(Chromosome = str_remove(Chromosome, "chr"))


#maf data.
GC01 = readRDS("data/GC01.rds")

GC01 = GC01 %>% mutate(Chromosome =  str_c("chr", Chromosome) ) %>%
  mutate(mutation_id = str_c(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")) %>%
  filter(mutation_id %in% ccf.data$mutation_id ) %>%
  mutate(Chromosome = str_remove(Chromosome, "chr"))

system("mkdir MesKit")

write.table(clinical.info, file = sprintf("MesKit/meskit.%s.clinical.txt", data.type), sep = "\t", quote = F, row.names = F)
write.table(GC01, file = sprintf( "MesKit/meskit.%s.mutation.txt", data.type), sep = "\t", quote = F, row.names = F)
write.table(ccf.data, file = sprintf("MesKit/meskit.%s.CCF.txt", data.type), sep = "\t", quote = F, row.names = F)




#between regions.
colorScale <- c("#3C5488FF", "#00A087FF", "#F39B7fFF",
                "#8491B4FF","#E64B35FF","#4DBBD5FF",
                "#E41A1C", "#377EB8", "#7F0000",
                "#35978f", "#FC8D62", "#2166ac",
                "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#8DD3C7", "#6E016B" ,
                "#BEBADA", "#e08214", "#80B1D3",
                "#d6604d", "#ffff99", "#FCCDE5",
                "#FF6A5A", "#BC80BD", "#CCEBC5" ,
                "#fb9a99", "#B6646A", "#9F994E", 
                "#7570B3" , "#c51b7d" ,"#66A61E" ,
                "#E6AB02" , "#003c30", "#666666")


cols_samples = setNames(
  colorScale[31:36],
  nm = c("Breast","Coad","Lung","OveryLM","OveryRM","UterusM")
)



############ ------------------------------------------------------------------

###### see driver mutations.

# Maf object with CCF information
data.type = "split1"
maf <- readMaf(mafFile = sprintf( "MesKit/meskit.%s.mutation.txt", data.type),
               ccfFile = sprintf("MesKit/meskit.%s.CCF.txt", data.type),
               clinicalFile  = sprintf("MesKit/meskit.%s.clinical.txt", data.type),
               refBuild = "hg19",
               ccf.conf.level = 0.90
               ) 

#maf1 is split1


#5. Mutational landscape

# Driver genes of CRC collected from [IntOGen] (https://www.intogen.org/search) (v.2020.2)
#driverGene.file <- system.file("extdata/", "IntOGen-DriverGenes_COREAD.tsv", package = "MesKit")
#driverGene <- as.character(read.table(driverGene.file)$V1)

driverGene = read.delim("data/IntOGen-Drivers-Cancer_Genes.tsv", header = T) %>%
  filter(CANCER_TYPE %in% c("BRCA","COREAD","LUAD", "LUSC") ) %>%
  pull(SYMBOL) %>% unique()

driverGene = driverGene[!driverGene %in% c("BRCA1")]

mut.class <- classifyMut(maf, class =  "SP", patient.id = 'Breast')
head(mut.class)

#plotMutProfile(maf, class =  "SP", geneList = driverGene, use.tumorSampleLabel = TRUE, removeEmptyCols = FALSE)

#plotMutProfile(maf, class =  "CS", geneList = driverGene, use.tumorSampleLabel = TRUE,  removeEmptyCols = FALSE)

pdf("MesKit/MutProfile.pdf", width = 10, height = 6)

plotMutProfile(maf, class =  "SPCS", geneList = driverGene, use.tumorSampleLabel = TRUE, removeEmptyCols = FALSE)

dev.off()


#Get the summary of each mutations.

#get data
maf_input = subMaf(maf, mafObj = FALSE, use.tumorSampleLabel = TRUE) 
#get mutation classifications.
maf_class = classifyMut(maf, patient.id = NULL, class = "SPCS", classByTumor = FALSE)


maf_merge = list()

for(i in names(maf_input) ){
  message(i)
  maf_merge[[i]] = maf_input[[i]] %>%
    mutate(Mut_ID = str_c(Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":" ) ) %>%
    left_join(
      maf_class[[i]]
    ) %>%
    dplyr::select(
      Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Mutation_Type, Patient_ID, Tumor_ID, Variant_Classification, VAF  
    )
    
}


maf_merge <- purrr::reduce(maf_merge, rbind)

cols_types = setNames(c("#00A087FF", "#3C5488FF", "#8491B4FF", "#F39B7FFF", "#E64B35FF", "#4DBBD5FF")
, nm = c("Public_Clonal","Public_Subclonal","Shared_Clonal","Shared_Subclonal","Private_Clonal","Private_Subclonal") )

pdf("MesKit/MutNumber.pdf", width = 10, height = 4)

  maf_merge %>%
    mutate(Mutation_Type = factor(Mutation_Type, levels = names(cols_types) ) ) %>%
    ggplot() + theme_classic() +
    geom_bar(aes(x = Tumor_Sample_Barcode , fill = Mutation_Type) ) +
    facet_grid( ~ Patient_ID, scales = "free_x", space = "free_x"  ) +
    scale_fill_manual(
      values = cols_types
    ) +
    labs(y = "Mutation Counts") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 0.5)
    )
  
dev.off()


maf_merge %>%
  group_by(
    Tumor_Sample_Barcode, Patient_ID
  ) %>%
  summarise(
    num = n()
  ) %>%
  group_by(Patient_ID) %>%
  summarise(
    mean = mean(num),
    max = max(num),
    min = min(num)
    
  )

#Patient_ID  mean   max   min
#<chr>      <dbl> <int> <int>
#1 Breast      151.   182   127
#2 Coad        280.   471   220
#3 Lung        174.   345   105
#4 OveryLM     237.   278   195
#5 OveryRM     202.   248   183
#6 UterusM     175    266   115



#Proportion of driver mutations across different groups.

mutation_group = data.frame(
  Mutation_Type = c("Public_Clonal","Public_Subclonal","Shared_Clonal","Shared_Subclonal","Private_Clonal","Private_Subclonal"),
  Type = c("Shared_Clonal",NA,"Shared_Clonal",NA,"Private_Clonal","Private_Subclonal")
)

maf_merge1 = maf_merge %>%
  mutate(
    is.driver = ifelse( (Hugo_Symbol %in% driverGene) & ( !Variant_Classification %in% c("Silent", "3'Flank", "IGR", "Intron", "RNA") ) , TRUE, FALSE)
  ) %>%
  left_join(mutation_group) %>%
  filter(!is.na(Type))

maf_merge1 %>%
  ggplot(aes(x = Type, y = VAF) ) + theme_classic() +
  geom_boxplot()

maf_merge1 %>%
  group_by(Type) %>%
  summarise(num = n(), median = median(VAF))


maf_merge_sum = maf_merge1 %>%
  filter( !is.na(Type) ) %>%
  group_by(Tumor_ID, is.driver, Type) %>%
  summarise(num = n()) %>%
  group_by( Tumor_ID, Type ) %>%
  mutate(num_total = sum(num)) %>%
  filter(is.driver) %>%
  mutate(prop = num/num_total) %>%
  mutate(Type = factor(Type, levels = c("Shared_Clonal","Private_Clonal","Private_Subclonal") )) 

pdf("MesKit/Driver.Prop.pdf", width = 8, height = 4)

maf_merge_sum %>%
  ggplot(aes(x = Tumor_ID , y = prop, fill = Type) ) + ggpubr::theme_pubr() +
  geom_bar(stat = "identity", position =  position_dodge(width = 0.90)) +
  labs(x = NULL, y = "Ratio of mutated driver mutations") +
  scale_fill_manual(
    values = cols_types[levels(maf_merge_sum$Type) ]
  )

dev.off()

#See Ka/Ks between different groups.

getKaKs = function(df, VAF_cutoff = 0.05){
  
  data(list = sprintf("submod_%s", "13r_3w"), package = "dndscv")
  
  mutations = df %>%
    filter(VAF >= VAF_cutoff) %>%
    select(c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2"))
  
  kaks = dndscv::dndscv( mutations = mutations,
                     max_muts_per_gene_per_sample = Inf,
                     max_coding_muts_per_sample = Inf,
                     #gene_list = intersect( Genes_Covered, genes),
                     sm = submod_13r_3w
  )
  
  kaks$globaldnds
}


maf_list =  maf_merge1 %>%
  split( paste( maf_merge1$Tumor_ID, maf_merge1$Type, sep = ":" ) )

maf_KaKs = lapply(maf_list, getKaKs)


KaKs_data = list()

for(i in names(maf_KaKs) ){
  KaKs_data[[i]] = maf_KaKs[[i]] %>%
    mutate(type = i)
}


KaKs_data = purrr::reduce(KaKs_data, rbind)
  
KaKs_data = KaKs_data %>%
  mutate( 
    Tumor_ID = mapply(function(x) x[1], str_split(type, ":") ),
    Type = mapply(function(x) x[2], str_split(type, ":") )
  )


pdf("MesKit/Driver.KaKs.pdf", width = 8, height = 4)

KaKs_data %>%
  filter(name %in% c("wall")) %>%
  mutate(name = factor(name, levels = c("wall") )) %>%
  mutate(Type = factor(Type, levels = c("Shared_Clonal","Private_Clonal","Private_Subclonal") ))  %>%
  ggplot(aes(x = Tumor_ID , y = mle, fill = Type) ) + ggpubr::theme_pubr() +
  geom_bar(stat = "identity", position =  position_dodge(width = 0.90)) +
  #geom_linerange(aes(ymin = cilow, ymax = cihigh), position =  position_dodge(width = 0.90) ) +
  geom_hline( yintercept = 1, linetype = 2, size = 1) +
  labs(x = NULL, y =  latex2exp::TeX("Dn/Ds ($\\omega_{all}$)") )  +
  #facet_grid( ~ name) +
  scale_fill_manual(
    values = cols_types[levels(maf_merge_sum$Type) ]
  ) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14)
  )
  
dev.off()

save(maf, driverGene, maf_merge, KaKs_data, colorScale, file = "MesKit/Plot.data.rda")


#Measurement of Neutral evolution.

library(neutralitytestr)

#tt = neutralitytest(VAFselection)
#plot_all(tt)

neutral = maf_merge %>% 
  #filter(!Mutation_Type %in% c("Public_Clonal") ) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    num = sum(VAF >= 0.1),
    R2 = as.numeric(neutralitytestr::neutralitytest(VAF, fmin = 0.1, fmax = 0.22  )$rsq[1]), 
    R2_pval = as.numeric( neutralitytestr::neutralitytest(VAF, fmin = 0.1, fmax = 0.22)$rsq[2])
  ) %>%
  mutate(
    Patient_ID = mapply(function(x) x[1],  str_split(Tumor_Sample_Barcode, "_") )
  )

pdf("MesKit/NeutralTest.pdf", width = 5, height = 3.5)

neutral %>%
  ggplot(aes(x = Patient_ID, y = R2, col = Patient_ID, shape = Patient_ID) ) + 
  ggpubr::theme_classic2() +
  geom_violin() +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = NULL, y = latex2exp::TeX("$R^{2}$")    ) +
  scale_color_manual( values = cols_samples ) +
  geom_hline( yintercept = 0.98, linetype = 2, size = 1) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 15),
    axis.text.x = element_text(hjust = 1, angle = 45)
  )

dev.off()



#6. Measurement of intra-tumor heterogeneity

#within tumors.
scores = mathScore(maf, min.vaf = 0.05) %>%
  purrr::reduce(rbind)

scores %>%
  ggplot(aes(x = Patient_ID, y = MATH_Score, col = Patient_ID) ) + theme_classic2() +
  geom_violin() +
  geom_jitter(width = 0.2) +
  theme(
    legend.position = "none"
  )


pdf("MesKit/ITH.MATHscores.pdf", width = 5, height = 3.5)

scores %>%
  ggplot(aes(x = Patient_ID, y = MATH_Score, col = Patient_ID, shape = Patient_ID) ) + theme_classic2() +
  geom_violin() +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = NULL, y = "MATH scores"  ) +
  scale_color_manual( values = cols_samples ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 15),
    axis.text.x = element_text(hjust = 1, angle = 45)
  )

dev.off()





scores.ccfs = ccfAUC(maf, use.tumorSampleLabel = TRUE, min.ccf = 0.1)

scores.ccfs = lapply(scores.ccfs, function(x) x$AUC.value ) %>%
  purrr::reduce(rbind)

pdf("MesKit/ITH.AUCccf.pdf", width = 5, height = 3.5)

scores.ccfs %>%
  ggplot(aes(x = Patient_ID, y = AUC, col = Patient_ID, shape = Patient_ID) ) + theme_classic2() +
  geom_violin() +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = NULL, y = "AUC (ccf)"  ) +
  scale_color_manual( values = cols_samples ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 15),
    axis.text.x = element_text(hjust = 1, angle = 45)
  )

dev.off()

######################################################################
## Between Regions different, Fst and Nei's distance
# calculate the Fst of brain metastasis from V402
scores.Fst = calFst(maf, plot = TRUE, use.tumorSampleLabel = TRUE, 
       withinTumor  = FALSE, number.cex = 10)

Fst.data = list()

for(i in names(scores.Fst) ){
  Fst = scores.Fst[[i]]$Fst.pair
  Fst = Fst[lower.tri(Fst)]
  
  Fst.data[[i]] = data.frame(
    Fst = Fst,
    tumor = i
  )
}

Fst.data = purrr::reduce(Fst.data, rbind)

pdf("MesKit/ITH.Fst.pdf", width = 5, height = 3.5)

Fst.data %>%
  ggplot(aes(x = tumor, y = Fst, col = tumor, shape = tumor) ) + theme_classic2() +
  geom_violin() +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = NULL, y = "Fst"  ) +
  scale_color_manual( values = cols_samples ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 15),
    axis.text.x = element_text(hjust = 1, angle = 45)
  )

dev.off()

##Nei's distance.

scores.Nei = calNeiDist(maf, plot = TRUE, use.tumorSampleLabel = TRUE, 
                    withinTumor  = FALSE, number.cex = 10)

Nei.data = list()

for(i in names(scores.Nei) ){
  Nei = scores.Nei[[i]]$Nei.dist
  Nei = Nei[lower.tri(Nei)]
  
  Nei.data[[i]] = data.frame(
    Nei = Nei,
    tumor = i
  )
}

Nei.data = purrr::reduce(Nei.data, rbind)

pdf("MesKit/ITH.Nei.distance.pdf", width = 5, height = 3.5)

Nei.data %>%
  ggplot(aes(x = tumor, y = Nei, col = tumor, shape = tumor) ) + theme_classic2() +
  geom_violin() +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = NULL, y = "Nei's distance"  ) +
  scale_color_manual( values = cols_samples ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 15),
    axis.text.x = element_text(hjust = 1, angle = 45)
  )

dev.off()

#########################################7. Metastatic routes inference

# Maf object with CCF information
data.type = "split1"
maf <- readMaf(mafFile = sprintf( "MesKit/meskit.%s.mutation.txt", data.type),
               ccfFile = sprintf("MesKit/meskit.%s.CCF.txt", data.type),
               clinicalFile  = sprintf("MesKit/meskit.%s.clinical.txt", data.type),
               refBuild = "hg19",
               ccf.conf.level = 0.95
) 


JSI.res <- calJSI(maf, patient.id = 'Met1', pairByTumor = TRUE, min.ccf = 0.02, 
                  use.adjVAF = TRUE, use.indel = FALSE, use.tumorSampleLabel = TRUE)
names(JSI.res)


ccf.list <- compareCCF(maf, patient.id = 'Met1', pairByTumor = TRUE, min.ccf = 0.02,
                       use.adjVAF = TRUE, use.indel = FALSE)

#Met1_C_OL <- ccf.list$Met1$`Coad-OveryLM`
# visualize via smoothScatter R package
#graphics::smoothScatter(matrix(c(Met1_C_OL[, 3], Met1_C_OL[, 4]),ncol = 2),
#                        xlim = c(0, 1), ylim = c(0, 1),
#                        colramp = colorRampPalette(c("white", RColorBrewer::brewer.pal(9, "BuPu"))),
#                        xlab = "P", ylab = "BM")

## show driver genes
#gene.idx <- which(Met1_C_OL$Hugo_Symbol %in% driverGene) 
#points(Met1_C_OL[gene.idx, 3:4], cex = 0.6, col = 2, pch = 2)
#text(Met1_C_OL[gene.idx, 3:4], cex = 0.7, pos = 1,
#     Met1_C_OL$Hugo_Symbol[gene.idx])
#title("V402 JSI = 0.0974026", cex.main = 1.5)

#maf_merge1
library(ggrepel)

#add driver information.

plotCCF = function(ccf.list , JSI.res, maf_merge1){
  
  plist = list()
  
  for(i in names(ccf.list) ){
    
    message(i)
    
    pairs = as.character( str_split(i, "-", simplify = T))
    
    JSI = JSI.res$JSI.pair[pairs[1],pairs[2]]
    
    #merge driver info
    Met1_C_OL = ccf.list[[i]] %>%
      left_join( maf_merge1 %>% 
                   mutate(Mut_ID = str_c(Chromosome, Start_Position, Reference_Allele , Tumor_Seq_Allele2, sep = ":") ) %>%
                   dplyr::select(Mut_ID, is.driver) %>%
                   unique.data.frame()
      ) 
    
    plist[[i]] =  Met1_C_OL %>%
      ggplot(aes_string(x = pairs[1] , y = pairs[2]) ) + theme_classic2() +
      stat_density2d(aes(fill = ..density..^2), geom = "tile", contour = FALSE, n = 100) +
      scale_fill_continuous(low = "white", high = "red") +
      geom_point(alpha = 0.5) +
      geom_point(col = "blue", data = Met1_C_OL %>% filter(is.driver), size = 3, shape = 5 ) +
      geom_label_repel(aes(label = Hugo_Symbol), data = Met1_C_OL %>% filter(is.driver),
                       nudge_x = 0.1, nudge_y = 0.1
      ) +
      labs(title = sprintf("JSI = %.3f", JSI)  ) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  return( plist)
  
}



plist =  plotCCF(ccf.list = ccf.list$Met1, JSI.res,  maf_merge1)


pdf("MesKit/JSI.Met1.Paired.pdf", width = 8*1.1, height = 6*1.1)

#cowplot::plot_grid(plotlist = plist, nrow = 2)

wrap_plots(plist) + plot_layout(nrow = 1)

dev.off()


###################################################################################

#references: Quantitative evidence for early metastatic seeding in colorectal cancer
#calculate the H(Lm and Lp) for the early or layer divergence.
# calHindex
#' @param PrimaryId Primary tumor IDs to indicate the primary-metastases relationships.
#' @param CCF_cutoff The minimal cutoffs for the present status.


calHindex = function(ccf.list , JSI.res, maf_merge1, PrimaryId = "Coad", CCF_cutoff = 0.25 ){
  
  plist = list()
  
  sample.pairs = names(ccf.list)[grepl(PrimaryId , names(ccf.list))]
  
  stats = list()
  
  
  for(i in sample.pairs ){
    
    nums = c()
    
    message(i)
    
    pairs = as.character( str_split(i, "-", simplify = T))
    
    JSI = JSI.res$JSI.pair[pairs[1],pairs[2]]
    
    #merge driver info
    Met1_C_OL = ccf.list[[i]] %>%
      left_join( maf_merge1 %>% 
                   mutate(Mut_ID = str_c(Chromosome, Start_Position, Reference_Allele , Tumor_Seq_Allele2, sep = ":") ) %>%
                   dplyr::select(Mut_ID, is.driver) %>%
                   unique.data.frame()
      ) 
    
    Met1_group = Met1_C_OL %>%
      dplyr::select( all_of(c(PrimaryId, pairs[!pairs %in% PrimaryId ]))) %>%
      setNames(nm = c("Primary","Metastases")) %>%
      #Set Lp and Lm  
      mutate(
          mut_group = ifelse(Primary>=CCF_cutoff & Metastases<CCF_cutoff, "Lp",
                             ifelse(Metastases >= CCF_cutoff & Primary<CCF_cutoff, "Lm",
                                    ifelse(Primary>=CCF_cutoff & Metastases >= CCF_cutoff, "La", "Ls")
                                    ) )) %>%
      mutate(
        S1 = ifelse(mut_group == "Lp" &  Primary >=0.1, 1, 0 ),
        S2 = ifelse(mut_group == "Lp" &  Primary >=0.2, 1, 0 ),
        S3 = ifelse(mut_group == "Lp" &  Primary >=0.4, 1, 0 ),
        S4 = ifelse(mut_group == "Lp" &  Primary >=0.6, 1, 0 ),
        S5 = ifelse(mut_group == "Lm" &  Metastases >=0.1, 1, 0 ),
        S6 = ifelse(mut_group == "Lm" &  Metastases >=0.2, 1, 0 ),
        S7 = ifelse(mut_group == "Lm" &  Metastases >=0.4, 1, 0 ),
        S8 = ifelse(mut_group == "Lm" &  Metastases >=0.6, 1, 0 ),
        S9 = ifelse( Metastases >0.6 & Primary <0.6 &  Primary> CCF_cutoff , 1, 0 )
      )
    
     nums[1] = sum(Met1_group$mut_group == "Lp")
     nums[2] = sum(Met1_group$mut_group == "Lm")
     nums[3] = sum(Met1_group$mut_group == "La")
     nums[4] = sum(Met1_group$mut_group == "Ls")
     
     nums = setNames(nums, c("Lp","Lm","La","Ls") )
     
     nums = c(nums,
              Met1_group %>%
       dplyr::select(all_of( paste0("S",1:9) ) ) %>%
       colSums())
       
     stats[[i]] = nums
     
     
     plist[[i]] =  Met1_C_OL %>%
       ggplot(aes_string(x = pairs[1] , y = pairs[2]) ) + theme_classic2() +
       stat_density2d(aes(fill = ..density..^2), geom = "tile", contour = FALSE, n = 100) +
       scale_fill_continuous(low = "white", high = "red") +
       geom_hline(yintercept = CCF_cutoff, linetype = 2, size = 0.5) +
       geom_vline(xintercept = CCF_cutoff, linetype = 2, size = 0.5) +
       geom_point(alpha = 0.5) +
       #geom_point(col = "blue", data = Met1_C_OL %>% filter(is.driver), size = 3, shape = 5 ) +
       #geom_label_repel(aes(label = Hugo_Symbol), data = Met1_C_OL %>% filter(is.driver),
       #                  nudge_x = 0.1, nudge_y = 0.1
       #) +
       labs(title = sprintf("JSI = %.3f\nH = %.3f", JSI, nums[2]/nums[1])  ) +
       theme(
         legend.position = "none",
         plot.title = element_text(hjust = 0.5)
       )
     
     
  }
  
  return( list(plist = plist,
          stats = stats
          )
  )
  
}



Hindex =  calHindex(ccf.list = ccf.list$Met1, JSI.res,  maf_merge1, 
                    PrimaryId = "Coad",
                    CCF_cutoff = 0.2 )


pdf("JSI.Met1.Paired-2.pdf", width = 8*1.1, height = 3.5*1.1)

wrap_plots(Hindex$plist) + plot_layout(nrow = 1)

dev.off()


stats = Hindex$stats %>%
  purrr::reduce(rbind) %>%
  as.data.frame() %>%
  mutate(Tumor = names(Hindex$stats)) 


write.table(stats,
            file = "MesKit/CRC.summary_statistics.txt",
            row.names = F, quote = F, sep = "\t"
            )

# see https://github.com/cancersysbio/SCIMET for parameter estimation.


########################################################################################################

neutralResult <- testNeutral(maf1, min.mut.count = 10, patient.id = 'Met1', use.tumorSampleLabel = TRUE)

neutralResult$neutrality.metrics
neutralResult$model.fitting.plot$Coad_1
neutralResult <- testNeutral(maf, min.mut.count = 10, patient.id = 'V402', use.tumorSampleLabel = TRUE)


#9. Phylogenetic tree visualization
phyloTree <- getPhyloTree(maf, patient.id = "Met1", method = "NJ", min.vaf = 0.05)

# A phylogenetic tree along with binary and CCF heatmap of mutations 
phylotree_V402 <- plotPhyloTree(phyloTree, use.tumorSampleLabel = TRUE)
binary_heatmap_V402 <- mutHeatmap(maf, min.ccf = 0.04, use.ccf = FALSE, patient.id = "Breast", use.tumorSampleLabel = TRUE)
CCF_heatmap_V402 <- mutHeatmap(maf, use.ccf = TRUE, patient.id = "Breast", 
                               min.ccf = 0.04, use.tumorSampleLabel = TRUE)
cowplot::plot_grid(phylotree_V402,
                   CCF_heatmap_V402, nrow = 1, rel_widths = c(1.2, 1))


###################################################################

library(BSgenome.Hsapiens.UCSC.hg19)

phyloTree <- getPhyloTree(maf, patient.id = NULL , method = "ML", min.ccf = 0.05)

mutClass <- mutTrunkBranch(phyloTree, CT = TRUE, plot = TRUE)
names(mutClass)

mutClass$mutTrunkBranch.res

mutClass$mutTrunkBranch.plot

trimatrix_V402 <- triMatrix(phyloTree, level = 1)

#plotMutSigProfile(trimatrix_V402)[[1]]

fit_V402 <- fitSignatures(trimatrix_V402, signaturesRef = "cosmic_v2", signature.cutoff = 0.2, 
                          associated = paste("Signature", c(1,6,14,15,29))  
                          )

#plotMutSigProfile(fit_V402)[[1]]

ComplexHeatmap::Heatmap( fit_V402$OveryLM$cosine.similarity, name = "Cosine similarity")


ComplexHeatmap::Heatmap( fit_V402$OveryLM$cosine.similarity, name = "Cosine similarity")



###############################################################################################
#see total mutations across different sites.

phyloTree <- getPhyloTree(maf, patient.id = NULL , method = "ML", min.ccf = 0.05)

trimatrix_V402 <- triMatrix(phyloTree, level = 1)

#fit_V402 <- fitSignatures(trimatrix_V402, signaturesRef = "cosmic_v2", signature.cutoff = 0.1, 
#                          associated = paste("Signature", c(1,6,15,4,18,24,14,29,11))  
#)

fit_V402 <- fitSignatures(trimatrix_V402, signaturesRef = "cosmic_v2", signature.cutoff = 0.2, 
                          #associated = paste("Signature", c(1,6,15,18,24, 4))  
)

#combined all samples.
col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue4", "white", "red3"))

pdf("MesKit/Signatures.pdf", width = 6, height = 4)

#lapply(fit_V402, function(x) x$cosine.similarity ) %>% purrr::reduce(rbind) %>%
#  ComplexHeatmap::Heatmap( name = "Cosine similarity", col = col_fun)


  sigContri = lapply(fit_V402, function(x) x$contribution.relative ) %>% purrr::reduce(rbind) 
  
  table(sigContri %>% colSums() >0.1)
  #FALSE  TRUE 
  #24     6  
  select.sig = colnames(sigContri)[sigContri %>% colSums() >0.1]
  
  ComplexHeatmap::Heatmap( 
    sigContri[, select.sig], 
    name = "Relative\nContribution", col = col_fun)

dev.off()
