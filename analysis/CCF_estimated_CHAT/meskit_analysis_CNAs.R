
library(MesKit)
library(tidyverse)

vardir = "MesKit//"

#for segmentation files from titanCNA.
read.seg = function(segfile , num.snps.cutoff = 100){
  
  seg = read.delim(file = segfile, header = T) %>%
    filter(Length.snp. > num.snps.cutoff ) %>%
    mutate(Chromosome  = str_replace(Chromosome , "chr","")) %>%
    mutate(Patient_ID = str_split(Sample, pattern = "_") %>% mapply(function(x) x[1], .) ) %>%
    select(Patient_ID, Sample, Chromosome, Start_Position.bp., End_Position.bp., Length.snp., Corrected_Copy_Number, Corrected_logR, logR_Copy_Number,  )
  
  colnames(seg) = c("Patient_ID","Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","num.snps","CopyNumber", "Log_Segment_Mean","Copy_Number_Abs")
  
  seg
  
}

segfile = list.files("TitanCNA/", pattern = "seg", full.names = T)

segCNA = c()

for(i in segfile){
  segCNA = rbind(segCNA,  read.seg(i))
}

write.table(segCNA, file = paste0(vardir, "/Three.TitanCNA.segCNA.txt"), quote = F, row.names = F, sep = "\t")


#for segmentation files from sequenza.

#***************************
#Filter for CNVs
#(1) Have at least 100 heterozygous loci within each segment;
#(2) Are at least 5Mb long;
#(3) Don't have extreme estimated copy numbers (e.g., total copy number >= 6);

read.seg = function(segfile , num.snps.cutoff = 100, seg.length =  5*10e5,
                    Tumor_Sample_Barcode = "tumor"){
  
  seg = read.delim(file = segfile, header = T) %>%
    mutate(length =  end.pos - start.pos) %>%
    filter(N.BAF > num.snps.cutoff ) %>%
    filter(length >=  seg.length) %>%
    mutate(chromosome  = str_replace(chromosome , "chr","")) %>%
    mutate(Tumor_Sample_Barcode = Tumor_Sample_Barcode) %>%
    mutate(Patient_ID =  str_split(Tumor_Sample_Barcode, pattern = "_") %>% mapply(function(x) x[1], .)  ) %>%
    select(Patient_ID, Tumor_Sample_Barcode, chromosome, start.pos, end.pos, N.BAF, CNt, A, B)
  
  colnames(seg) = c("Patient_ID","Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","num.snps","CopyNumber", "CopyNumber_A","CopyNumber_B")
  
  seg
  
}


segfile = list.files("/data1/qingjian/Project/ThreePrimary/Results/sequenza/Seg/", pattern = "_segments.txt", full.names = T)

sampelid = list.files("/data1/qingjian/Project/ThreePrimary/Results/sequenza/Seg/", pattern = "_segments.txt", full.names = F) %>% 
  str_replace(pattern = "_segments.txt", replacement = "") 

segCNA = c()

for(i in 1:length(segfile)){
  segCNA = rbind(segCNA,  read.seg(segfile[i], Tumor_Sample_Barcode = sampelid[i]) )
}

write.table(segCNA, file = paste0(vardir, "/Three.sequenza.segCNA.txt"), quote = F, row.names = F, sep = "\t")


# 5.2 segFiles --------------------------------------------------



Tumors =c("Breast", "Lung","UterusM","Coad","OveryLM","OveryRM")
# Read segment file
segCN <- paste0(vardir, "/Three.TitanCNA.segCNA.txt")
segCN.sequenza <- paste0(vardir, "/Three.sequenza.segCNA.txt")


# Reading gistic output files
all.lesions <- "MesKit/PANC_all_lesions.conf_99.txt"
amp.genes <- "MesKit/PANC_amp_genes.conf_99.txt"
del.genes <- "MesKit/PANC_del_genes.conf_99.txt"

#seg <- readSegment(segFile  = segCN, gisticAllLesionsFile = all.lesions,
#                   gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes)

#seg <- readSegment(segFile  = segCN)

#seg$Lung = seg.tmp$Lung 
#plotCNA(seg, patient.id = Tumors )


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#for Sequenza
seg.sequenza <- readSegment(segFile  = segCN.sequenza,
                            gisticAllLesionsFile = all.lesions,
                            gisticAmpGenesFile = amp.genes, 
                            gisticDelGenesFile = del.genes,
                            txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                            gistic.qval = 1e-8
                            )

#seg.sequenza <- readSegment(segFile  = segCN.sequenza
#)

pdf("MesKit/MutCNAsProfile.pdf", width = 10, height = 8)
plotCNA(seg.sequenza, patient.id = Tumors, chrSilent = c("X", "Y"),
        #showGene = TRUE,
        #showCytoband = TRUE
        )

dev.off()

############################################################################################

#### Using HeatmapAnnotation to plot

# set certain colors
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

library(ComplexHeatmap)

arms = read.delim("Gistics/broad_values_by_arm.txt", header = T, sep = "\t")
rownames(arms) = arms$Chromosome.Arm
arms$Chromosome.Arm = NULL

arms = t(arms)

#annotation for clinical information.
Patient = factor( mapply(function(x) x[1], strsplit(rownames(arms), "_") ) )

ra = rowAnnotation(
  Sites = Patient,
  col = list(
    Sites = setNames( colorScale[31:36], 
                      nm =  unique( Patient ))
  ),
  #show_legend = c("Sites" = FALSE),
  
  #annotation_legend_param = list( grid_width = unit(3, "mm") ),
  
  #control width of anno
  simple_anno_size = unit(2, "mm"),
  #show the title of anno
  show_annotation_name = FALSE
) 

#column annotation for colnames.http://192.168.120.51:8787/graphics/plot_zoom_png?width=1083&height=816

ba = columnAnnotation (foo = anno_mark(at = grep("q", colnames(arms) ) , 
                                   labels = c(seq(1:22), "X"),
                                   side = "column"
                                   )
                       )

mat = Heatmap( arms,
         cluster_rows = FALSE, cluster_columns = FALSE,
         right_annotation  = ha,
         bottom_annotation = ba,
         row_split = Patient,
         name = "mat",
         show_column_names = FALSE
         )


pdf("Gistics/Broad_copy_change_sequenza.pdf", width = 10, height = 8)

mat

for(i in grep("q", colnames(arms) )[1:22] ){
  for(j in 1:6){
    decorate_heatmap_body("mat",{
      i = i
      x = i/ncol(arms)
      grid.lines( c(x, x), c(0, 1), gp = gpar(lwd = 2, lty = 2, col = "grey70"))
    },
    row_slice = j
    )
  }
}

dev.off()


#add dash lines.
decorate_heatmap_body("mat",{
  i = c(2)
  x = i/ncol(arms)
  grid.lines( c(x, x), c(0, 1), gp = gpar(lwd = 2, lty = 2))
},
row_slice = 1
)

list_components()
