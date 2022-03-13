
library(devtools)

use_gpl_license(version = 3, include_future = TRUE)

library(usethis)



#add packages.
use_package("tidyverse", type = "depends")


document()
load_all()


library(ggtree)
library(phangorn)
library(treeio)

library(MPTevol)


#For CNAs trees.


dist = system.file(package="MPTevol", "extdata", "tree_final.dist")

dist = "tmp/analysis/medicc/results/Pan/tree_final.dist"

grp <- list(NORMAL = "NORMAL",
            Breast =  paste0("Breast_", 1:5),
            Coad = paste0("Coad_", 1:5),
            Lung  = paste0("Lung_", 1:5),
            OveryLM = paste0("OveryLM_", 1:5),
            OveryRM = paste0("OveryRM_", 1:6),
            UterusM = paste0("UterusM_", c(1:7))
)


plotCNAtree(dist = dist,
            grp = grp
            )



#add vignette
usethis::use_vignette("plotTrees", title = "plotTrees - add p-value in the boxplot")


########################################################################################


#build a data structure.

library(MesKit)

#build a maf data: see MesKit

data.type = "split1"

maf <- readMaf(mafFile = sprintf( "tmp/analysis/CCF_estimated_CHAT/MesKit/meskit.%s.mutation.txt", data.type),
               ccfFile = sprintf("tmp/analysis/CCF_estimated_CHAT/MesKit/meskit.%s.CCF.txt", data.type),
               clinicalFile  = sprintf("tmp/analysis/CCF_estimated_CHAT/MesKit/meskit.%s.clinical.txt", data.type),
               refBuild = "hg19",
               ccf.conf.level = 0.90
)


driverGene = read.delim("tmp/analysis/CCF_estimated_CHAT/data/IntOGen-Drivers-Cancer_Genes.tsv", header = T) %>%
  filter(CANCER_TYPE %in% c("BRCA","COREAD","LUAD", "LUSC") ) %>%
  pull(SYMBOL) %>% unique()

mut.class <- classifyMut(maf, class =  "SP", patient.id = 'Breast')

plotMutProfile(maf, class =  "SPCS", geneList = driverGene, use.tumorSampleLabel = TRUE, removeEmptyCols = FALSE)


#get data
maf_input = subMaf(maf, mafObj = FALSE, use.tumorSampleLabel = TRUE)
#get mutation classifications.
maf_class = classifyMut(maf, patient.id = NULL, class = "SPCS", classByTumor = FALSE)



################################################################

#see clinical targetable sites


###################################################################


readSegment



library(CNAqc)

plotCNAProfile

cnv.inputA = "/data1/qingjian/Project/ThreePrimary/Results/sequenza/Seg/Breast_1_segments.txt"

cnvA = read.delim(cnv.inputA, header = T, sep = "\t")

cnaqc.x = readCNAProfile(
  TumorName = "Breast_1",
  maf = maf$Breast,
  seg = cnvA,
  purity = 0.9
)






