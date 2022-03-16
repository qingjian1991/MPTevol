
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


######################################################################


phyloTree <- getPhyloTree(maf, patient.id = "Breast", method = "NJ", min.vaf = 0.02,
                          bootstrap.rep.num = 1000)



#######################################################################

# Maf object with CCF information
data.type = "split1"
maf <- readMaf(mafFile = sprintf( "tmp/analysis/CCF_estimated_CHAT/MesKit/meskit.%s.mutation.txt", data.type),
               ccfFile = sprintf("tmp/analysis/CCF_estimated_CHAT/MesKit/meskit.%s.CCF.txt", data.type),
               clinicalFile  = sprintf("tmp/analysis/CCF_estimated_CHAT/MesKit/meskit.%s.clinical.txt", data.type),
               refBuild = "hg19",
               ccf.conf.level = 0.95
)


cal = calRoutines(maf = maf,
          patient.id = "Met1",
          PrimaryId = "Coad",
          pairByTumor = TRUE,
          use.tumorSampleLabel = TRUE,
          subtitle = "both"
)


plot_grid(plotlist = cal$Met1$plist, nrow = 1)


## add driver information

maf_driver = data.frame(
  Mut_ID = c("5:112170777:CAGA:-","1:147092680:-:C"),
  is.driver = c(TRUE,TRUE)
)


cal = calRoutines(maf = maf,
                patient.id = "Met1",
                PrimaryId = "Coad",
                pairByTumor = TRUE,
                use.tumorSampleLabel = TRUE,
                subtitle = "both",
                maf_drivers = maf_driver
)

plot_grid(plotlist = cal$Met1$plist, nrow = 1)


#############################################################################################

#build clonal evolution tree

mutdata.vaf1 = readRDS("tmp/analysis/CCF_estimated_CHAT/data/mutdata.vaf1.rds")
new.label = readRDS(file = "tmp/analysis/CCF_estimated_CHAT/data/new.label.rds")

mutdata = mutdata.vaf1 %>%
  left_join(new.label)

mutdata[, 2:34] = mutdata[, 2:34]*100

mutdata$cluster = mutdata$Clone2

#Visualizing the variant clusters
set.colors = c("#C6C6C6", "#6FDCBF","#5AA36E","#E99692","#B4D985","#EA7D23","#E53CFB","#4B85B7","#8439FA","#BD8736","#B3B371","#A7C7DE","#EE97FC","#57C222","#BFABD0","#44589B", "#794C18", RColorBrewer::brewer.pal(n = 10, name = "Paired") )

clone.colors = set.colors[ 1:length(unique(mutdata$cluster))]
cluster.col.name  = "cluster"

head(mutdata)



GC01 = readRDS("tmp/analysis/CCF_estimated_CHAT/data/GC01.rds")

GC01 = GC01 %>% mutate(Chromosome =  str_c("chr", Chromosome) ) %>%
  mutate(mutid = str_c(Chromosome, Start_Position, Reference_Allele, sep = ":"))


drivers = read.delim("tmp/analysis/CCF_estimated_CHAT/data/IntOGen-Drivers-Cancer_Genes.tsv", header = T) %>%
  filter(CANCER_TYPE %in% c("BRCA","COREAD","LUAD", "LUSC") ) %>%
  pull(SYMBOL) %>% unique()

mutinfo = GC01 %>%
  select(Chromosome, Start_Position, End_Position, Hugo_Symbol, Variant_Classification, Protein_Change, Reference_Allele, Tumor_Seq_Allele2 ) %>%
  unique.data.frame() %>%
  mutate(#Chromosome = str_c("chr", Chromosome),
    Protein_Change = ifelse(Variant_Classification == "Splice_Site", "splice", Protein_Change),
    Protein_Change = ifelse(Variant_Classification == "3'UTR", "3UTR", Protein_Change),
    Protein_Change = ifelse(Variant_Classification == "5'UTR", "5UTR", Protein_Change),
    gene_site = str_c(Hugo_Symbol, "_" ,Protein_Change),
    is.driver = ifelse( (Hugo_Symbol %in% drivers) & ( !Variant_Classification %in% c("Silent", "3'Flank", "IGR", "Intron", "RNA") ) , TRUE, FALSE)  ) %>%

  mutate( mutid  = str_c(Chromosome , Start_Position , Reference_Allele , sep = ":") )


mutdata = right_join(mutinfo, mutdata ) %>%
  arrange(cluster) %>%
  as.data.frame()

#add ref and var columns for merge data.

sampleNames = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                 paste0("Coad_", 1:5), paste0("OveryLM_", 1:5),
                 paste0("OveryRM_", 1:6), paste0("UterusM_", c(1:7)) )

var.data = round(mutdata[, sampleNames])
ref.data = 100 - var.data

colnames(var.data) = paste0(colnames(var.data), ".var")
colnames(ref.data) = paste0(colnames(ref.data), ".ref")

mutdata = cbind(mutdata, var.data, ref.data)

message("<2> Analysis using clonevol")


output = "Met2"



pp <- plot.variant.clusters(mutdata,
                            show.cluster.size = F,
                            show.cluster.label= F,
                            cluster.col.name = cluster.col.name,
                            vaf.col.names = vaf.col.names,
                            violin = FALSE,
                            box = TRUE,
                            jitter = TRUE,
                            jitter.shape = 1,
                            variant.class.col.name =cluster.col.name,
                            jitter.color = clone.colors,
                            jitter.size = 1.2,
                            jitter.alpha = 1,
                            jitter.width = 0.2,
                            jitter.center.method = "median",
                            jitter.center.size = 1,
                            jitter.center.color = "darkgray",
                            jitter.center.display.value = "none",
                            display.plot = F,
                            horizontal = T,
                            order.by.total.vaf = F,
                            highlight = 'is.driver',
                            highlight.shape = 21,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.size = 2.5,
                            highlight.note.col.name = NULL,
                            highlight.note.size = 2,
                            highlight.note.color = "blue",
                            highlight.note.angle = 0,
                            founding.cluster = 1,
                            ccf = show.ccf
)


pdf( sprintf( "%s/%s.box.pdf", output, output), width = 2*length(pp), height = 4)
ggpubr::ggarrange(plotlist = pp, ncol = length(pp), align  = "h" )
dev.off()


output = "Met2"

saveRDS(mutinfo, file = sprintf("%s/%s.mutinfo.rds", output))

sampleNames = c( paste0("Coad_", 1:5), paste0("OveryLM_", 1:5),
                 paste0("OveryRM_", 1:6),paste0("UterusM_", c(1, 3)) )


sample.groups = mapply(function(x) x[1], strsplit(vaf.col.names, "_")  )
names(sample.groups) = vaf.col.names

vaf.col.names = sampleNames

show.ccf = F

sel = c(1,6,11)

y = inferClonalTrees( project.names  = output ,
                      variants = mutdata,
               ccf.col.names = vaf.col.names[sel],
               sample.groups = sample.groups[sel],
               cancer.initiation.model = "monoclonal",
               founding.cluster = 1 , ignore.clusters = 4 , cluster.col.name = "cluster",
               subclonal.test.model = "non-parametric",
               sum.p = 0.01, alpha = 0.05, weighted = FALSE,
               plot.pairwise.CCF  = FALSE,
               highlight.note.col.name = NULL,
               highlight = "is.driver", highlight.CCF = TRUE
)

pdf(file = sprintf("%s/%s.trees.pdf", output, output), width = 6, height = 6)

plot.all.trees.clone.as.branch(y, branch.width = 0.5,
                               node.size = 2, node.label.size = 0.5,
                               tree.rotation	=180,
                               angle = 20
)

dev.off()


samples = names(y$models)

times = clonevol2timescape(results = y, samples = names(y$models)  )


# run timescape
i = 1
timescape(clonal_prev = times$clonal_prev[[i]],
          tree_edges = times$tree_edges[[i]],
          clone_colours = times$clone_colours[[i]],
          genotype_position = "centre",
          xaxis_title = NULL
)











pp = plotVafCluster(
  variants = variants,
  cluster.col.name = "cluster",
  vaf.col.names = vaf.col.names[1:6],
  highlight = "is.driver",
  highlight.note.col.name = "gene_site",
  box = F,
  violin = TRUE
)

pp = plotVafCluster(
  variants = variants,
  cluster.col.name = "cluster",
  vaf.col.names = vaf.col.names[1:6],
  #highlight = "is.driver",
  #highlight.note.col.name = "gene_site",
  box = F,
  violin = TRUE
)








