
library(devtools)

use_gpl_license(version = 3, include_future = TRUE)

library(usethis)


add_

#add packages.
use_package("tidyverse", type = "depends")

use_package("clonevol", type = "Imports" )

document()
load_all()

use_vignette("MPTevol")


library(ggtree)
library(phangorn)
library(treeio)

library(MPTevol)


#read samples distances.
#This dist file is the output of MEDICC
dist = system.file(package="MPTevol", "extdata", "tree_final.dist")

#set group information
group <- list(NORMAL = "NORMAL",
            Breast =  paste0("Breast_", 1:5),
            Coad = paste0("Coad_", 1:5),
            Lung  = paste0("Lung_", 1:5),
            OveryLM = paste0("OveryLM_", 1:5),
            OveryRM = paste0("OveryRM_", 1:6),
            UterusM = paste0("UterusM_", c(1:7))
)

#set group colors
group.colors = setNames( set.colors(n = length(group)), nm = names(group) )

#built trees
tree = plotCNAtree(dist = dist,
            group = group,
            group.colors = group.colors
            )

tree$plot


viewTrees(phyloTree = tree$phyloTree,
          tree.format = "list",
          group = group,
          group.colors = group.colors
          )




#add vignette
usethis::use_vignette("plotTrees", title = "plotTrees - add p-value in the boxplot")


########################################################################################


#build a data structure.

library(MesKit)

#build a maf data: see MesKit

data.type = "split1"

maf <- readMaf(mafFile = system.file(package="MPTevol", "extdata", sprintf("meskit.%s.mutation.txt", data.type)),
               ccfFile = system.file(package="MPTevol", "extdata", sprintf("meskit.%s.CCF.txt", data.type)),
               clinicalFile = system.file(package="MPTevol", "extdata", sprintf("meskit.%s.clinical.txt", data.type)),
               refBuild = "hg19",
               ccf.conf.level = 0.95
)


driverGene = read.delim(system.file(package="MPTevol", "extdata", "IntOGen-Drivers-Cancer_Genes.tsv"), header = T) %>%
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

sites = getClinSites(maf, Patient_ID = "Breast")

DT::datatable(sites)

sites = getClinSites(maf)


###################################################################

#calKaKas

kaks = calKaKs(maf, patient.id = "Breast", class = "SP", parallel = TRUE, vaf_cutoff = 0.05)

kaks = calKaKs(maf, patient.id = "Breast", class = "CS", parallel = TRUE, vaf_cutoff = 0.05)

kaks = calKaKs(maf, class = "SP", parallel = TRUE, vaf_cutoff = 0.05)


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


phyloTree <- getPhyloTree(maf, patient.id = "Met1", method = "NJ", min.vaf = 0.02,
                          bootstrap.rep.num = 1000)

viewTrees(phyloTree)

group <- list(
              Coad = paste0("Coad_", 1:5),
              OveryLM = paste0("OveryLM_", 1:5),
              OveryRM = paste0("OveryRM_", 1:6),
              UterusM = paste0("UterusM_", c(1,3))
)

viewTrees(phyloTree = phyloTree,
          group = group
          )




phyloTree <- getPhyloTree(maf, patient.id = "Met1", method = "NJ", min.vaf = 0.02,
                          bootstrap.rep.num = 1000)


plotMutTree(maf, patient.id = "Met1")
mutTrees = plotMutTree(maf, patient.id = "Met1", group = group, title = "CRC Met")
mutTrees1 = plotMutTree(maf, patient.id = "Met1", group = group, title = "CRC Met", method = "MP")


mutTrees$plot

compareTree(mutTrees$phyloTree,
            mutTrees1$phyloTree,
            plot = TRUE
            )





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



#######################################################################################################


#save(variants, file = "variants.rda")
#save(variants.ref, file = "variants.ref.rda")


# load data
data("variants", package = "MPTevol")
data("variants.ref", package = "MPTevol")

vaf.col.names = c( paste0("Coad_", 1:5), paste0("OveryLM_", 1:5),
                 paste0("OveryRM_", 1:6),paste0("UterusM_", c(1, 3)) )

sample.groups = mapply(function(x) x[1], strsplit(vaf.col.names, "_")  )
names(sample.groups) = vaf.col.names

cluster.col.name = "cluster"

clones.number = 10
clone.colors = set.colors(10)

#Check data.
pp = plotVafCluster(
  variants = variants,
  cluster.col.name = "cluster",
  vaf.col.names = vaf.col.names[ c(1,5,6, 10, 11, 16)],
  highlight = "is.driver",
  highlight.note.col.name = "gene_site",
  box = TRUE,
  violin = FALSE
)

pp

#check cluster changes.
plot.cluster.flow(variants,
                  cluster.col.name = cluster.col.name,
                  vaf.col.names = vaf.col.names,
                  sample.names = vaf.col.names,
                  colors = set.colors(clones.number),
                  y.title = "Variant Allele Frequency %"
) +
  theme(axis.text.x = element_text(angle = 90) )


###### inferring subclonal structures.

#sel = c(1,5,6, 10, 11, 16)
sel = 1:18

y = inferClonalTrees( project.names  = "Met" ,
                      variants = variants.ref,
                      ccf.col.names = vaf.col.names[sel],
                      sample.groups = sample.groups[sel],
                      cancer.initiation.model = "monoclonal",
                      founding.cluster = 1 , ignore.clusters = 4 ,
                      cluster.col.name = "cluster",
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


# We combined MRS samples into one samples.

#Merge multiple samples into a single sample. This can be used to merge multi-region samples to one sample representing the tumor where the regions are taken from. This functions required the ref and var columns. So, we used the variants.ref

#Note: CCF should be in the range between 0-100

#merge coad
sel = 1:5
y.merge = merge.samples(y, samples = vaf.col.names[sel],
                        new.sample = "Coad", new.sample.group = "Coad",
                        ref.cols = str_c(vaf.col.names[sel], ".ref"),
                        var.cols = str_c(vaf.col.names[sel], ".var") )
#merge overyLM
sel = 6:10
y.merge = merge.samples(y.merge, samples = vaf.col.names[sel],
                        new.sample = "OveryLM", new.sample.group = "OveryLM",
                        ref.cols = str_c(vaf.col.names[sel], ".ref"),
                        var.cols = str_c(vaf.col.names[sel], ".var") )
#merge overyRM
sel = c(11:13, 14, 15, 16)
y.merge = merge.samples(y.merge, samples = vaf.col.names[sel],
                        new.sample = "OveryRM", new.sample.group = "OveryRM",
                        ref.cols = str_c(vaf.col.names[sel], ".ref"),
                        var.cols = str_c(vaf.col.names[sel], ".var") )

#merge UterusM
sel = c(17, 18)
y.merge = merge.samples(y.merge, samples = vaf.col.names[sel],
                        new.sample = "UterusM", new.sample.group = "UterusM",
                        ref.cols = str_c(vaf.col.names[sel], ".ref"),
                        var.cols = str_c(vaf.col.names[sel], ".var") )

#####################################################################################
#View clonal evolving.



samples = names(y.merge$models)

times = tree2timescape(results = y.merge, samples = names(y.merge$models)  )


# run timescape
i = 1
timescape(clonal_prev = times$clonal_prev[[i]],
          tree_edges = times$tree_edges[[i]],
          clone_colours = times$clone_colours[[i]],
          genotype_position = "centre",
          xaxis_title = NULL
)

#save svg file plot
rsvg::rsvg_pdf( svg = sprintf("%s/%s.model%s.svg", project.names, project.names, i),
                file = sprintf("%s/%s.model%s.pdf", project.names, project.names, i)
)

#######################################################################################################

#Running MEDICC

folder = "/data1/qingjian/Rproject/Three/medicc/Seg.new1"
segfiles = list.files( folder, pattern = ".txt", full.names = T)
sampleid = list.files( folder, pattern = ".txt") %>% stringr::str_remove("_segments.txt")


#Running for Breast

library(IRanges)

seg = splitSegment(
  segfiles = segfiles[1:5],
  sampleid = sampleid[1:5],
  project.names = "Breast",
  out.dir = "medicc/Breast",
  N.baf = 30,
  cnv_min_length = 1e+05,
  max_CNt = 15,
  minLength = 1e+05,
  maxCNV = 4
)


#Running for Met

library(IRanges)

seg = splitSegment(
  segfiles = segfiles[c(6:10,16:27, 29)],
  sampleid = sampleid[c(6:10,16:27, 29)],
  project.names = "Met",
  out.dir = "medicc/Met",
  N.baf = 30,
  cnv_min_length = 1e+05,
  max_CNt = 15,
  minLength = 1e+05,
  maxCNV = 4,
  medicc.py = NULL,
  python = "python"
)

## Running Medicc.

#medicc=/data1/soft/medicc/medicc.py
#python=/data1/qingjian/anaconda3/envs/pyclone/bin/python
#nohup $python $medicc Breast/Breast.descr.txt results/Breast.out -v > results/Breast.runinfo.txt &


plotMediccSeg = function(file){
  major = read.table(file, header = T)

  major = major %>%
    mutate(seq = str_c(seqnames, start, end, sep = "_")) %>%
    column_to_rownames(var = "seq") %>%
    mutate(seqnames = NULL, start = NULL, end = NULL, width = NULL)

  Heatmap(major,
          row_names_gp = gpar(fontsize = 6)
  )

}




plotMediccSeg("medicc/Breast/Breast.major.txt")
plotMediccSeg("medicc/Breast/Breast.minor.txt")

plotMediccSeg("medicc/Met/Met.major.txt")
plotMediccSeg("medicc/Met/Met.minor.txt")


#Plot CNA-based trees from the MEDICC.


#read samples distances.
#This dist file is the output of MEDICC
dist = "medicc/Breast.run/pairwise_summed.dist"

#set group information
group <- list(
              Breast =  paste0("Breast_", 1:5)
)

#set group colors
#group.colors = setNames( set.colors(n = length(group), rev = T), nm = names(group) )
group.colors = setNames( "#9F994E", nm = names(group) )

#built trees
tree = plotCNAtree(dist = dist,
                   group = group,
                   group.colors = group.colors
)

tree$plot

pdf("medicc/Breast.CNAs.Trees.pdf", width = 5.5, height = 4.5)
tree$plot
dev.off()


#Plot CNA-based trees from the MEDICC.


#read samples distances.
#This dist file is the output of MEDICC
dist = "medicc/Met.run/tree_final.dist"

#set group information
#set group information
group <- list(
              Coad = paste0("Coad_", 1:5),
              OveryLM = paste0("OveryLM_", 1:5),
              OveryRM = paste0("OveryRM_", 1:6),
              UterusM = paste0("UterusM_", c(1,3))
)

#set group colors
group.colors = setNames( c("#7570B3","#E6AB02","#003C30","#666666") , nm = names(group) )

#built trees
tree = plotCNAtree(dist = dist,
                   group = group,
                   group.colors = group.colors
)


pdf("medicc/Met.CNAs.Trees.pdf", width = 9.7, height = 7.5)
tree$plot
dev.off()


x = tree$phyloTree$tree %>% as.tibble()


tidytree::child(x, 4)

viewTrees

######################################################################

dist = "../medicc/results-new/pan1/pairwise_summed.dist"

#built trees
tree = plotCNAtree(dist = dist,
                   group = NULL,
                   group.colors = NULL
)

# we want to cut the trees.

# we plot the node ids in the trees.

ggtree(tree$phyloTree$tree, size=1) +
  geom_tiplab(size=4) +
  geom_nodelab(aes(label = round(node)) )

tree.info = tree$phyloTree$tree %>%
  as_tibble()


KaKs_data = kaks$Breast$KaKs_data

KaKs_data %>%
  filter(name %in% c("wall")) %>%
  mutate(name = factor(name, levels = c("wall") )) %>%
  #mutate(Type = factor(Type, levels = c("Shared_Clonal","Private_Clonal","Private_Subclonal") ))  %>%
  ggplot(aes(x = Tumor_ID , y = mle, fill = Type) ) + ggpubr::theme_pubr() +
  geom_bar(stat = "identity", position =  position_dodge(width = 0.90)) +
  #geom_linerange(aes(ymin = cilow, ymax = cihigh), position =  position_dodge(width = 0.90) ) +
  geom_hline( yintercept = 1, linetype = 2, size = 1) +
  labs(x = NULL, y =  latex2exp::TeX("Dn/Ds ($\\omega_{all}$)") )  +
  scale_fill_manual(values = set.colors(length(unique(KaKs_data$Type))) ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14)
  )



#rename the IDs.

samples.IDs = data.frame(
  Tumor_Sample_Barcode = c(
  paste0("Coad_", 1:5), paste0("OveryLM_", 1:5),
  paste0("OveryRM_", 1:6), paste0("UterusM_", 1:7),
  paste0("Breast_", 1:5), paste0("Lung_", 1:5)
  ),

  Tumor_Sample_Barcode1 = c(
    paste0("READ_", 1:5), paste0("OvaryLM_", 1:5),
    paste0("OvaryRM_", 1:6), paste0("UterusM_", 1:7),
    paste0("BRCA_", 1:5), paste0("LNET_", 1:5)
  ),

  Patient_ID1 = c(
    rep("READ", 5), rep("OvaryLM", 5),
    rep("OvaryRM", 6), rep("UterusM", 7),
    rep("BRCA", 5), rep("LNET", 5)
  )

)


data.type = "split"


#(1) CCF
ccf = read.table( sprintf( "inst/extdata/meskit.%s.CCF.txt", data.type), header = T)

head(ccf)

ccf = ccf %>% left_join(samples.IDs) %>%
  mutate(Tumor_Sample_Barcode = Tumor_Sample_Barcode1,
         Patient_ID = Patient_ID1,
         Tumor_Sample_Barcode1 = NULL,
         Patient_ID1 = NULL
         )

write.table(ccf, file = sprintf( "inst/extdata/meskit.%s.CCF.txt", data.type), row.names = F, sep = "\t", quote = F)


#(2) mutations.

mutation = read.delim( sprintf( "inst/extdata/meskit.%s.mutation.txt", data.type), header = T, sep = "\t")

mutation = mutation %>%
  left_join(samples.IDs) %>%
  mutate(
    Tumor_Sample_Barcode = Tumor_Sample_Barcode1,
    Tumor_Sample_Barcode1 = NULL,
    Patient_ID1 = NULL
  )

write.table(mutation, file = sprintf( "inst/extdata/meskit.%s.mutation.txt", data.type), row.names = F, sep = "\t", quote = F)

#(3) clinical
data.type = "split"

clinical = read.delim( sprintf( "inst/extdata/meskit.%s.clinical.txt", data.type), header = T, sep = "\t")

clinical = clinical %>%
  left_join(samples.IDs) %>%
  mutate(
    Tumor_Sample_Barcode = Tumor_Sample_Barcode1,
    Tumor_Sample_Barcode1 = NULL,
    Patient_ID1 = NULL
  )

write.table(clinical, file = sprintf( "inst/extdata/meskit.%s.clinical.txt", data.type), row.names = F, sep = "\t", quote = F)




#rename the IDs.

samples.IDs = data.frame(
  Tumor_Sample_Barcode = c(
    paste0("Coad_", 1:5), paste0("OveryLM_", 1:5),
    paste0("OveryRM_", 1:6), paste0("UterusM_", 1:7),
    paste0("Breast_", 1:5), paste0("Lung_", 1:5)
  ),

  Tumor_Sample_Barcode1 = c(
    paste0("READ_", 1:5), paste0("OvaryLM_", 1:5),
    paste0("OvaryRM_", 1:6), paste0("UterusM_", 1:7),
    paste0("BRCA_", 1:5), paste0("LNET_", 1:5)
  ),

  Patient_ID1 = c(
    rep("Met1", 5), rep("Met1", 5),
    rep("Met1", 6), "Met1", "Uterus" ,"Met1", rep("Uterus", 4),
    rep("BRCA", 5), rep("LNET", 5)
  )

)


data.type = "split1"

#(1) CCF
ccf = read.table( sprintf( "inst/extdata/meskit.%s.CCF.txt", data.type), header = T)

head(ccf)

ccf = ccf %>% left_join(samples.IDs) %>%
  mutate(Tumor_Sample_Barcode = Tumor_Sample_Barcode1,
         Patient_ID = Patient_ID1,
         Tumor_Sample_Barcode1 = NULL,
         Patient_ID1 = NULL
  )

write.table(ccf, file = sprintf( "inst/extdata/meskit.%s.CCF.txt", data.type), row.names = F, sep = "\t", quote = F)


data("variants", package = "MPTevol")
data("variants.ref", package = "MPTevol")


clumns = colnames(variants)

clumns = str_replace_all(clumns, "Breast", "BRCA")
clumns = str_replace_all(clumns, "Lung", "LNET")
clumns = str_replace_all(clumns, "Coad", "READ")
clumns = str_replace_all(clumns, "Overy", "Ovary")

colnames(variants) = clumns


clumns = colnames(variants.ref)

clumns = str_replace_all(clumns, "Breast", "BRCA")
clumns = str_replace_all(clumns, "Lung", "LNET")
clumns = str_replace_all(clumns, "Coad", "READ")
clumns = str_replace_all(clumns, "Overy", "Ovary")

colnames(variants.ref) = clumns

save(variants, file = "data/variants.rda")
save(variants.ref, file = "data/variants.ref.rda")



plot.all.trees.clone.as.branch()


MPTevol::plot.all.trees.clone.as.branch()




clonevol::plot.all.trees.clone.as.branch()



clonevol:::plotBellsCells(y, out.pdf.file = "test.pdf")

y$models$READ_1

segs = read.table("inst/extdata/meskit.sequenza.CNAs.txt", header = T, sep = "\t")

head(segs)


segs = segs %>%
  filter(
  ! Tumor_Sample_Barcode %in% c( paste0("LNET_", 1:5), paste0("UterusM_", c(2,4:7) ) )
)

write.table(segs, file = "inst/extdata/meskit.sequenza.CNAs.txt", sep = "\t", quote = F, row.names = F)



















