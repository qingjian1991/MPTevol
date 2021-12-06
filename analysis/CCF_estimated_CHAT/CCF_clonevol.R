
#Analysis CCF_estimated and using Clonevol to analysis the data.

library(clonevol)
library(ggrepel)

#Using CCF
#mutdata is in CCF_plot.

#> head(mutdata, 2)
#mutid Breast_1 Breast_2 Breast_3 Breast_4 Breast_5 Coad_1 Coad_2 Coad_3 Coad_4 Coad_5 Lung_1 Lung_2
#1 chr10:87898775:C        0        0        0     0.00        0      0      1      1   0.73   1.00      0   0.14
#2 chr12:56492633:A        0        0        0     0.19        0      0      1      1   1.00   0.84      0   0.00
#Lung_3 Lung_4 Lung_5 OveryLM_1 OveryLM_2 OveryLM_3 OveryLM_4 OveryLM_5 OveryRM_1 OveryRM_2 OveryRM_3 OveryRM_4
#1      0      0   0.31         1         1         1         1         1         1      0.38         1         1
#2      0      0   0.00         1         1         1         1         1         1      1.00         1         1
#OveryRM_5 OveryRM_6 UterusM_1 UterusM_2 UterusM_3 UterusM_4 UterusM_5 UterusM_6 UterusM_7 sciClone kmeans
#1      0.58      0.44      0.87      0.25      1.00         0         0         0         0        1      2
#2      1.00      1.00      0.99      0.00      0.45         0         0         0         0        1      2

mutdata.vaf1 = readRDS("data/mutdata.vaf1.rds")
new.label = readRDS(file = "data/new.label.rds")



mutdata = mutdata.vaf1 %>%
  left_join(new.label)

mutdata[, 2:34] = mutdata[, 2:34]*100

mutdata$cluster = mutdata$Clone2

#Visualizing the variant clusters
set.colors = c("#C6C6C6", "#6FDCBF","#5AA36E","#E99692","#B4D985","#EA7D23","#E53CFB","#4B85B7","#8439FA","#BD8736","#B3B371","#A7C7DE","#EE97FC","#57C222","#BFABD0","#44589B", "#794C18", RColorBrewer::brewer.pal(n = 10, name = "Paired") )

clone.colors = set.colors[ 1:length(unique(mutdata$cluster))]
cluster.col.name  = "cluster"


##########################################################Driver mutations
#Stomach drivers

GC01 = readRDS("data/GC01.rds")

GC01 = GC01 %>% mutate(Chromosome =  str_c("chr", Chromosome) ) %>%
  mutate(mutid = str_c(Chromosome, Start_Position, Reference_Allele, sep = ":")) 


drivers = read.delim("data/IntOGen-Drivers-Cancer_Genes.tsv", header = T) %>%
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


############################################################################################################################

message("<2> Analysis using clonevol")

output = "Met2"

saveRDS(mutinfo, file = sprintf("%s/%s.mutinfo.rds", output))

sampleNames = c( paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                 paste0("OveryRM_", 1:6),paste0("UterusM_", c(1, 3)) )


vaf.col.names = sampleNames

show.ccf = F

mutdata


pp <- plot.variant.clusters(mutdata,  
                            show.cluster.size = F,
                            show.cluster.label= F,
                            cluster.col.name = cluster.col.name, 
                            vaf.col.names = vaf.col.names,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.shape = 1,
                            variant.class.col.name =cluster.col.name,
                            #vaf.limits = 
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

#if(show.ccf){max.limits = 180 }else{max.limits = 80}
max.limits = 110

pp[[1]] = pp[[1]] +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(position = "right", limits  = c(0, max.limits)) +
  theme(axis.text.x  = element_text(angle = 270, size = 8) )

for(ii in 2:length(pp)){
  pp[[ii ]] = pp[[ii]] +
    scale_y_continuous(position = "right", limits  = c(0, max.limits)) +
    theme(axis.text.x  = element_text(angle = 270, size = 8 ),
          axis.text.y = element_blank()
    )
}


#annotate the driver genes.

if(any(mutdata$is.driver) ){
  
  labels = pp[[1]]$data
  #select colors for annotation.
  clone.colors.sel =unique( clone.colors[ labels$cluster[labels$is.driver]])
  
  pp1 = labels %>% filter(is.driver) %>%
    mutate(cluster_1 = factor(cluster)) %>%
    ggplot(aes(y = cluster, x = 1, label = gene_site)) + theme_classic() +
    geom_point(aes( color = cluster_1), size = 3)+
    geom_text_repel(
      nudge_x      = 0.15,
      direction    = "y",
      hjust        = 0,
      segment.size = 0.2,
      size = 3,
    ) +
    ylim(0,  length(clone.colors) +1 ) +
    xlim(1, 0.8) +
    scale_color_manual( values= clone.colors.sel, guide=FALSE) +
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title   = element_text(hjust = 0.5)
    )
  
  pp[[ length(pp) +1 ]] = pp1
  
}

pdf( sprintf( "%s/%s.box.pdf", output, output), width = 2*length(pp), height = 4)
ggpubr::ggarrange(plotlist = pp, ncol = length(pp), align  = "h" )
dev.off()

#Plotting mean/median of clusters across samples (cluster flow)
pdf( sprintf( "%s/%s.flow.pdf", output, output) , width=8, height=5 )

plot.cluster.flow(mutdata,
                  cluster.col.name = cluster.col.name,
                  vaf.col.names = vaf.col.names,
                  sample.names =vaf.col.names,
                  colors = clone.colors,
                  y.title = "Variant Allele Frequency %"
) +
  theme(axis.text.x = element_text(angle = 90) )

dev.off()

#PLot Pairs

#plot.pairwise(mutdata, col.names = vaf.col.names,
#              group.col.name = cluster.col.name,
#              onePage = FALSE,
#              multiPages = TRUE,
#              out.prefix = sprintf("%s.pair", output),
#              yMaxSmall = 100, xMaxSmall = 110,
#              colors = clone.colors)


##(4.2) Clonal Structure inference and visualization

source("data/auxiliary_functions.R")

sample.groups = mapply(function(x) x[1], strsplit(vaf.col.names, "_")  )
names(sample.groups) = vaf.col.names


# Note:
# we modify the frequency of C1 in Coad_1(same with Coad_2) and C4 in OveryLM_3(same with Coad_5), thus all cluster and samples are included in the finial results.

#mutdata[mutdata$cluster == "4", "OveryLM_3"] = mutdata[mutdata$cluster == "4", "Coad_5"]
#mutdata[mutdata$cluster == "1", "Coad_1"] = mutdata[mutdata$cluster == "1", "Coad_2"]



#sel = c(2:5, 6:10, 11:13, 15, 16, 17, 18) 
sel = c(1, 2:5, 6:10, 11:13, 14, 15, 16, 17, 18) 


y = plot.tree( output  = output , mutdata = mutdata, 
               ccf.col.names = vaf.col.names[sel],
               sample.groups = sample.groups[sel],
               cancer.initiation.model = "monoclonal",
               founding.cluster = 1 , ignore.clusters = 4 , cluster.col.name = "cluster",
               subclonal.test.model = "non-parametric", 
               #scale.monoclonal.cell.frac = FALSE,
               sum.p = 0.01, alpha = 0.05, weighted = FALSE,
               plt.pairwise = F,
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


plot.clonal.models(y,
                   # box plot parameters
                   box.plot = FALSE,
                   fancy.boxplot = FALSE,
                   
                   # bell plot parameters
                   clone.shape = "bell",
                   bell.event = TRUE,
                   bell.event.label.color = "blue",
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   
                   # node-based consensus tree parameters
                   merged.tree.plot = FALSE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = "circle",
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = FALSE,
                   
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = "black",
                   clone.grouping = "horizontal",
                   
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   
                   # output figure parameters
                   out.dir = sprintf("%s/%s", output, "merge"),
                   out.format = "pdf",
                   overwrite.output = TRUE,
                   width = 9,
                   #height = ,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(4,2))


save(mutdata, clone.colors, cluster.col.name , output, sampleNames, y, sample.groups, file = sprintf("%s/clonevol.data.rda", output) )





## Merge multi samples into a meta sample

#Merge multiple samples into a single sample. This can be used to merge multi-region samples to one sample representing the tumor where the regions are taken from

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


plot.clonal.models(y.merge,
                   # box plot parameters
                   box.plot = FALSE,
                   fancy.boxplot = FALSE,
                   
                   # bell plot parameters
                   clone.shape = "bell",
                   bell.event = TRUE,
                   bell.event.label.color = "blue",
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = "circle",
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = FALSE,
                   
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = "black",
                   clone.grouping = "horizontal",
                   
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   
                   # output figure parameters
                   out.dir = sprintf("%s/%s", output, "merge1"),
                   out.format = "pdf",
                   overwrite.output = TRUE,
                   width = 9,
                   #height = ,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(4,4,4))


#merge MRS data into one.
save(y.merge , file = sprintf("%s/%s.y.merge.rda", output, output))



