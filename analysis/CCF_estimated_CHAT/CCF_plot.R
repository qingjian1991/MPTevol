#analysis and plot CCF data.

#sampABlist is in CCF_adjust

sampABlist = readRDS("data/sampABlist.rds")


GC01 = readRDS("data/GC01.rds")

GC01 = GC01 %>% mutate(Chromosome =  str_c("chr", Chromosome) ) %>%
  mutate(mutation_id = str_c(Chromosome, Start_Position, Reference_Allele, sep = ":")) 
#get all adjust maf.
ccf.data.dscat = sampABlist$entire %>%
  mutate(mutation_id = str_c(chr , from , ref , sep = ":")) %>%
  dplyr::select(mutation_id, ends_with("ccf") & !starts_with("merge") ) %>%
  column_to_rownames(var = "mutation_id")

colnames(ccf.data.dscat) = str_remove( colnames(ccf.data.dscat), "ccf" )

#Filtering mutations which max CCF across samples <= 0.2
ccf.data.dscat.row = apply(ccf.data.dscat, 1, max)
ccf.data.dscat = ccf.data.dscat[ ccf.data.dscat.row>=0.2, ]

ccf.data.dscat[ccf.data.dscat >1] = 1




#Plot Heatmaps

library(RColorBrewer)
library(pheatmap)

setcolor = colorRampPalette( brewer.pal(n = 7, name = "PuRd") )(7)
pheatmap::pheatmap(ccf.data.dscat, show_rownames = F, color = setcolor )


####################################################################################################
# building phylogenetic tree based on dist between samples. -------------------------------------------------------
####################################################################################################

library(ape)
library(phangorn)
source("phyloTreeObject.R") #get trees from binary data
library(ggtree)



#CCF >=0.2 and Binary
ccf.data.binary = ccf.data.dscat

ccf.data.binary[ccf.data.binary >= 0.05] = 1
ccf.data.binary[ccf.data.binary < 0.05] = 0

treesinfo = getTrees( ccf.data.binary, method = "MP", bootstrap.rep.num = 2000)


#combined bootstrap
tr = treesinfo$tree
bp2 <- data.frame(node=1:Nnode(tr) + Ntip(tr), 
                  bootstrap = treesinfo$bootstrap.value)
tree <- dplyr::full_join(tr, bp2, by="node")

ggtree(tree) + geom_tiplab() +
  geom_nodelab((aes(label = bootstrap)))


ggtree(tree, size=1) +
  geom_tiplab(size=4) +
  geom_treescale(fontsize=6, linesize=1, offset=1) +
  geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
                 shape=21, size=2.5) +
  theme_tree(legend.position=c(0.8, 0.2))+
  scale_fill_manual(values=c("white", "grey", "black"),guide="legend",
                    name="Bootstrap Percentage(BP)",
                    breaks=c("(90,100]", "(70,90]", "(0,70]"),
                    labels=expression(BP>=90,70 <= BP * " < 90", BP < 70))  


#add groups

grp <- list(NORMAL = "NORMAL",
            Breast =  paste0("Breast_", 1:5),
            Coad = paste0("Coad_", 1:5),
            Lung  = paste0("Lung_", 1:5),
            OveryLM = paste0("OveryLM_", 1:5),
            OveryRM = paste0("OveryRM_", 1:6),
            UterusM = paste0("UterusM_", c(1:7))
)

pdf("Phylotree_all_ccf_0.02.pdf", width = 6, height = 8)

p_iris <- ggtree(tree, size=1) +
  geom_tiplab(size=4) +
  geom_treescale(fontsize=4, linesize=1) 

groupOTU(p_iris, grp, 'Sites') + aes(color=Sites)+
  geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
                 shape=21, size=2.5) +
  scale_colour_manual(
    values  = colorScale[29:36]
  )+
  scale_fill_manual(values=c("black", "grey", "white"),guide="legend",
                    name="Bootstrap Percentage(BP)",
                    breaks=c("(90,100]", "(70,90]", "(0,70]"),
                    labels=expression(BP>=90,70 <= BP * " < 90", BP < 70)) +
  theme_tree(legend.position=c(0.8, 0.25))

dev.off()

# cluster mutations -------------------------------------------------------

# to get the cluster, we using the sciclone to cluster mutations.

library(sciClone)

#filter mutations with CCF >= 0.2 in 2 or more samples.
mutNums = apply(ccf.data.dscat, 1,  function(x) sum(x >=0.2) )

ccf.data.info = ccf.data.dscat[mutNums >=3 ,]

dim(ccf.data.dscat)
dim(ccf.data.info)

pheatmap::pheatmap(ccf.data.info, show_rownames = F, color = setcolor, cutree_rows = 12)

#Removing some artificial mutations.
hc = hclust(dist(ccf.data.info))
hccut <- cutree(hc, k = 12)

table(hccut)

hc.cc =  ccf.data.info[which(hccut == 6 | hccut == 8|  hccut == 10|  hccut == 12 ), ]
ccf.data.info[which(hccut == 6 ), ]
ccf.data.info[which(hccut == 9), ]

#removed some false positive mutations.
ccf.data.info = ccf.data.info[ -which(hccut == 6 | hccut == 8|  hccut == 10|  hccut == 12), ]


pheatmap::pheatmap(ccf.data.info, show_rownames = F, color = setcolor, cutree_rows = 10)


#formate mutations for sciclone analysis.
vafslist = list()

for(i in colnames(ccf.data.info)  ){
  
  mutations = c()
  
  message( sprintf("Sample is %s", i ))
  
  for(j in rownames(ccf.data.info)){
    #for mutations
    tmp = GC01 %>%
      filter( mutation_id  == j) %>%
      dplyr::select(Chromosome, Start_Position) %>%
      unique.data.frame()
    
    tmp1 = GC01 %>%
      filter(mutation_id == j & Tumor_Sample_Barcode == i)
    
    if( nrow(tmp1) ==1 ){
      mutations = rbind(mutations,
                        c( tmp1$Chromosome, tmp1$Start_Position, 50, tmp1$Alt_allele_depth, ccf.data.info[j, i])  )
    }else{
      mutations = rbind(mutations, c(tmp$Chromosome, tmp$Start_Position, 50, 0, ccf.data.info[j, i] ))
    }
    
  }
  
  mutations = mutations %>%
    as.data.frame() %>%
    setNames( nm = c("chr","from","NR","NV","VAF") ) %>%
    mutate(chr = as.character(chr),
           from = as.numeric(from),
           NR = as.numeric(NR),
           NV = as.numeric(NV),
           VAF = as.numeric(VAF)
    )
  
  vafslist[[i]] = mutations
  
}



annotation_vcfs = vafslist[[1]] %>%
  select(chr, from) %>%
  mutate(gene_name = rownames(ccf.data.info) )


#cluster
sc = sciClone(vafs=vafslist, 
              sampleNames= names(vafslist),
              minimumDepth=0,
              doClusteringAlongMargins=FALSE,
              annotation = annotation_vcfs,
              maximumClusters=20)

table(sc@clust$cluster.assignments)


output = "all"

saveRDS(sc, file = sprintf("%s.sc.maxC20.cluster13.rds", output) )

writeClusterTable(sc, sprintf("%s.sc.maxC20.cluster13", output)  )
writeClusterSummaryTable(sc, sprintf("%s.sc.maxC20.cluster13.summary", output))


scidata = read.table(file = sprintf("%s.sc.maxC20.cluster13", output, output), sep = "\t", stringsAsFactors = F, header = T) %>%
  select(cluster, gene_name) 

annotation_row = data.frame(
  sciClone = factor(scidata$cluster)
)
rownames(annotation_row) = scidata$gene_name



mutcolors = c("#7fc97f","#fdc086", "#E64B35", "#82166E",
              "#B77B42","#6349B7","#D5017D","#B77562",
              "#88A4FF", "#439F18", "#971D37","#8C9F3C",
              "#6FDCBF","#5AA36E", "#E7371E")

#Visualizing the variant clusters
set.colors = c("#C6C6C6", "#6FDCBF","#5AA36E","#E99692","#B4D985","#EA7D23","#E53CFB","#4B85B7","#8439FA","#BD8736","#B3B371","#A7C7DE","#EE97FC","#57C222","#BFABD0","#44589B", "#794C18", RColorBrewer::brewer.pal(n = 10, name = "Paired") )


ann_colors = list(
  sciClone = setNames( mutcolors[ c(1:4, 6:12, 14,15) ] , nm = 1:13 ),
  kmeans = setNames( mutcolors[ c(1:4,6:11)] , nm = 1:10 )
)


# Specify colors

#removing C13. C13 is spared, removing this sample.

ccf.data.info.subset = ccf.data.info[ annotation_row %>% 
                                        arrange(sciClone) %>% 
                                        filter(sciClone != 13) %>% 
                                        rownames() ,]


pdf("all.sciclone.heatmap.only.pdf", width = 8 , height = 8)
pheatmap::pheatmap(ccf.data.info.subset, 
                   show_rownames = F, color = setcolor, 
                   annotation_row = annotation_row,
                   annotation_colors = ann_colors,
                   cluster_rows = FALSE )
dev.off()


#####clustring mutations with K-means.
#### Using K-means and sciclone to monitor the clustering

cl <- kmeans( ccf.data.info.subset , 10)


annotation_row1 = annotation_row %>% arrange(sciClone) %>% filter(sciClone != 13) %>%
  rownames_to_column(var  = "Genes") %>%
  mutate(kmeans = factor(  as.numeric( cl$cluster ) ))
rownames(annotation_row1) = annotation_row1$Genes
annotation_row1$Genes = NULL

ann_colors = list(
  sciClone = setNames( mutcolors[c(1:4, 6:10)] , nm = 1:12 ),
  kmeans = setNames( mutcolors[c(1:4,6:11)] , nm = 1:10 )
)


pdf("all.kmeans.heatmap.pdf", width = 8 , height = 8)
pheatmap::pheatmap(ccf.data.info.subset[ annotation_row1 %>% arrange(kmeans) %>% rownames() , ], 
                   show_rownames = F, color = setcolor, 
                   annotation_row = annotation_row1,
                   annotation_colors = ann_colors,
                   cluster_rows = FALSE )
dev.off()


pdf("all.sciclone.heatmap.pdf", width = 8 , height = 8)
pheatmap::pheatmap(ccf.data.info.subset, 
                   show_rownames = F, color = setcolor, 
                   annotation_row = annotation_row1,
                   annotation_colors = ann_colors,
                   cluster_rows = FALSE )
dev.off()


###################################################################################################

## analysis in clonevol.

##Note:
### There are three different data, including CCF, mafa and mafc. 
#(1) CCF is the adjust CCF
#(2) mafa is the adjust maf 
#(3) mafc is the original vaf.

#ccf.data.info.subset is for  running clonevol 
mutdata.ccf = ccf.data.info.subset %>%
  rownames_to_column("mutid") %>%
  full_join( annotation_row1  %>% rownames_to_column("mutid") )

# mutdata is for clonevol.

#as the highly frequency of CCFs, we use the raw vafs to clonal order analysis.

#get all adjust maf.
mutdata.vaf = sampABlist$entire %>%
  mutate(mutid = str_c(chr , from , ref , sep = ":")) %>%
  select(mutid, ends_with("mafa") & !starts_with("merge") ) %>%
  right_join( annotation_row1  %>% rownames_to_column("mutid") ) 

colnames(mutdata.vaf) = str_remove( colnames(mutdata.vaf), "mafa" )


#Original VAFs.

#get all adjust maf.
mutdata.vaf1 = sampABlist$entire %>%
  mutate(mutid = str_c(chr , from , ref , sep = ":")) %>%
  select(mutid, ends_with("mafc") & !starts_with("merge") ) %>%
  right_join( annotation_row1  %>% rownames_to_column("mutid") ) 

colnames(mutdata.vaf1) = str_remove( colnames(mutdata.vaf1), "mafc" )

saveRDS(mutdata.ccf, file = "data/mutdata.ccf.rds")
saveRDS(mutdata.vaf, file = "data/mutdata.vaf.rds")
saveRDS(mutdata.vaf1, file = "data/mutdata.vaf1.rds")

####################################

#Now, we should  run the code for clonevol analysis first. Then run code below.

#adjust plot using complexheatmap

#mutinfo: mutinfo is in clonevol, Running clonevol first.

library(ComplexHeatmap)

mutdata.ccf = readRDS(file = "data/mutdata.ccf.rds")
mutdata.vaf = readRDS(file = "data/mutdata.vaf.rds")
mutdata.vaf1 = readRDS(file = "data/mutdata.vaf1.rds")


mutinfo = readRDS(file = "data/mutinfo.rds")



#reorder some data.
# 
# 
 C5.subclones = mutdata.ccf %>%
   filter(sciClone == 5) %>%
   filter(Coad_5 >= 0.20) %>%
   pull(mutid)
#   group_by(kmeans) %>%
#   summarise(
#     OR1 = median(OveryRM_1),
#     OL1 = median(OveryLM_1),
#     OR2 = median(OveryRM_2),
#     OL2 = median(OveryLM_2),
#     OR3 = median(OveryRM_3),
#     OL3 = median(OveryLM_3),
#     OR4 = median(OveryRM_4),
#     OL4 = median(OveryLM_4),
#   )


mutdata.ccf = mutdata.ccf %>%
  mutate(
    Clone1 = sciClone,
    
    Clone1 = ifelse(sciClone == 4 & kmeans == 8,
                    5, ifelse(sciClone == 4 & kmeans != 8, 1, Clone1) ),
    
    Clone1 = ifelse(
      sciClone == 3, 5, Clone1
    ),
    
    Clone1 = ifelse(
      sciClone == 2, 1, Clone1
    ),
    Clone1 = ifelse(mutid %in% C5.subclones, 4, Clone1)
    
  ) %>%
  arrange(
    Clone1
  ) 
  
  
new.label = mutdata.ccf %>%
  select(mutid, sciClone, kmeans, Clone1)

new.label.unique = new.label %>%
  select(Clone1) %>%
  arrange(
    Clone1
  ) %>%
  unique.data.frame() %>%
  mutate(Clone2 = 1:nrow(.)) %>%
  mutate(
    Clone2 = as.factor(Clone2)
  )

new.label = left_join(
  new.label, new.label.unique
)

sampleNames = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                 paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                 paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7)) )

heatdata = right_join(mutinfo, mutdata.ccf) %>%
  left_join(new.label) %>%
  arrange(Clone1)

#Annotation rows clustering.
ra = rowAnnotation(
  #sciClone = heatdata$sciClone,
  #kmens = heatdata$kmeans,
  sciClone1 = heatdata$Clone2,
  col = list(
    #sciClone = setNames( set.colors[1:  length( unique(heatdata$sciClone ) ) ], 
    #                     nm =  unique( heatdata$sciClone )),
    sciClone1 = setNames( set.colors[1: length( unique(heatdata$Clone2 ) ) ], 
                         nm =  unique( heatdata$Clone2 ))
    #kmeans = setNames( set.colors[1:10] , nm = 1:10 )
    
  ),
  #show_legend = c("Sites" = FALSE),
  
  #annotation_legend_param = list( grid_width = unit(3, "mm") ),
  
  #control width of anno
  simple_anno_size = unit(2, "mm"),
  #show the title of anno
  show_annotation_name = FALSE
) 


#Annotation rows driver genes

heatdata_label = heatdata %>%
  mutate(id = 1:nrow(.)) %>%
  filter(is.driver) %>%
  select(gene_site, id)

ba = rowAnnotation(foo = anno_mark(at =  heatdata_label$id, 
                                   labels = heatdata_label$gene_site,
                                   side = "column",
                                   labels_rot = 20,
                                   labels_gp = gpar(fontsize  = 8)
)
)


heatdata %>%
  select( all_of(sampleNames) ) %>%
  as.matrix() %>%
  Heatmap(show_row_names = FALSE,
          cluster_rows = FALSE,
          col = setcolor ,
          name = "CCF",
          right_annotation  = ra,
          left_annotation = ba,
          column_names_gp = gpar(fontsize = 10)
  )

pdf("all.heatmap.Annotations-new.pdf", width = 8, height = 8)

heatdata %>%
  select( all_of(sampleNames) ) %>%
  as.matrix() %>%
  Heatmap(show_row_names = FALSE,
          cluster_rows = FALSE,
          col = setcolor ,
          name = "CCF",
          right_annotation  = ra,
          left_annotation = ba,
          column_names_gp = gpar(fontsize = 10)
  )

dev.off()

saveRDS(new.label, file = "data/new.label.rds")



