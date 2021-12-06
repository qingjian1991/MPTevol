#Parse the medicc to the ggtrees.

library(ggtree)
library(treeio)
library(phangorn)

library(MEDICCquant)

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

#The output of the medicc is the trees with phyloXML formate. In this code, we read phyloXML trees by the ggtree and visualization the trees based on CNV distances.

#ggtree: read tree formate.
#https://yulab-smu.top/treedata-book/chapter1.html
#mtree = read.phyloxml("results-new/panCancer-New/tree_final.xml")

#data = ntree@data

#mtree = read.phyloxml("results-new/panCancer-New/tree_fitch_nc.xml")

#mtree = read.newick("results-new/panCancer-New/tree_final.new")

mtree = read.newick("results/Pan/tree_final.new")

#combined bootstrap
#tr = treesinfo$tree
#bp2 <- data.frame(node=1:Nnode(mtree) + Ntip(mtree), 
#                  bootstrap =bootstrap.value[2:33])
#tree <- dplyr::full_join(mtree, bp2, by="node")




#plot(midpoint(mtree))

#ggtree(mtree)+ geom_tiplab()

#ggtree(mtree, branch.length="none")

#ggtree(mtree, ladderize=FALSE)

#ggtree(mtree, layout="rectangular")

#data =  mtree@data



#ggtree(mtree, layout="slanted")

#ggtree(mtree) + theme_tree2()

#ggtree(mtree) + geom_tiplab() +
#  geom_nodelab((aes(label = bootstrap)))


ggtree(mtree, size=1) +
  geom_tiplab(size=4) +
  geom_treescale(fontsize=6, linesize=1, offset=1) 


#add groups

grp <- list(NORMAL = "diploid",
            Breast =  paste0("Breast_", 1:5),
            Coad = paste0("Coad_", 1:5),
            Lung  = paste0("Lung_", 1:5),
            OveryLM = paste0("OveryLM_", 1:5),
            OveryRM = paste0("OveryRM_", 1:6),
            UterusM = paste0("UterusM_", c(1:7))
)

p_trees <- ggtree(mtree, size=1) +
  geom_tiplab(size=4) +
  geom_treescale(fontsize=4, linesize=1) 

p1 = groupOTU(p_trees, grp, 'Sites') + aes(color=Sites)+
  #geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
  #               shape=21, size=2.5) +
  scale_colour_manual(
    values  = colorScale[29:36]
  )+
  #geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
  #               shape=21, size=2.5) +
  scale_fill_manual(values=c("black", "grey", "white"),guide="legend",
                    name="Bootstrap Percentage(BP)",
                    breaks=c("(90,100]", "(70,90]", "(0,70]"),
                    labels=expression(BP>=90,70 <= BP * " < 90", BP < 70)) +
  theme_tree(legend.position=c(0.8, 0.25))

pdf("pan1.medicc.pdf", width = 8, height = 8)
p1
dev.off()



##################################################################################################

#functions to get the bootstrap values.

medicc.resample.distance.matrices = function(D, niter=100) {
  result = list()
  pb=txtProgressBar(min=0,max=niter,style=3)
  for (s in 1:niter) {
    setTxtProgressBar(pb,s)
    Dnew = as.matrix(D)
    for (i in 1:nrow(Dnew)) {
      for (j in 1:i) {
        Dnew[i,j] = rnorm(1, mean=Dnew[i,j], sd=sqrt(Dnew[i,j]))
        Dnew[j,i] = Dnew[i,j]
      }
    }
    Dnew[Dnew<0]=0
    Dnew=round(Dnew)
    Dnew = as.dist(Dnew)
    result[[s]] = Dnew
  }
  return(result)
}


bootstrap.trees = function(dist, bootstrap.rep.num = 1000, title = "Cancer"){
  
  #using NJ to create a new tree with bootstrap values.
  
  bootstrap.rep.num = 1000
  
  getTrees = function(D){
    matTree <- nj(D)
    root_num <- which(matTree$tip.label == "diploid")
    matTree <- root(matTree, root_num)
    matTree
  }
  
  D = medicc.read.distance.matrix(dist)
  colnames(D) =rownames(D)
  
  
  matTree = getTrees(D)
  plot(matTree)
  #bootstrap
  #resampled trees.
  resampled = medicc.resample.distance.matrices(D, bootstrap.rep.num)
  
  bootstrap.value = prop.clades(matTree,  lapply(resampled, getTrees) , rooted = is.rooted(matTree))/bootstrap.rep.num*100
  
  #bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e){nj(dist.gene(e))},B = bootstrap.rep.num,quiet = TRUE,rooted = TRUE)/(bootstrap.rep.num)*100
  
  plot( matTree, main = title)
  nodelabels(bootstrap.value)
  
  #for MesKit to plot trees.
  matTree$tip.label[which( matTree$tip.label == "diploid") ] = "NORMAL"
  
  phyloTree = list(
    tree = matTree,
    bootstrap.value = bootstrap.value[1:(length(bootstrap.value)-1)],
    patientID = title
  )
  
  phyloTree
  
}


###############################################################################

#running

#get trees and using NJ to get the bootstraps.
phyloTree = bootstrap.trees(dist = "results/Pan/tree_final.dist",
                             title = "Three"
)


#replot the data


#add groups

mtree = phyloTree$tree
bootstrap.value = phyloTree$bootstrap.value

plot( mtree)
nodelabels(bootstrap.value)

#set outgroup and removing the Normal
#mtree = root(mtree, outgroup = "NORMAL")
#mtree= drop.tip(mtree, "NORMAL")


#combined bootstrap
bp2 <- data.frame(node=1:(Nnode(mtree)-1) + Ntip(mtree),
                  bootstrap = bootstrap.value  )
mtree <- dplyr::full_join(mtree, bp2, by="node")



grp <- list(NORMAL = "NORMAL",
            Breast =  paste0("Breast_", 1:5),
            Coad = paste0("Coad_", 1:5),
            Lung  = paste0("Lung_", 1:5),
            OveryLM = paste0("OveryLM_", 1:5),
            OveryRM = paste0("OveryRM_", 1:6),
            UterusM = paste0("UterusM_", c(1:7))
)

p_trees <- ggtree(mtree, size=1) +
  geom_tiplab(size=4) +
  geom_treescale(fontsize=4, linesize=1) 

p1 = groupOTU(p_trees, grp, 'Sites') + aes(color=Sites)+
  geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
                 shape=21, size=2.5) +
  #geom_nodelab( aes(label = round(bootstrap)), hjust  = -0.2, size = 3.5) +
  scale_colour_manual(
    values  = colorScale[29:36]
  )+
  geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
                 shape=21, size=2.5) +
  scale_fill_manual(values=c("black", "grey", "white"),guide="legend",
                    name="Bootstrap Percentage(BP)",
                    breaks=c("(90,100]", "(70,90]", "(0,70]"),
                    labels=expression(BP>=90,70 <= BP * " < 90", BP < 70)) +
  theme_tree(legend.position=c(0.8, 0.25))

p1

pdf("pan1.medicc-add-bootstrap.pdf", width = 8, height = 8)
p1
dev.off()
