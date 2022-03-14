
#' viewTrees
#'
#' we used the ggtree to visualize the trees.
#'
#' @param phyloTree phyloTree: The tree is in Parenthetic Format. see \code{\link[read.tree]{subMaf}}.
#'
#' @param tree.format the format of tree, S4 or list. Default is S4.
#' @param normal.node the sample name of normal sample in the tree.
#' @param group  a list that used to indicate the sample groups
#' @param showBootstrap whether showing the bootstrap values. Default is TRUE.
#' @param title: title of the plot.
#' @param hexpand_ratio: hexpand ratio. see \code{\link[ggtree]{hexpand}}
#'
#' @examples
#'
#' #This dist file is the output of MEDICC
#' dist = system.file(package="MPTevol", "extdata", "tree_final.dist")
#'
#' #plot CNA trees without colored samples.
#' plotCNAtree(dist = dist)
#'
#' #create a list to indicate the sample groups.
#' grp <- list(NORMAL = "NORMAL",
#'     Breast =  paste0("Breast_", 1:5),
#'     Coad = paste0("Coad_", 1:5),
#'     Lung  = paste0("Lung_", 1:5),
#'     OveryLM = paste0("OveryLM_", 1:5),
#'     OveryRM = paste0("OveryRM_", 1:6),
#'     UterusM = paste0("UterusM_", c(1:7))
#' )
#'
#' plotCNAtree(dist = dist, grp = grp)
#'
#'
#'
#' @return a ggtree object
#' @export
#'

viewTrees = function(phyloTree,
                     tree.format = "S4",
                     normal.node = "NORMAL",
                     group = NULL,
                     title = "Cancer",
                     showBootstrap = TRUE,
                     hexpand_ratio = 0.3
){

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


  if(tree.format == "S4"){
    mtree = phyloTree@tree
  }else{
    mtree = phyloTree$tree
  }

  #all.length = mtree$edge.length
  root.length = rev(mtree$edge.length)[1]

  #set outgroup and removing the Normal

  mtree = root(mtree, outgroup = normal.node)
  mtree= drop.tip(mtree, normal.node)

  #grp <- list(
  #  ACA =  mtree$tip.label[grepl("Aca", mtree$tip.label)] ,
  #  NEC =  mtree$tip.label[grepl("Nec", mtree$tip.label)]
  #)

  #combined bootstrap values.

  if(showBootstrap){

    bootstrap.value = ifelse( tree.format == "S4", phyloTree@bootstrap.value, phyloTree$bootstrap.value)

    bp2 <- data.frame(node=1:Nnode(mtree) + Ntip(mtree),
                      bootstrap = bootstrap.value)
    mtree <- dplyr::full_join(mtree, bp2, by="node")
  }


  p_trees <- ggtree(mtree, size=1) +
    geom_tiplab(size=4) +
    #add scale bars.
    geom_treescale(fontsize=4, linesize=1, x = 0.1) +
    #set root length
    geom_rootedge(rootedge = root.length, size=1, colour = "grey40") +
    hexpand(hexpand_ratio, direction = 1)

  if(showBootstrap){

    p_trees = p_trees +
      geom_nodepoint(aes(fill=cut(bootstrap, c(0, 70, 90, 100))),
                     shape=21, size=2 ) +
      geom_nodelab( aes(label = round(bootstrap)), hjust  = -0.2, size = 3.5) +
      labs(title = title) +
      scale_fill_manual(values=c("black", "grey", "white"),guide="legend",
                        name="Bootstrap Percentage(BP)",
                        breaks=c("(90,100]", "(70,90]", "(0,70]"),
                        labels=expression(BP>=90,70 <= BP * " < 90", BP < 70))
  }else{
    p_trees = p_trees +
      labs(title = title)

  }

  #color the groups
  if(!is.null(group)){

    #check samples ids between trees and grp
    if( !identical(sort(purrr::reduce(group, c)), sort( mtree@phylo$tip.label) ) ){
      stop("the samplenames in grp were not identical to sample names in the tree")
    }

    p_trees = groupOTU(p_trees, group, 'Sites') + aes(color=Sites)+
      scale_colour_manual(
        values  = colorScale[1:length(group)]
      )
  }

    p_trees = p_trees +
      theme_tree(
        plot.title = element_text(hjust = 0.5)
      )

  p_trees
}











