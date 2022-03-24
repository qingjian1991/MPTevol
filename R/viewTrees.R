#' Visualize the trees
#'
#' @param phyloTree phyloTree: The tree is in Parenthetic format.
#' @param tree.format the format of tree, S4 or list. Default is S4.
#' @param normal.node the sample name of normal sample in the tree.
#' @param group  a list that used to indicate the sample groups.
#' @param group.colors an array indicates the colors of sample groups.
#' @param showBootstrap whether showing the bootstrap values. Default is TRUE.
#' @param title title of the plot.
#' @param hexpand_ratio hexpand ratio. see \code{\link[ggtree]{hexpand}}
#'
#' @examples
#' # This dist file is the output of MEDICC
#' dist <- system.file(package = "MPTevol", "extdata", "tree_final.dist")
#'
#' # plot CNA trees without colored samples.
#' plotCNAtree(dist = dist)
#'
#' # create a list to indicate the sample groups.
#' grp <- list(
#'   NORMAL = "NORMAL",
#'   Breast = paste0("Breast_", 1:5),
#'   Coad = paste0("Coad_", 1:5),
#'   Lung = paste0("Lung_", 1:5),
#'   OveryLM = paste0("OveryLM_", 1:5),
#'   OveryRM = paste0("OveryRM_", 1:6),
#'   UterusM = paste0("UterusM_", c(1:7))
#' )
#'
#' plotCNAtree(dist = dist, grp = grp)
#' @return a ggtree object
#' @export
viewTrees <- function(phyloTree,
                      tree.format = "S4",
                      normal.node = "NORMAL",
                      group = NULL,
                      group.colors = NULL,
                      title = "Cancer",
                      showBootstrap = TRUE,
                      hexpand_ratio = 0.3) {
  if (tree.format == "S4") {
    mtree <- phyloTree@tree
  } else {
    mtree <- phyloTree$tree
  }

  # all.length = mtree$edge.length
  root.length <- rev(mtree$edge.length)[1]

  # set outgroup and removing the Normal
  mtree <- treeio::root(mtree, outgroup = normal.node)
  mtree <- treeio::drop.tip(mtree, normal.node)

  # grp <- list(
  #  ACA =  mtree$tip.label[grepl("Aca", mtree$tip.label)] ,
  #  NEC =  mtree$tip.label[grepl("Nec", mtree$tip.label)]
  # )

  # combined bootstrap values.

  if (showBootstrap) {
    bootstrap.value <- ifelse(tree.format == "S4",
      phyloTree@bootstrap.value,
      phyloTree$bootstrap.value
    )

    bp2 <- data.frame(
      node = 1:treeio::Nnode(mtree) + treeio::Ntip(mtree),
      bootstrap = bootstrap.value
    )
    mtree <- dplyr::full_join(mtree, bp2, by = "node")
  }

  p_trees <- ggtree::ggtree(mtree, size = 1) +
    ggtree::geom_tiplab(size = 4) +
    # add scale bars.
    ggtree::geom_treescale(fontsize = 4, linesize = 1, x = 0.1) +
    # set root length
    ggtree::geom_rootedge(rootedge = root.length, size = 1, colour = "grey40") +
    ggtree::hexpand(hexpand_ratio, direction = 1)

  if (showBootstrap) {
    p_trees <- p_trees +
      ggtree::geom_nodepoint(ggplot2::aes(fill = cut(bootstrap, c(0, 70, 90, 100))),
        shape = 21, size = 2
      ) +
      ggtree::geom_nodelab(ggplot2::aes(label = round(bootstrap)), hjust = -0.2, size = 3.5) +
      ggplot2::labs(title = title) +
      ggplot2::scale_fill_manual(
        values = c("black", "grey", "white"), guide = "legend",
        name = "Bootstrap Percentage(BP)",
        breaks = c("(90,100]", "(70,90]", "(0,70]"),
        labels = expression(BP >= 90, 70 <= BP * " < 90", BP < 70)
      )
  } else {
    p_trees <- p_trees +
      ggplot2::labs(title = title)
  }

  # color the groups
  if (!is.null(group)) {
    # check samples ids between trees and grp
    if (!identical(sort(purrr::reduce(group, c)), sort(mtree@phylo$tip.label))) {
      stop("the samplenames in grp were not identical to sample names in the tree")
    }

    p_trees <- ggtree::groupOTU(p_trees, group, "Sites")

    # get levels
    Sites <- levels(p_trees$data$Sites)

    if (!is.null(group.colors)) {
      if (!identical(sort(names(group.colors)), sort(levels(p_trees$data$Sites)))) {
        # get levels that were not in the group.colors
        Site1 <- setdiff(levels(p_trees$data$Sites), names(group.colors))

        Site.colors <- setNames(
          c(set.colors(n = length(Site1), rev = T), group.colors),
          nm = c(Site1, names(group.colors))
        )

        Site.colors <- Site.colors[levels(p_trees$data$Sites)]
      } else {
        Site.colors <- group.colors[levels(p_trees$data$Sites)]
      }
    } else {
      Site.colors <- set.colors(n = length(Sites))
    }

    p_trees <- p_trees + ggplot2::aes(color = Sites) +
      ggplot2::scale_colour_manual(
        values = Site.colors
      )
  }

  p_trees <- p_trees +
    ggtree::theme_tree(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  p_trees
}
