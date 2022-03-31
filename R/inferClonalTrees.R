#' inferClonalTrees
#'
#' This function have two main modules,
#' including inferring the clonal trees and plot the clonal models.
#' @details
#' Inferring the clonal trees is the central process in clonal construction.
#' However, users always find that it is difficult to build clonal trees.
#' Therefore, we should check the cluster structures before building clonal trees.
#' Here are some suggestions about building clonal trees.
#' 1. chose the optimal clustering methods. Before do mutation clustering.
#' We should removing the low-quality mutations.
#' The indels are suggested to be removed.
#' The mutations in the LOH regions are suggested to be removed.
#' The mutations in the cnv-regions are should be carefully checked.
#' 2. chose the right founding cluster.
#' 3. ignore some false-negative clusters. For some clusters,
#' especially clusters that have low vafs in all samples,
#' were probably false-positive clusters. Removing clusters that having too few mutations.
#' 4. try different cutoffs. The two parameters **sum.p**
#' and **alpha** are used to determine whether a cluster is in a sample.
#' A relaxed cutoffs (small values of the two parameters) enables more
#' clusters are though to be present in the sample.
#' @param project.names the project names used in the output.
#' @param variants data frame of the variants.
#' At least cluster column and VAF or CCF columns are required.
#' Cluster column should contain cluster identities as continuous integer values
#' starting from 1.
#' @param vaf.col.names the column names of samples containing VAF.
#' @param ccf.col.names the column names of samples containing CCF.
#' Note: either setting **vaf.col.names** or **ccf.col.names**.
#' @param sample.groups indicate the samples groups.
#' An example is setNames(c("Primary","Primary","Met","Met"), nm = c("P1","P2","M1","M1") )
#' @param founding.cluster the name of founding clones, one of the most important parameters. For most of circumstances, the founding cluster is the cluster with the highest average CCF cluster.
#' @param ignore.clusters the clusters that ignores to analysis.
#' For some clusters, especially clusters that have low vafs in all samples,
#' were probably false-positive clusters.
#' @param cluster.col.name the column names that containing cluster information.
#' @param clone.colors setting clone colors.
#' @param subclonal.test.model What model to use when generating the bootstrap values
#' are: c('non-parametric', 'normal', 'normal-truncated', 'beta', 'beta-binomial').
#' (Default = "non-parametric")
#' @param cancer.initiation.model cancer evolution model to use, c('monoclonal', 'polyclonal').
#' Monoclonal model assumes the orginal tumor (eg. primary tumor) arises from a
#' single normal cell; polyclonal model assumes the original tumor can
#' arise from multiple cells (ie. multiple founding clones).
#' In the polyclonal model, the total VAF of the separate founding clones
#' must not exceed 0.5.
#' @param sum.p min probability that the cluster is non-negative in a sample(Default = 0.05).
#' @param alpha alpha level in confidence interval estimate for the cluster (Default = 0.05).
#' @param weighted weighted model (default = FALSE)
#' @param consensus.tree whether build the consensus tree (Default = TRUE).
#' @param plot.models whether plot the models (Default = TRUE).
#' @param plot.pairwise.CCF whether plot pairwise CCF comparison (Default = FALSE).
#' @param highlight column name to indicate whether highlight the sites (TRUE or FALSE).
#' @param highlight.note.col.name highlight context.
#' @param highlight.CCF highlight is CCF or VAF (Default = FALSE).
#' @import clonevol
#'
#' @export
inferClonalTrees <- function(project.names,
                             variants,
                             vaf.col.names = NULL,
                             ccf.col.names = NULL,
                             sample.groups = NULL,
                             founding.cluster = 1, ignore.clusters = NULL,
                             cluster.col.name = "cluster",
                             clone.colors = NULL,
                             subclonal.test.model = "non-parametric",
                             cancer.initiation.model = "monoclonal",
                             sum.p = 0.05, alpha = 0.05,
                             weighted = FALSE,
                             consensus.tree = TRUE,
                             plot.models = TRUE, # plot.clonal.models
                             plot.pairwise.CCF = F,
                             highlight.note.col.name = NULL,
                             highlight = "is.driver",
                             highlight.CCF = FALSE) {
  variants <- as.data.frame(variants)

  if (is.null(clone.colors)) {
    # Visualizing the variant clusters
    set.colors <- c(
      "#C6C6C6", "#6FDCBF", "#5AA36E", "#E99692",
      "#B4D985", "#EA7D23", "#E53CFB", "#4B85B7",
      "#8439FA", "#BD8736", "#B3B371", "#A7C7DE",
      "#EE97FC", "#57C222", "#BFABD0", "#44589B",
      "#794C18",
      RColorBrewer::brewer.pal(n = 10, name = "Paired")
    )

    clone.colors <- set.colors[1:length(unique(variants[, cluster.col.name]))]
  }

  if (plot.pairwise.CCF) {
    clonevol::plot.pairwise(
      data = variants,
      col.names = vaf.col.names,
      group.col.name = cluster.col.name,
      onePage = FALSE,
      multiPages = TRUE,
      out.prefix = sprintf("%s/%s.pair", project.names, project.names),
      yMaxSmall = 100, xMaxSmall = 120,
      colors = clone.colors
    )
  }
  message("Main Function: Inferring clonal models")

  y.B <- clonevol::infer.clonal.models(
    variants = variants,
    cluster.col.name = cluster.col.name,
    ccf.col.names = ccf.col.names,
    vaf.col.names = vaf.col.names,
    sample.groups = sample.groups,
    cancer.initiation.model = cancer.initiation.model,
    subclonal.test = "bootstrap",
    subclonal.test.model = subclonal.test.model,
    num.boots = 1000,
    founding.cluster = founding.cluster,
    cluster.center = "median",
    ignore.clusters = ignore.clusters,
    clone.colors = clone.colors,
    min.cluster.vaf = 0.01,
    merge.similar.samples = F,
    weighted = weighted,
    # min probability that CCF(clone) is non-negative
    sum.p = sum.p,
    # alpha level in confidence interval estimate for CCF(clone)
    alpha = alpha
  )

  if (consensus.tree) {
    message("convert.consensus.tree.clone.to.branch")
    # Converting node-based trees to branch-based trees
    y.B <- clonevol::convert.consensus.tree.clone.to.branch(
      y.B,
      cluster.col = cluster.col.name,
      branch.scale = "sqrt"
    )
  }

  if (plot.models) {
    message("plot.clonal.models")
    clonevol::plot.clonal.models(y.B,
                                 # box plot parameters
                                 box.plot = TRUE,
                                 fancy.boxplot = TRUE,
                                 fancy.variant.boxplot.highlight = "is.driver",
                                 fancy.variant.boxplot.highlight.size = 2.5,
                                 fancy.variant.boxplot.highlight.shape = 21,
                                 fancy.variant.boxplot.highlight.color = "blue",
                                 fancy.variant.boxplot.highlight.fill.color = "green",
                                 fancy.variant.boxplot.highlight.note.col.name = highlight.note.col.name,
                                 fancy.variant.boxplot.highlight.note.color = "blue",
                                 fancy.variant.boxplot.highlight.note.size = 2,
                                 fancy.variant.boxplot.jitter.alpha = 1,
                                 fancy.variant.boxplot.jitter.center.color = "grey50",
                                 fancy.variant.boxplot.base_size = 12,
                                 fancy.variant.boxplot.plot.margin = 1,
                                 fancy.variant.boxplot.vaf.suffix = ".VAF",
                                 fancy.variant.boxplot.founding.cluster = founding.cluster,
                                 fancy.variant.boxplot.ccf = highlight.CCF,

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
                                 merged.tree.clone.as.branch = TRUE,
                                 mtcab.event.sep.char = ",",
                                 mtcab.branch.text.size = 1,
                                 mtcab.branch.width = 0.3,
                                 mtcab.node.size = 3,
                                 mtcab.node.label.size = 1,
                                 mtcab.node.text.size = 1.5,
                                 # cellular population parameters
                                 cell.plot = TRUE,
                                 num.cells = 100,
                                 cell.border.size = 0.25,
                                 cell.border.color = "black",
                                 clone.grouping = "horizontal",

                                 # meta-parameters
                                 scale.monoclonal.cell.frac = TRUE,
                                 show.score = FALSE,
                                 cell.frac.ci = TRUE,
                                 disable.cell.frac = FALSE,
                                 # output figure parameters
                                 out.dir = sprintf("%s", project.names),
                                 out.format = "pdf",
                                 overwrite.output = TRUE,
                                 width = 15,
                                 # height = 4,
                                 # vector of width scales for each panel from left to right
                                 panel.widths = c(3, 4, 2, 4, 4)
    )


    pdf(
      file = sprintf("%s/%s.trees.pdf", project.names, project.names),
      width = 6, height = 6
    )
    clonevol::plot.all.trees.clone.as.branch(y.B,
                                             branch.width = 0.5,
                                             node.size = 2, node.label.size = 0.5,
                                             tree.rotation = 180
    )
    dev.off()
  }

  return(y.B)
}
