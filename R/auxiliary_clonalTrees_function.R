# save the functions related to clonal trees.

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

#' plotVafCluster
#'
#' Plot variant clustering in each sample by using combination of box,
#' violin and jitter plots.
#'
#' @param variants data frame of the variants.
#' At least cluster column and VAF or CCF columns are required.
#' Cluster column should contain cluster identities as continuous integer values
#' starting from 1.
#' @param cluster.col.name the column names that containing cluster
#' information (Default = "cluster").
#' @param vaf.col.names the column names of samples containing VAF.
#' @param violin whether plotting violin (Default = FALSE).
#' @param box whether plotting box (Default = TRUE).
#' @param jitter whether plotting jitter plot (Default = TRUE).
#' @param founding.cluster the name of founding clones, one of the most important parameters. For most of circumstances, the founding cluster is the cluster with the highest average CCF cluster.
#' @param clone.colors setting clone colors.
#' @param highlight column name to indicate whether highlight the sites (TRUE or FALSE).
#' @param highlight.note.col.name highlight context.
#' @param output.file the output file name (Default = NULL)
#' @import clonevol
#' @return a ggplot object
#' @export
plotVafCluster <- function(variants,
                           cluster.col.name = "cluster",
                           vaf.col.names,
                           clone.colors = NULL,
                           violin = FALSE,
                           box = TRUE,
                           jitter = TRUE,
                           founding.cluster = 1,
                           output.file = NULL,
                           highlight = NULL,
                           highlight.note.col.name = NULL) {
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

  pp <- clonevol::plot.variant.clusters(variants,
    show.cluster.size = F,
    show.cluster.label = F,
    cluster.col.name = cluster.col.name,
    vaf.col.names = vaf.col.names,
    violin = violin,
    box = box,
    jitter = jitter,
    jitter.shape = 1,
    variant.class.col.name = cluster.col.name,
    # vaf.limits =
    jitter.color = clone.colors,
    jitter.size = 1.2,
    jitter.alpha = 1,
    jitter.width = 0.2,
    jitter.center.method = "median",
    jitter.center.size = 1,
    jitter.center.color = "darkgray",
    jitter.center.display.value = "none",
    display.plot = FALSE,
    horizontal = TRUE,
    order.by.total.vaf = F,
    highlight = highlight,
    highlight.shape = 21,
    highlight.color = "blue",
    highlight.fill.color = "green",
    highlight.size = 2.5,
    highlight.note.col.name = NULL,
    highlight.note.size = 2,
    highlight.note.color = "blue",
    highlight.note.angle = 0,
    founding.cluster = founding.cluster,
    ccf = FALSE
  )

  # add annotation of driver genes.
  if (!is.null(highlight)) {
    if (any(mutdata[, highlight])) {
      labels <- pp[[1]]$data
      # select colors for annotation.
      clone.colors.sel <- unique(clone.colors[labels[, cluster.col.name][labels[, highlight]]])

      pp1 <- labels %>%
        dplyr::filter(is.driver) %>%
        dplyr::mutate(cluster_1 = factor(cluster)) %>%
        ggplot2::ggplot(ggplot2::aes(y = cluster, x = 1, label = gene_site)) +
        ggplot2::theme_classic() +
        ggplot2::geom_point(ggplot2::aes(color = cluster_1), size = 3) +
        ggrepel::geom_text_repel(
          nudge_x = 0.15,
          direction = "y",
          hjust = 0,
          segment.size = 0.2,
          size = 3,
        ) +
        ggplot2::ylim(0, length(clone.colors) + 1) +
        ggplot2::xlim(1, 0.8) +
        ggplot2::scale_color_manual(values = clone.colors.sel, guide = "none") +
        ggplot2::theme(
          axis.line = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5)
        )

      pp[[length(pp) + 1]] <- pp1
    }
  }

  if (!is.null(output.file)) {
    pdf(output.file, width = 2 * length(pp), height = 4)
    ggpubr::ggarrange(plotlist = pp, ncol = length(pp), align = "h")
    dev.off()
  }

  ggpubr::ggarrange(plotlist = pp, ncol = length(pp), align = "h")
}

#' tree2timescape
#'
#' This function generates the input of timescape to visual the fisher plot of
#' clonal evolution by using the results of [inferClonalTree()].
#'
#' @param results the clonal trees that generated by [inferClonalTree()].
#' @param samples the samples to show in the fisher plot.
#'
#' @import clonevol
#' @export
tree2timescape <- function(results, samples = NULL) {
  if (is.null(samples)) {
    samples <- names(results$models)
  } else {
    if (!all(samples %in% names(results$models))) {
      stop("check input samplesNames : ", samples[!samples %in% names(results$models)])
    }
  }

  # store the clonevol results in a list
  res <- list(
    samples = samples, clonevol.clone.names = NULL, clonevol.clone.colors = NULL,
    timepoints = seq(1, length(samples)), num.models = nrow(results$matched$index),
    parents = list(), cell.fractions = list(), all = list()
  )

  clonevol.clone.names <- NULL
  clone.nums <- NULL
  clonevol.clone.colors <- NULL

  # create the needed inputs to fishplot
  for (i in 1:nrow(results$matched$index)) {
    vv <- NULL
    for (s in samples) {
      v <- results$models[[s]][[results$matched$index[i, s]]]
      # if (rescale){v = rescale.vaf(v)}
      v <- clonevol:::rescale.vaf(v)
      v <- v[, c("lab", "vaf", "parent", "color")]

      ## scale vaf and make cell.frac
      max.vaf <- max(v$vaf)
      scale <- 0.5 / max.vaf * 2 * 100
      v$vaf <- v$vaf * scale
      v$vaf[v$vaf > 100] <- 100 # safeguard against rounding error making some vaf slightly > 100

      colnames(v) <- c("clone", s, "parent", "color")
      v <- v[!is.na(v$parent) & v$clone != "0", ]
      if (is.null(vv)) {
        vv <- v
      } else {
        vv <- merge(vv, v, all = TRUE)
      }
    }
    for (s in samples) {
      vv[is.na(vv[[s]]), s] <- 0
    }
    vv <- vv[order(as.integer(vv$clone)), ]
    vv$parent[vv$parent == "-1"] <- 0
    rownames(vv) <- vv$clone

    ## fishplot requires clones to be named in sequential order. Do that, but
    ## store the clonevol-generated names and colors for pass-through
    if (is.null(clone.nums)) {
      clone.nums <- c(0, seq(1, nrow(vv)))
      names(clone.nums) <- c(0, vv$clone)

      clonevol.clone.names <- names(clone.nums)
      names(clonevol.clone.names) <- as.character(clone.nums)
      res$clonevol.clone.names <- clonevol.clone.names[-1]

      clonevol.clone.colors <- c("white", vv$color)
      names(clonevol.clone.colors) <- as.character(clone.nums)
      res$clonevol.clone.colors <- clonevol.clone.colors[-1]
    }
    vv$clone <- clone.nums[vv$clone]
    vv$parent <- clone.nums[vv$parent]

    par <- vv$parent
    frac <- vv[, samples]
    res$parents[[i]] <- par
    res$cell.fractions[[i]] <- as.matrix(frac)
    res$all[[i]] <- vv
  }

  ##############################
  # for timescape input

  times <- list(
    clonal_prev = list(),
    tree_edges = list(),
    clone_colours = list()
  )

  for (i in 1:res$num.models) {
    # re-set ancestor clonal prev

    # get sum of prev of certain clone.
    clonal_prev_ancestor <- res$all[[i]] %>%
      dplyr::mutate(clone = rownames(.)) %>%
      data.table::melt(
        id.vars = c("clone", "parent"),
        measure.vars = samples
      ) %>%
      dplyr::group_by(variable, parent) %>%
      dplyr::summarise(sumvalue = sum(value)) %>%
      dplyr::rename(clone = parent) %>%
      dplyr::mutate(clone = as.character(clone))

    # prev = curent - ancestor
    times$clonal_prev[[i]] <- res$all[[i]] %>%
      dplyr::mutate(clone = rownames(.)) %>%
      data.table::melt(
        id.vars = c("clone", "parent"),
        measure.vars = samples
      ) %>%
      dplyr::left_join(clonal_prev_ancestor) %>%
      dplyr::mutate(
        sumvalue = ifelse(is.na(sumvalue), 0, sumvalue),
        value1 = value - sumvalue,
        # set value1 = 0 if value1 <=0
        value1 = ifelse(value1 < 0, 0, value1)
      ) %>%
      dplyr::select(clone, variable, value1) %>%
      dplyr::rename(
        clone_id = clone,
        timepoint = variable,
        clonal_prev = value1
      ) %>% # set arrange of samples.
      dplyr::mutate(timepoint = factor(timepoint, levels = samples)) %>%
      dplyr::arrange(timepoint) %>%
      dplyr::mutate(timepoint = as.character(timepoint))

    # re-mapping clone ids.
    cloneNames <- setNames(
      c(0, rownames(res$all[[i]])),
      nm = c(0, res$all[[i]]$clone)
    )

    times$tree_edges[[i]] <- data.frame(
      source = cloneNames[as.character(res$parents[[i]])],
      target = rownames(res$all[[i]])
    ) %>%
      filter(source != "0")

    times$clone_colours[[i]] <- data.frame(
      clone_id = rownames(res$all[[i]]),
      colour = res$all[[i]]$color
    )
  }

  times
}
