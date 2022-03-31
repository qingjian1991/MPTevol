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
