


# To DO


batchClonalTrees <- function(mutdata) {
  message("1. Pre-analysis data checking")


  message("1.1 Plot variants of clusters in each sample")

  cluster.col.name <- "cluster"
  vaf.col.names

  clone.colors


  pp <- plot.variant.clusters(mutdata,
    show.cluster.size = F,
    show.cluster.label = F,
    cluster.col.name = cluster.col.name,
    vaf.col.names = vaf.col.names,
    violin = FALSE,
    box = FALSE,
    jitter = TRUE,
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
    display.plot = F,
    horizontal = T,
    order.by.total.vaf = F,
    highlight = "is.driver",
    highlight.shape = 21,
    highlight.color = "blue",
    highlight.fill.color = "green",
    highlight.size = 2.5,
    highlight.note.col.name = NULL,
    highlight.note.size = 2,
    highlight.note.color = "blue",
    highlight.note.angle = 0,
    founding.cluster = 1,
    ccf = show.ccf
  )
}
