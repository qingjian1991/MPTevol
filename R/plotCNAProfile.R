#' Visualiza CNA profile
#'
#' This function plots the allele-specific CNAs of multiple-samples.
#' See [readCNAProfile()] for examples.
#'
#' @param cnaqc.list cnaqc.list
#' @param min_length_show the minimal length of CNVs to show.
#'
#' @import CNAqc
#'
#' @export
plotCNAProfile <- function(cnaqc.list, min_length_show = 1e5) {
  L <- x
  Ln <- names(L)
  if (is.null(Ln)) {
    Ln <- paste0("Sample ", 1:length(L))
  }
  KARYO_colors <- CNAqc:::get_karyotypes_colors(NULL) # NOTE ::: is an invalid operation when submitted as CRAN/Bioc Package

  # KARYO_colors[1] = "white"
  # KARYO_colors$"3:1" = ""

  KARYO_colors <- setNames(
    c("white", "steelblue", "darkblue", "turquoise4", "#F3BA45", "#F7BCB4", "#EF7969"),
    nm = c("1:1", "1:0", "0:0", "2:0", "2:1", "2:2", "3:1")
  )

  calls <- lapply(Ln, function(s) {
    W <- L[[s]]$cna %>%
      dplyr::mutate(label = paste(Major, minor,
        sep = ":"
      ), CN = minor + Major, sample = s) %>%
      dplyr::select(chr, from, to, label, CN, sample)
    CNAqc:::relative_to_absolute_coordinates(L[[s]], W)
  })
  calls_flat <- suppressWarnings(Reduce(
    function(x, y) {
      dplyr::full_join(x,
        y,
        by = c("chr", "from", "to", "label", "CN", "sample")
      )
    },
    calls
  ) %>% dplyr::mutate(label = ifelse(label %in% names(KARYO_colors),
    label, "other_AMP"
  )))
  KARYO_colors <- c(KARYO_colors, other_AMP = "#9A2414")
  chromosomes <- calls_flat$chr %>% unique()
  reference_genome <- CNAqc:::get_reference(L[[1]]$reference_genome) %>%
    dplyr::filter(chr %in% chromosomes)
  low <- min(reference_genome$from)
  upp <- max(reference_genome$to)

  bl_genome <- suppressMessages(
    blank_genome1( # Where is blank_genome1
    ref = L[[1]]$reference_genome,
    chromosomes = chromosomes,
    label_chr = NA
  ) + ggplot2::labs(x = "", y = ""))

  seg_id <- pio:::nmfy(Ln, seq_along(Ln))
  calls_flat$sample_id <- seg_id[calls_flat$sample]

  calls_flat <- calls_flat %>%
    dplyr::filter(label != "1:1") %>%
    dplyr::filter(to - from >= min_length_show)

  bl_genome +
    ggplot2::geom_segment(data = calls_flat,
                          ggplot2::aes(x = from, xend = to, y = sample_id,
                                       yend = sample_id, color = label), size = 5) +
    ggplot2::scale_color_manual(values = KARYO_colors[2:length(KARYO_colors)]) +
    ggplot2::coord_polar(theta = "x", clip = "off") +
    ggplot2::guides(color = ggplot2::guide_legend("Karyotype", row = 1)) +
    ggplot2::ylim(-4, max(seg_id) + 1) +
    ggplot2::labs(title = "Comparative CNA",
                  subtitle = paste0("Tracks: ", paste(Ln, collapse = ", "))) +
    ggplot2::theme(
      legend.key.height = ggplot2::unit(0.1, "cm"), axis.text.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(size = 0.3)
    )
}
