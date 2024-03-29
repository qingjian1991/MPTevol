#' Split the segment regions into several parts
#'
#' According to their shared status.
#' The function first obtains the common shared regions across samples.
#' The corresponding A allele and B allele are output as the format requirements of MEDICC.
#'
#' @param segfiles The allele-specific copy number alterations files generated by **sequenza**.
#' @param sampleid the corresponding sample ids.
#' @param out.dir output dir.
#' @param project.names the project names used in the output.
#' @param N.baf quality control for the sequenza output.
#' @param cnv_min_length quality control for the sequenza output.
#' @param max_CNt quality control for the sequenza output.
#' @param minLength output control: the min length of CNVs to output.
#' @param maxCNV output control: the max CNV to output. When the raw CNV greater than maxCNV, then its value was set to maxCNV.
#' @param medicc.py the position of meidcc.py.
#' @param python the position of meidcc.py.
#'
#' @details
#' This function takes the **sequenza** results as the input and outputs
#' the format requirements of MEDICC.
#'
#' @import GenomicRanges
#' @import ComplexHeatmap
#' @export
splitSegment <- function(segfiles,
                         sampleid,
                         project.names = "tumor",
                         out.dir = "data",
                         N.baf = 30, cnv_min_length = 1e5, max_CNt = 15,
                         minLength = 1e5, maxCNV = 4,
                         medicc.py = "medicc.py",
                         python = "python") {
  seglist <- list()
  # read segs
  for (i in 1:length(sampleid)) {
    seg <- read.delim(file = segfiles[i], header = T, stringsAsFactors = F) %>%
      # dplyr::filter(chromosome == "chr1") %>%
      dplyr::mutate(sample = sampleid[i])

    seglist[[i]] <- seg
  }

  # removing low-confidence regions
  # seglist = base::Reduce(rbind, seglist) %>%
  #  dplyr::filter(N.BAF >= N.baf &  (end.pos - start.pos) >=  cnv_min_length & CNt <= max_CNt)

  # change low-confidence regions into 2,1,1
  seglist <- base::Reduce(rbind, seglist) %>%
    dplyr::mutate(keep = ifelse(N.BAF >= N.baf &
      (end.pos - start.pos) >= cnv_min_length &
      CNt <= max_CNt, TRUE, FALSE)) %>%
    dplyr::mutate(
      CNt = ifelse(keep, CNt, 2),
      A = ifelse(keep, A, 1),
      B = ifelse(keep, B, 1)
    )

  gseg <- GenomicRanges::GRanges(
    seqnames = seglist$chromosome,
    ranges = IRanges::IRanges(seglist$start.pos, seglist$end.pos),
    strand = "+"
  )

  # metadata columns can be added to a GRanges object
  GenomicRanges::mcols(gseg) <- seglist

  # split regions into small regions.
  segdis <- GenomicRanges::disjoin(gseg)
  # add region infor to data
  GenomicRanges::mcols(segdis) <- data.frame(segdis)

  # set the minimal site of each segs.
  # minLength = 1e5

  # get overlaps regions
  overlaps <- GenomicRanges::findOverlaps(segdis, gseg)

  # combined information
  merge <- cbind(
    GenomicRanges::mcols(segdis[queryHits(overlaps), ]),
    GenomicRanges::mcols(gseg[subjectHits(overlaps), ])
  ) %>%
    base::as.data.frame() %>%
    # filter seg length
    dplyr::filter(width >= minLength) %>%
    dplyr::select(seqnames, start, end, width, CNt, A, B, sample)

  # removing neutral region. CNt = 2 and A = 1 and B =1.
  merge_summary <- merge %>%
    dplyr::group_by(seqnames, start, end, width, CNt, A, B) %>%
    dplyr::summarise(num = dplyr::n()) %>%
    dplyr::filter(num == length(sampleid))

  merge <- dplyr::left_join(merge, merge_summary) %>%
    dplyr::filter(is.na(num)) %>%
    dplyr::mutate(num = NULL)

  # major info: A
  merge_A <- merge %>%
    dplyr::select(seqnames, start, end, width, A, sample) %>%
    tidyr::spread(key = sample, value = A, fill = 1)
  # change cnv number when it >=10
  tmp <- merge_A[, 5:ncol(merge_A)]
  tmp[tmp >= maxCNV] <- maxCNV
  merge_A <- cbind(merge_A[, 1:4], tmp)

  # minor info: B
  merge_B <- merge %>%
    dplyr::select(seqnames, start, end, width, B, sample) %>%
    tidyr::spread(key = sample, value = B, fill = 1)

  tmp <- merge_B[, 5:ncol(merge_B)]
  tmp[tmp >= maxCNV] <- maxCNV
  merge_B <- cbind(merge_B[, 1:4], tmp)

  if (!file.exists(out.dir)) {
    system(paste0("mkdir ", out.dir))
  }

  file <- sprintf("%s/%s.descr.txt", out.dir, project.names)

  # output
  write.table(merge_A, file = sprintf("%s/%s.major.txt", out.dir, project.names), quote = F, row.names = F, sep = "\t")
  write.table(merge_B, file = sprintf("%s/%s.minor.txt", out.dir, project.names), quote = F, row.names = F, sep = "\t")

  write.fasta(merge_A = merge_A, major = "major", out.dir = out.dir, project.names = project.names)
  write.fasta(merge_A = merge_B, major = "minor", out.dir = out.dir, project.names = project.names)

  # plot the heatmaps to provide better quality controls
  major <- merge_A %>%
    dplyr::mutate(seq = stringr::str_c(seqnames, start, end, sep = "_")) %>%
    tibble::column_to_rownames(var = "seq") %>%
    dplyr::mutate(seqnames = NULL, start = NULL, end = NULL, width = NULL)

  minor <- merge_B %>%
    dplyr::mutate(seq = stringr::str_c(seqnames, start, end, sep = "_")) %>%
    tibble::column_to_rownames(var = "seq") %>%
    dplyr::mutate(seqnames = NULL, start = NULL, end = NULL, width = NULL)

  plist <- list()
  plist$major <- ComplexHeatmap::Heatmap(major,
    row_names_gp = grid::gpar(fontsize = 6)
  )

  plist$minor <- ComplexHeatmap::Heatmap(minor,
    row_names_gp = grid::gpar(fontsize = 6)
  )

  message(
    sprintf("nohup %s %s %s/%s.descr.txt %s/%s.run -v >%s.run.info.txt &", python, medicc.py, out.dir, project.names, out.dir, project.names, project.names)
  )

  return(list(
    major = merge_A,
    minor = merge_B,
    plist = plist
  ))
}

#' write.fasta
#'
#' Prepare the formate of MEDICC input.
#'
#' @param project.names the project names used in the output.
write.fasta <- function(merge_A, major = "major", out.dir = "data", project.names = "tumor") {
  system(sprintf("mkdir %s", out.dir))
  chrs <- unique(merge_A$seqnames)[!(unique(merge_A$seqnames) %in% c("chrX", "chrY"))]
  num <- ncol(merge_A)

  for (i in chrs) {
    merge_B <- subset(merge_A, seqnames == i)
    text <- rep("0", 2 * (num - 3))
    text[1] <- ">diploid"
    text[2] <- paste0(rep(1, nrow(merge_B)), collapse = "")

    for (j in 5:num) {
      text[2 * (j - 3) - 1] <- sprintf(">%s", colnames(merge_B)[j])
      text[2 * (j - 3)] <- merge_B[, j] %>% stringr::str_c(collapse = "")
    }
    write.table(text, file = sprintf("%s/%s_%s_%s.fasta", out.dir, project.names, major, i), quote = F, col.names = F, row.names = F)
  }

  # write desc
  chrinfo <- rep("0", length(chrs))
  for (i in 1:length(chrs)) {
    chrinfo[i] <- sprintf(
      "%s %s %s", chrs[i],
      sprintf("%s_%s_%s.fasta", project.names, "major", chrs[i]),
      sprintf("%s_%s_%s.fasta", project.names, "minor", chrs[i])
    )
  }
  write.table(chrinfo, file = sprintf("%s/%s.descr.txt", out.dir, project.names), quote = F, col.names = F, row.names = F)
}
