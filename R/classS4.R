# Define the S4 Class

#' Segment Class
#' @slot data `data.table` of segment file containing CNA information.
#' @slot sample.inof `data.frame` of sample information per patient.
#' @slot ref.build human reference genome version. Default 'hg19'. Optional: 'hg18' or 'hg38'.
#' @slot allele Indicate whether this is allele-specific CNAs. Default: TRUE.
#' @export

Segment <- setClass(
  Class = "Seg",
  slots = c(
    data = "data.table",
    sample.info = "data.frame",
    ref.build = "character",
    allele = "character"
  )
)
