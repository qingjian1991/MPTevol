% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calPropDriver.R
\name{calPropDriver}
\alias{calPropDriver}
\title{calPropDriver calculates the proportions of driver mutation for mutation sets.}
\usage{
calPropDriver(
  maf,
  patient.id = NULL,
  driverGene,
  class = "SP",
  classByTumor = FALSE,
  vaf.cutoff = 0.01,
  silent.columns = NULL
)
}
\arguments{
\item{maf}{Maf or MafList object generated by \code{readMaf()} function}

\item{patient.id}{Select the specific patients. Default \code{NULL}, all patients are included.}

\item{driverGene}{The driver Gene names (Genes Symbols)}

\item{class}{The class which would be represented.
"SP" (Shared pattern: Public/Shared/Private), other options: "CS" (Clonal status: Clonal/Subclonl)
and "SPCS". see \code{\link[MesKit:classifyMut]{MesKit::classifyMut()}}.}

\item{classByTumor}{Logical (Default: \code{FALSE}). Classify mutations based on "Tumor_ID".}

\item{vaf.cutoff}{Removing mutations of low variant allele frequency (VAF).}

\item{silent.columns}{The \code{Variant_Classification} field in the MAF files that indicates the silent mutations. Defaults: c("Silent", "3'Flank", "IGR", "Intron", "RNA")}
}
\description{
The mutations are classified by \code{classifyMut()} internally.
}
\examples{

# Get the driver gene.
driverGene <- read.delim(system.file(package = "MPTevol", "extdata", "IntOGen-Drivers-Cancer_Genes.tsv"), header = T) \%>\%
filter(CANCER_TYPE \%in\% c("BRCA", "COREAD", "LUAD", "LUSC")) \%>\%
  pull(SYMBOL) \%>\%
  unique()

prop = calPropDriver(maf, patient.id = "BRCA", driverGene = driverGene)

prop$BRCA$plot


}
