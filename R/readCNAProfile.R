#' Read CNA Profiles
#'
#' We used a CNAqc object, containing a set of mutations, CNA calls and tumor purity values. The CNAqc was used to deal with the allele-specific CNAs.
#' @param Patient_ID Patient_ID: select the specific patients. IF not indicate, the input is Maf and seg, or the input is MafList and segList.
#' @param maf Maf or MafList.
#' @param seg seg or seglist.
#' @param purity: purity information for each samples.
#' @param ref: human reference genome version. Default 'hg19'. Optional: 'hg18' or 'hg38'.
#' @details
#' This code reads the CNA Profiles for each patient. The tumor names of  maf and seg are required to match each other.
#'
#'
#' @return cnaqc.list for cnaqc initiation.
#'
#' @export

readCNAProfile = function(
  maf,
  seg,
  Patient_ID = NULL,
  purity = 1,
  ref = "hg19"
){

  if(!is.null(Patient_ID)){
    maf = maf[[Patient_ID]]
    seg = seg[[Patient_ID]]
  }

  #read mutations.
  snvs = maf@data %>%
    mutate( DP = Ref_allele_depth + Alt_allele_depth) %>%
    dplyr::rename(
      chr = Chromosome , from = Start_Position , to = End_Position ,
      ref = Reference_Allele ,
      alt = Tumor_Seq_Allele2 ,
      NV = Alt_allele_depth
    ) %>%
    dplyr::select( chr, from, to,  ref, alt, DP, NV, VAF, dplyr::everything())


  #read CNAs.

  cna = seg %>%
    dplyr::rename(
      chr = Chromosome,
      from = Start_Position,
      to = End_Position,
      Major = Major_CN,
      minor = Minor_CN
    ) %>%
    dplyr::select(chr, from, to, Major, minor,dplyr::everything())

  if(!identical( unique(sort(snvs$Tumor_Sample_Label)),
                 unique(sort(cna$Tumor_Sample_Label))) ){
    stop("The Tumor_Sample_Labels of SNV and CNA are not identical")
  }


  TumorSamples = unique(sort(snvs$Tumor_Sample_Label))

  cnaqc.list = list()

  #Get purity info
  if( is.null(names(purity)) ){
    purity = setNames(
      rep(purity, length(TumorSamples)),
      nm = TumorSamples
    )

  }


  for(i in TumorSamples){
    cnaqc.list[[i]] = init(
      snvs %>% filter(Tumor_Sample_Label == i),
      cna %>% filter(Tumor_Sample_Label == i),
      purity[i],
      ref = ref
    )
  }


  return(cnaqc.list)


}


