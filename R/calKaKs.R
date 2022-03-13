
#' calKaKs
#'
#' calculate the Ka/Ks of each group.
#'
#' @param maf: maf data
#' @param VAF_cutoff: Removing mutations with low variant allele frequency (VAF).
#' @param class: The class which would be represented.
#' @param parallel: Whether calculate the Ka/Ks in parallel? Default: TRUE
#' "SP" (Shared pattern: Public/Shared/Private),
#' other options: "CS" (Clonal status: Clonal/Subclonl) and "SPCS". see MesKit classifyMut.
#'
#'
#' @export

calKaKs = function(maf,
                   vaf_cutoff = 0.05,
                   class = "SP",
                   parallel = TRUE
                   ){

  # To do: be careful about the samples and tumors.

  #Get the mutation groups.

  maf_input = subMaf(maf, mafObj = FALSE, use.tumorSampleLabel = TRUE)

  #get mutation classifications.
  maf_class = classifyMut(maf, patient.id = NULL, class = class , classByTumor = FALSE)


  #Merge the maf input and mutation class
  maf_merge = list()

  for(i in names(maf_input) ){
    message(i)
    maf_merge[[i]] = maf_input[[i]] %>%
      mutate(Mut_ID = str_c(Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":" ) ) %>%
      left_join(
        maf_class[[i]]
      ) %>%
      dplyr::select(
        Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Mutation_Type, Patient_ID, Tumor_ID, Variant_Classification, VAF
      )

  }

  maf_merge <- purrr::reduce(maf_merge, rbind)


  #cols_types = setNames(c("#00A087FF", "#3C5488FF", "#8491B4FF", "#F39B7FFF", "#E64B35FF", "#4DBBD5FF"),
  #                      nm = c("Public_Clonal","Public_Subclonal","Shared_Clonal","Shared_Subclonal","Private_Clonal","Private_Subclonal") )


  maf_list =  maf_merge %>%
    split( paste( maf_merge$Tumor_ID, maf_merge$Mutation_Type, sep = ":" ) )


  # We use the easypar to do parallel calculations.

  if(parallel){
    maf_KaKs = easypar::run(FUN = getKaKs,
                        PARAMS = lapply(1:length(maf_list), function(x)
                          list(df = maf_list[[x]], VAF_cutoff = VAF_cutoff) ) ,
                        parallel = TRUE,
                        outfile = NULL,
                        export = NULL,
                        packages = "tidyverse",
                        filter_errors = FALSE
    )

  }else{
    maf_KaKs = lapply(maf_list, getKaKs)
  }

  names(maf_KaKs) = names(maf_list)


  KaKs_data = list()

  for(i in names(maf_KaKs) ){
    KaKs_data[[i]] = maf_KaKs[[i]] %>%
      mutate(type = i)
  }

  KaKs_data = purrr::reduce(KaKs_data, rbind)

  KaKs_data = KaKs_data %>%
    mutate(
      Tumor_ID = mapply(function(x) x[1], str_split(type, ":") ),
      Type = mapply(function(x) x[2], str_split(type, ":") )
    )


  KaKs_data %>%
    filter(name %in% c("wall")) %>%
    mutate(name = factor(name, levels = c("wall") )) %>%
    #mutate(Type = factor(Type, levels = c("Shared_Clonal","Private_Clonal","Private_Subclonal") ))  %>%
    ggplot(aes(x = Tumor_ID , y = mle, fill = Type) ) + ggpubr::theme_pubr() +
    geom_bar(stat = "identity", position =  position_dodge(width = 0.90)) +
    #geom_linerange(aes(ymin = cilow, ymax = cihigh), position =  position_dodge(width = 0.90) ) +
    geom_hline( yintercept = 1, linetype = 2, size = 1) +
    labs(x = NULL, y =  latex2exp::TeX("Dn/Ds ($\\omega_{all}$)") )  +
    #facet_grid( ~ name) +
    #scale_fill_manual(
    #  values = cols_types[levels(maf_merge_sum$Type) ]
    #) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14)
    )

  return(
    KaKs_data
  )


}




#' getKaKs
#'
#' See Ka/Ks between different groups.
#'
#' @param df: data. six columns are required to calculate the Ka/Ks, including "Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2" and "VAF"
#' @param VAF_cutoff: VAF cutoff. Removing mutations with low variant allele frequency (VAF).
#'
#' @details
#'
#' We use the dndscv to calculate the Ka/Ks.
#'
#'
#' @export

getKaKs = function(df, VAF_cutoff = 0.05){

  data(list = sprintf("submod_%s", "13r_3w"), package = "dndscv")

  mutations = df %>%
    filter(VAF >= VAF_cutoff) %>%
    select(c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2"))

  kaks = dndscv::dndscv( mutations = mutations,
                         max_muts_per_gene_per_sample = Inf,
                         max_coding_muts_per_sample = Inf,
                         #gene_list = intersect( Genes_Covered, genes),
                         sm = submod_13r_3w
  )

  kaks$globaldnds

}


