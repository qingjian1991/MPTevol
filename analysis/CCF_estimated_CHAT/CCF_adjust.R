#adjust CCF using CHAT( From Huzheng)

#Code is from "MRS_analysis_run.R"

library(tidyverse)
library(flow)
library(neutralitytestr)
library(patchwork)


source("data/MRS_analysis_m1.R")


#Purity information.
purity_table = read.table("data/purity.summary.txt", header = T)

#for mutations.
GC01 = readRDS("data/GC01.rds")


cancers = c("Breast","Lung","Coad","OveryLM","OveryRM","UterusM")

sampleNames = c( paste0("Breast_", 1:5), paste0("Lung_", 1:5),
                 paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), 
                 paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7)) )

sampleGroups = list(
            Breast =  paste0("Breast_", 1:5),
            Coad = paste0("Coad_", 1:5),
            Lung  = paste0("Lung_", 1:5),
            OveryLM = paste0("OveryLM_", 1:5),
            OveryRM = paste0("OveryRM_", 1:6),
            UterusM = paste0("UterusM_", c(1:7)),
            Met1 = c( paste0("Coad_", 1:5), 
                      paste0("OveryLM_", 1:5), paste0("OveryRM_", 1:6), 
                      paste0("UterusM_", c(1,3))
                      ),
            entire = c(
              paste0("Breast_", 1:5), paste0("Coad_", 1:5),paste0("Lung_", 1:5),
              paste0("OveryLM_", 1:5),paste0("OveryRM_", 1:6),paste0("UterusM_", c(1:7))
            )
)




###############################################################################################################
##
##Running for single Nec and Aca.
##

#' runSampAB
#' Running for each samples.

runSampAB = function(GC01_mutations, samples, new_samples, normalproportion, 
                     titanPath="/data1/qingjian/Project/QiuMZ/TitanCNA/scripts/snakemake/results/titan/hmm/Seg/",
                     seg.formate = "titan",
                     outputfolder = sprintf("MRS_analysis_%s/", seg.formate),
                     binw = 0,
                     plot.hist = TRUE,
                     t = 0.45
                     
                     ){
  if( !file.exists(outputfolder)){
    system( sprintf("mkdir %s", outputfolder) )
  }
  
  
  #(1) get input
  message("(1) get input")
  sampAB = trans2sampAB(GC01_mutations, 
                        samples = samples , 
                        new_samples = NULL)
  
  #normalproportion = setNames( 1- purity_table[ samples, "titan.adj"] ,nm = samples)
  
  
  #(2) get CCF
  message("(2) get CCFs")
  sampAB = adjust.ccf.titan.multi(sampAB = sampAB, normalproportion = normalproportion, samples = samples , t =  t, titanPath = titanPath,  seg.formate = seg.formate                               
  )
  #(3) classification public or privates.
  if( !is.null(new_samples) & !( all(samples == new_samples) ) ){
    sampAB = sampAB2newNames(sampAB, samples = samples, new_samples = new_samples)
  }else{
    new_samples = samples
  }
  

  sampAB = pubOrSub(sampAB , samples = new_samples)
  
  #filter absent.
  sampAB = sampAB %>%
    filter(!pubOrSub %in% c("absent", "unknown") ) 
  #mutate(pubOrSubadj = ifelse( grepl(pubOrSub, pattern = ",") , "shared" , pubOrSub) )
  
 
  
  if(plot.hist){
    message("(4) get the parameters and plot.")
    #(4) Get the parameters and plot.
    multilist = list()
    multinames = c()
    for(i in 1:(length(new_samples)-1) ){
      for(j in (i+1):length(new_samples)){
        sn1 = new_samples[i]
        sn2 = new_samples[j]
        message(sn1,".Vs.",sn2)
        multinames = paste( sn1, sn2, sep = ".Vs.")
        sampABstat = plotRes.multi.pdf(sampAB, sampName = paste( sn1, sn2, sep = ".Vs.") ,
                                       outputfolder = outputfolder,
                                       pdf = TRUE,
                                       main=mapply( function(x) x[1] , str_split(new_samples[1], "_")), 
                                       sn1 = paste0(new_samples[i], "mafa"), 
                                       sn2 = paste0(new_samples[j], "mafa"),
                                       minAF = 0.01,
                                       binw= binw, widthadj=0, heightadj=0
        )
        
        sampABsub = sampAB[ union(sampABstat$subArow, sampABstat$subBrow) %>% union(sampABstat$pubTrow ), ]
        
        sampABsub = sampABsub %>%
          filter(pubOrSub %in% c("public", 
                                 paste0("private=", c(sn1, sn2)),  
                                 sprintf("private=%s,%s", sn1, sn2),
                                 sprintf("private=%s,%s", sn2, sn1) )
          ) %>%
          mutate(pubOrSubadj = ifelse( grepl(pubOrSub, pattern = ",") , "shared" , pubOrSub) )
        
        #pdf( sprintf("MRS_analysis/%s.ccf.pdf", simple), width = 4, height =  3.5)
        pccf = plt.ccf.pairs(sampABsub, 
                             sn1 = paste0(new_samples[1], "mafa"), 
                             sn2 = paste0(new_samples[2], "mafa"),
                             pob = "pubOrSubadj"
        )
        multilist[[multinames]] = list( sampAB = sampABsub, sampABstat = sampABstat, pccf = pccf )
        
      }
    }
  }else{
    return(sampAB)
  }
  
  return(multilist)
  
}


######################################################

formate = "sequenza"

#set Path of Titan and sequenza.
if(formate != "sequenza"){
  #for tian Path
  titanPath="/data1/qingjian/Project/QiuMZ/TitanCNA/scripts/snakemake/results/titan/hmm/Seg/" 
  titanfiles = list.files(path = titanPath)
}else{
  #for sequenza Path
  titanPath="/data1/qingjian/Project/ThreePrimary/Results/sequenza/Seg" 
  titanfiles = list.files(path = titanPath)
}


singleRuns = list()

sampABlist = list()

for(cancer in names(sampleGroups)){

  message(cancer)

  samples = sampleGroups[[cancer]]
  
  #set purity info
  if(formate == "sequenza"){
    normalproportion = setNames( 1- purity_table[ purity_table$sample %in%  samples, "CCF.manual"] ,nm = samples) #Using sequenza purity info
  }else if(formate == "titan"){
    normalproportion = setNames( 1- purity_table[ purity_table$sample %in%  samples, "CCF.manual"] ,nm = samples)
  }
    
  #change mutation formate for analysis.
  GC01_mutations = GC01 %>%
    filter(Tumor_Sample_Barcode %in% samples) %>%
    mutate(
      Chromosome = str_c("chr", Chromosome),
      DP = Ref_allele_depth + Alt_allele_depth,
      NV = Alt_allele_depth) %>%
    pivot_wider(
     id_cols = c("mutid","Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Reference_Allele","Tumor_Seq_Allele2" ),
     names_from = c("Tumor_Sample_Barcode"),
     values_from = c("VAF", "DP","NV"),
     names_glue = "{.value}.{Tumor_Sample_Barcode}",
     values_fill  = 0
    )
    
  colnames(GC01_mutations)[1:9] = c("mutid", "Hugo_Symbol", "chr", "from", "to", "Variant_Classification", "Variant_Type", "ref", "alt")
  
  GC01_mutations = as.data.frame(GC01_mutations)
  ##Set minial Depth is 50 for absent mutations.
  
  clumns.dp = colnames(GC01_mutations)[startsWith(colnames(GC01_mutations), "DP")  ]
  
  for(i in clumns.dp ){
    GC01_mutations[, i] = ifelse( (GC01_mutations[,i]) == 0, 50, GC01_mutations[,i])
  }
  
    #for sequenza
    #singleRuns[[cancer]] =  runSampAB(GC01_mutations, samples, new_samples = samples, normalproportion, titanPath = titanPath, seg.formate = formate, t = 1)
    
  sampABlist[[cancer]] =  runSampAB(GC01_mutations, samples, new_samples = samples, normalproportion, titanPath = titanPath, seg.formate = formate, t = 1, plot.hist = FALSE)
  
}

############################################################

saveRDS(sampABlist, file = "data/sampABlist.rds")

