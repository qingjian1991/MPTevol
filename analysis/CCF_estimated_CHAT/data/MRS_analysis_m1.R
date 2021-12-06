# The main code is from <<Sun, R., et al., Between-region genetic divergence reflects the mode and tempo of tumor evolution. Nat Genet, 2017. 49(7): p. 1015-1024.>>


# plot multi sample mutation Table
script.dir <- dirname(sys.frame(1)$ofile)
library(caTools)
library(KernSmooth)
library(RColorBrewer)


#' trans2sampAB
#' trans non_tail_mutations into samAB
trans2sampAB = function(non_tail_mutations, samples, new_samples =NULL){
  
  sampAB = non_tail_mutations
  
  sampledb = c()
  if(is.null(new_samples)){
    sampledb = data.frame(samples = samples, samples1 = samples)
  }else{
    sampledb = data.frame(samples = samples, samples1 = new_samples)
  }
  
  #reshape data.
  for(i in 1:nrow(sampledb)){
    
    sn = sampledb$samples[i]
    sn1 = sampledb$samples1[i]
    
    if( sprintf("VAF.%s", sn) %in% colnames(sampAB)  ){
      sampAB[,paste(sn1, "mafc", sep = "")] = sampAB[, sprintf("VAF.%s", sn) ]
      sampAB[,paste(sn1, "refc", sep = "")] = sampAB[, sprintf("DP.%s", sn)] - sampAB[, sprintf("NV.%s", sn)]
      sampAB[,paste(sn1, "altc", sep = "")] = sampAB[, sprintf("NV.%s", sn)]
      sampAB[,paste(sn1, "d", sep = "")] = sampAB[, sprintf("DP.%s", sn)]
      sampAB[,paste(sn1, "mafa", sep = "")] = 0
      sampAB[,paste(sn1, "ccf", sep = "")] = 0
      sampAB[,paste(sn1, "ccfSD", sep = "")] = 0
    }else{
      message( sprintf("Data without samples: %s", sn) )
    }
    
  }
  sampAB[,"mergeMAFA"] = NA
  sampAB[,"mergeCCF"] = NA
  sampAB[,"mergeCCFsd"] = NA
  
  sampAB
}


#' sampAB2newNames
#' change sampAB into new names
sampAB2newNames = function(sampAB, samples, new_samples){
  #reshape data.
  for(i in 1:length(samples)){
    
    sn = samples[i]
    sn1 = new_samples[i]
    
    colms = c("mafc","refc","altc","d","mafa","ccf","ccfSD","pu","pa","sAGP","nt","nb","seg")
    
    if( sprintf("VAF.%s", sn) %in% colnames(sampAB)  ){
      
      for(j in colms){
        sampAB[,paste(sn1, j, sep = "")] = sampAB[,paste(sn, j, sep = "")]
        sampAB[,paste(sn, j, sep = "")] = NULL
      }
    }else{
      message( sprintf("Data without samples: %s", sn) )
    }
  }
  sampAB
}



#grep mutation table for a bunch of samples from the output mutation table (d is produced from a tsv file) produced by VAP
getSampMutMulti <- function(samples, normal, d, cmedianTh, original) {
  rindexTF = vector()
  cindex = vector()
  for (si in 1:length(samples)) {
    if (length(rindexTF) > 0) {
      rindexTF = rindexTF | grepl(paste(samples[si],"\\[",sep=""), d$somatic)
    } else {
      rindexTF = grepl(paste(samples[si],"\\[",sep=""), d$somatic)
    }
    cindex = append(cindex, match(paste(samples[si], "maf", sep=""), colnames(d)))
  }
  rindex = which(rindexTF)
  cindex = as.vector(sapply(cindex, function(x){c(x-1,x,x+1)}))
  
  ncindex = match(paste(normal, "maf", sep=""), colnames(d))
  ncindex = as.vector(sapply(ncindex, function(x){c(x,x+1)}))  #normal samples maf and depth indexes
  
  dronindex = match("dron", colnames(d))
  gnindex = match("geneName", colnames(d))
  glindex = match("geneLoc", colnames(d))
  gfindex = match("functionalClass", colnames(d))
  caddindex = match("CADD_phred", colnames(d))
  gerpindex = match("GERP_RS", colnames(d))
  siftindex = match("SIFT_score", colnames(d))
  polyphenindex = match("Polyphen2_HVAR_pred", colnames(d))
  somindex = match("somatic", colnames(d))
  
  if (! is.na(caddindex)) {
    res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,caddindex,gerpindex,siftindex,polyphenindex,dronindex,somindex,ncindex)]
  } else {
    res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,dronindex,somindex,ncindex)]
  }
  
  resColnames = colnames(res)
  resAdd = t(data.frame(apply(res, 1, function(x, cindex, original, resColnames) {     #x is every row
    maxMaf = 0
    maxTlod = 0
    ssb = 0                    #combined strand bias check
    ssbc = 0
    ssbp = 0
    refdav = 0
    altdav = 0
    for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample get maxMaf
      ss = strsplit(as.character(x[j+1]), "\\|")
      mafTmp = as.numeric(ss[[1]][1])
      oriTlodTmp = 0
      if ( grepl("\\|", as.character(x[j])) ) {
        ori = strsplit(as.character(x[j]), "\\|")
        oriLod = strsplit(ori[[1]][3], ",")            
        oriTlodTmp = as.numeric(oriLod[[1]][1])
      }
      if (mafTmp > maxMaf) {
        maxMaf = mafTmp    #define max maf
      }
      if (original > 0 & oriTlodTmp > maxTlod) {
        maxTlod = oriTlodTmp
      }
    }
    resVector = vector()
    for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample
      sampleName = resColnames[j]
      sampleNameMaf = paste(sampleName, "mafc", sep="")
      sampleNameMafa = paste(sampleName, "mafa", sep="")
      sampleNameCcf = paste(sampleName, "ccf", sep="")
      sampleNameCcfSD = paste(sampleName, "ccfSD", sep="")
      sampleNameAlt = paste(sampleName, "altc", sep="")
      sampleNameRef = paste(sampleName, "refc", sep="")
      oriTlod = 0
      if ( grepl("\\|", as.character(x[j])) ) {
        ori = strsplit(as.character(x[j]), "\\|")
        oriLod = strsplit(ori[[1]][3], ",")
        oriTlod = as.numeric(oriLod[[1]][1])
      }
      ss = strsplit(as.character(x[j+1]), "\\|")
      mafTmp = as.numeric(ss[[1]][1])
      cmeme = strsplit(ss[[1]][3], ",")
      endBias = as.numeric(ss[[1]][2])
      strandBiases = strsplit(ss[[1]][4], ",")
      strandBias = as.numeric(strandBiases[[1]][1])
      strandBiasRef = 0
      strandBiasFisherP = -1
      if (length(strandBiases[[1]]) > 1){
        strandBiasRef = as.numeric(strandBiases[[1]][2])
        strandBiasFisherP = as.numeric(strandBiases[[1]][3])
      }
      mappingBias = as.numeric(ss[[1]][5])
      cmedianav = as.numeric(cmeme[[1]][2])
      cmemeSum = sum(as.numeric(cmeme[[1]]))
      
      #decide a b allele count
      refnow = round(as.numeric(x[j+2])*(1-mafTmp))
      altnow = round(as.numeric(x[j+2])*mafTmp)
      
      if ( mafTmp > 0 ) {
        ssb = ssb + strandBias*altnow
        ssbp = ssbp + strandBiasFisherP*altnow
        refdav = refdav + refnow
        altdav = altdav + altnow
        ssbc = ssbc + 1
      }
      
      #decide mafNow
      if ( mafTmp == 0 ) {
        mafNow = mafTmp
      } else if (endBias < 0.9 & ((strandBias != 0 & strandBias != 1) | (strandBiasFisherP > 0.7 & refnow >= 10 & altnow >= 5 & mafTmp >= 0.1)) & mappingBias < 0.8 & cmemeSum < 5.2 & cmedianav < cmedianTh) {  
        mafNow = mafTmp
        if (original > 0) {   #if oriLod
          if (oriTlod < original) {
            if (!(maxMaf > 0.1 & (original == 0 | (original > 0 & maxTlod >= original)))) {   #not called
              mafNow = 0
            }
          }
        }
      } else {
        if (maxMaf > 0.2 & (original == 0 | (original > 0 & maxTlod >= original)) & mafTmp >= 0.01) {
          if (cmedianav < 4) {
            mafNow = mafTmp
          } else {
            mafNow = 0
          }
        } else if (maxMaf > 0.05 & (original == 0 | (original > 0 & maxTlod < original)) & mafTmp >= 0.04) {
          if (cmedianav == 1 & mappingBias == 0) {
            mafNow = mafTmp
          } else {
            mafNow = 0
          }
        } else {
          mafNow = 0
        }
      }
      resVector = c(resVector, c(mafNow, 0, 0, 0, refnow, altnow))
      names(resVector)[(length(resVector)-5):length(resVector)] = c(sampleNameMaf,sampleNameMafa,sampleNameCcf,sampleNameCcfSD,sampleNameRef,sampleNameAlt)
    } #for each sample
    
    ssb = ssb/altdav
    ssbp = ssbp/altdav
    altdav = altdav/ssbc
    refdav = refdav/ssbc
    if (((ssb >= 0.95 | ssb <= 0.05) & ssbp < 0.1)| ssbp < 0.001) {                          #multiple sample strand bias
      #if ((ssb >= 0.95 | ssb <= 0.05)| ssbp < 0.001) {                                        #multiple sample strand bias 
      for (j in seq(6,(3+length(cindex)), by=3)) {
        sampleName = resColnames[j]
        sampleNameMaf = paste(sampleName, "mafc", sep="")
        resVector[sampleNameMaf] = 0
      }
    }
    
    totalMaf = sum(as.numeric(resVector[grepl("mafc", names(resVector))]))
    totalAlt = sum(as.numeric(resVector[grepl("altc", names(resVector))]))
    totalRef = sum(as.numeric(resVector[grepl("refc", names(resVector))]))
    mergeMAFC = round(totalAlt/(totalAlt+totalRef), 4)
    mergeMAFA = mergeMAFC
    mergedCCF = mergeMAFC
    mergedCCFsd = mergeMAFC
    resVector = c(resVector, totalMaf, totalAlt, totalRef, mergeMAFC, mergeMAFA, mergedCCF, mergedCCFsd)
    names(resVector)[(length(resVector)-6):length(resVector)] = c("totalMaf", "totalAlt", "totalRef", "mergeMAFC", "mergeMAFA", "mergeCCF", "mergeCCFsd")
    
    resVector                #return the result
    
  }, cindex=cindex, original=original, resColnames=resColnames)))
  
  
  
  res = cbind(res, resAdd)
  res = res[which(res$totalMaf != 0),]
  return(res)
}


#adjust CCF titan for multi samples
adjust.ccf.titan.multi <- function(sampAB, samples, normalproportion , t, titanPath="./titan/",  seg.formate = "titan" , correctColname=FALSE, overadj=1.6, sigTh=0.9) {
  
  titanfiles = list.files(path = titanPath)
  
  purities = vector()
  
  for (i in 1:length(samples)) {
    sn = samples[i]
    message(sn)
    
    #cnv.inputA  = paste(titanPath, sn, "_cluster1.titan.ichor.seg.txt", sep="")
    cnv.inputA  = paste(titanPath, titanfiles[ grepl(pattern = sn, x = titanfiles)]  , sep="/")
    message(cnv.inputA)
    cnvA = read.delim(cnv.inputA)
    cnvA$normalproportion = normalproportion[sn] #********************
    
    #formate is sequenza, not titan.
    if(seg.formate != "titan"){
      
      cnvA =  cnvA %>%
        dplyr::rename(
          Chromosome = chromosome,
          Start = start.pos,
          End = end.pos,
          Copy_Number = CNt,
          MinorCN = B,
          
        ) %>%
        mutate(Cellular_Prevalence = 1
               )
      
    }
    
    #cnvA$Cellular_Prevalence = ifelse(is.na(cnvA$Cellular_Prevalence), 1, cnvA$Cellular_Prevalence)
    
    cnvA = cnvA[which(!is.na(cnvA$Cellular_Prevalence)),]                #skip NA
    
    
    #cnvA = cnvA[which(cnvA$num.mark > 9),]                              #skip two few marks
    cnvA$nt = cnvA$Copy_Number
    if ("MinorCN" %in% colnames(cnvA) ) {
      cnvA$nb = cnvA$MinorCN
    } else {
      cnvA$nb = partialRound(cnvA$copynumber*(
        1-((cnvA$allelicratio - cnvA$normalproportion*0.5)/(1-cnvA$normalproportion) - (1-cnvA$cellularprevalence)*0.5)
        /cnvA$cellularprevalence))
    }
    
    cnvSeqNames = cnvA$Chromosome
    if (!grepl("chr",cnvA$Chromosome[1])){
      cnvSeqNames = paste("chr",cnvA$Chromosome,sep="")
    }
    snvSeqNames = sampAB$chr
    if (!grepl("chr",sampAB$chr[1])){
      snvSeqNames = paste("chr",sampAB$chr,sep="")
    }
    cnvRangeA = GRanges(seqnames = cnvSeqNames, ranges = IRanges(cnvA$Start, end=cnvA$End), strand=rep('+',dim(cnvA)[1]))
    snvRange = GRanges(seqnames = snvSeqNames, ranges = IRanges(sampAB$from, end=sampAB$to), strand=rep('+',dim(sampAB)[1]))
    foA = findOverlaps(snvRange, cnvRangeA)
    
    #pa,pu,nt,nb,seg
    message("info table building")
    queHits = queryHits(foA)
    subHits = subjectHits(foA)
    infos = t(sapply(1:dim(sampAB)[1], function(x, queHits, subHits) {                           
      if (x %in% queHits){
        queMIndex = match(x, queHits)
        subMIndex = subHits[queMIndex]
        c(cnvA$Cellular_Prevalence[subMIndex],     #pa
          cnvA$nt[subMIndex],                     #nt
          cnvA$nb[subMIndex],                     #nb
          subMIndex)                              #seg
      } else {
        c(0,2,1,0)
      }}, queHits = queHits, subHits = subHits))
    infos = data.frame(infos)
    message("info table built")
    
    pa1 = infos[,1]
    nt1 = infos[,2]
    nb1 = infos[,3]
    seg1 = infos[,4]
    con1 = as.numeric(cnvA$normalproportion[1])          #pu1 = 1-cnvA$normalproportion[1]
    pu1 = 1 - (con1 + (1-max(as.numeric(pa1)))*(1-con1))
    
    if(pu1 == 0){
      pu1 = 1-cnvA$normalproportion[1]
    }
    
    sAGP = pa1*(1-con1)
    message(pu1)
    purities = append(purities, pu1)
    names(purities)[length(purities)] = sn
    
    #sampAB = data.frame(sampAB, pu=pu1, pa=pa1, sAGP=sAGP, nt=nt1, nb=nb1, seg=seg1)
    sampAB = cbind(sampAB, data.frame(pu=pu1, pa=pa1, sAGP=sAGP, nt=nt1, nb=nb1, seg=seg1) )
    colnames(sampAB)[(dim(sampAB)[2]-5):dim(sampAB)[2]] = paste(sn,colnames(sampAB)[(dim(sampAB)[2]-5):dim(sampAB)[2]], sep="")
    if ( correctColname == TRUE ) {
      colnames(sampAB) = gsub("\\.","-",colnames(sampAB))
    }
  }
  
  
  for( i in 1:dim(sampAB)[1]) {  # rescale the maf and calculate CCF
    
    if (i %% 1000 == 0) {
      message(i)
    }
    
    foundSites = 0             # count how many sites found
    depthTotal = 0
    mafaTotal = 0
    ccfTotal = 0
    ccfsdTotal = 0
    for (j in 1:length(samples)) {
      sn = samples[j]
      maf1 = as.numeric(sampAB[i, match(paste( sn, "mafc", sep=""), colnames(sampAB))])
      if (maf1 > t) {
        foundSites = foundSites+1
      }
    }
    
    for (j in 1:length(samples)) {
      sn = samples[j]
      
      pa1 = as.numeric(sampAB[i, match(paste(sn, "pa", sep=""), colnames(sampAB))])
      nt1 = as.numeric(sampAB[i, match(paste(sn, "nt", sep=""), colnames(sampAB))])
      nb1 = as.numeric(sampAB[i, match(paste(sn, "nb", sep=""), colnames(sampAB))])
      maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
      refc1 = as.numeric(sampAB[i, match(paste(sn, "refc", sep=""), colnames(sampAB))])
      altc1 = as.numeric(sampAB[i, match(paste(sn, "altc", sep=""), colnames(sampAB))])
      pu1 = as.numeric(sampAB[i, match(paste(sn, "pu", sep=""), colnames(sampAB))])       #cell purity
      #if (nt1 > 0) pu1 = nt1*pu1/(nt1*pu1+2*(1-pu1))                                     #effective purity
      sAGP = as.numeric(sampAB[i, match(paste(sn, "sAGP", sep=""), colnames(sampAB))])    #segmental aneu- ploid genome proportion
      if (nt1 > 0 & nt1 != 2) pu1 = nt1*sAGP/(nt1*sAGP+2*(1-sAGP))                                   #effective purity
      
      if (maf1 > 0) {
        if (maf1 > t & foundSites >= 2) {
          CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"unknown",overadj=overadj,sigTh=sigTh)
          sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
          sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
          sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
        } else {
          CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"late",overadj=overadj,sigTh=sigTh)
          sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
          sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
          sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
        }
      }
      
      #for merged MAF
      dep1 = as.numeric(sampAB[i, match(paste(sn, "d", sep=""), colnames(sampAB))])
      mafa1 = as.numeric(sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))])
      ccf1 = as.numeric(sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))])
      ccfsd1 = as.numeric(sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))])
      depthTotal = depthTotal + dep1
      mafaTotal = mafaTotal + dep1*mafa1
      ccfTotal = ccfTotal + dep1*ccf1
      ccfsdTotal = ccfsdTotal + dep1*ccfsd1
      #for merged MAF
      
    }     #for each sample
    sampAB[i, match("mergeMAFA", colnames(sampAB))] = round(mafaTotal/depthTotal, 5)
    sampAB[i, match("mergeCCF", colnames(sampAB))] = round(ccfTotal/depthTotal, 5)
    sampAB[i, match("mergeCCFsd", colnames(sampAB))] = round(ccfsdTotal/depthTotal, 5)
  }
  return(sampAB)
}


#plotting pairwise comparisons, and return ITH stats.
plotRes.multi.pdf <- function(sampAB, sampName, main=sampName, sn1n="", sn2n="", sn1, sn2, minAF, ratio=1, plotAF=TRUE, pdf=TRUE,
                              outputfolder = "./",
                              alpha=1, binw=0, widthadj=0, heightadj=0, nohistlegend = FALSE, ssAF=0, pob="pubOrSub") {
  
  sn1s = gsub("mafc","",sn1)
  sn1s = gsub("mafa","",sn1s)
  sn1s = gsub("ccf","",sn1s)
  sn2s = gsub("mafc","",sn2)
  sn2s = gsub("mafa","",sn2s)
  sn2s = gsub("ccf","",sn2s)
  
  nb1i = match(paste(sn1s, "nb", sep=""),colnames(sampAB))
  nb2i = match(paste(sn2s, "nb", sep=""),colnames(sampAB))
  dp1i = match(paste(sn1s, "d", sep=""),colnames(sampAB))
  dp2i = match(paste(sn2s, "d", sep=""),colnames(sampAB))
  
  #sampAB = sampAB[which(grepl(paste(sn1s,"\\[",sep=""), sampAB$somatic) | grepl(paste(sn2s,"\\[",sep=""), sampAB$somatic)),]
  sampAB = sampAB[which((sampAB[,nb1i] > 0 & sampAB[,nb2i] > 0) | (sampAB[,nb1i] == 0 & sampAB[,nb2i] == 0)),]                 #remove different LOH locations
  
  maf1Index = match(sn1, colnames(sampAB))
  maf2Index = match(sn2, colnames(sampAB))
  mafa1Index = match(paste(sn1s,"mafa",sep=""), colnames(sampAB))    #for mafa
  mafa2Index = match(paste(sn2s,"mafa",sep=""), colnames(sampAB))    #for mafa
  
  #For recording AUC
  AUCstuff = subclonalMut(sampAB, sn1s, sn2s, minAF, ssAF=ssAF, pob=pob)
  
  #check depth power to reject a presence of a mutation
  #depthPowerKeep <- as.vector(apply(sampAB[15,], 1, function(x,mafa1i,mafa2i,dp1i,dp2i) {
  #  if( (as.numeric(x[mafa1i]) == 0) & (as.numeric(x[mafa2i]) == 0) ){
  #    return(FALSE)
  # }
  #  if(as.numeric(x[mafa1i]) == 0){
  #    vaf = as.numeric(x[mafa2i])
  #    if (vaf > 1 | vaf < 0){FALSE}  else if (pbinom(0,as.integer(x[dp1i]),vaf) < 0.05){TRUE} else {FALSE}}
  #  else if (as.numeric(x[mafa2i]) == 0){vaf = as.numeric(x[mafa1i])
  #  if (vaf > 1 | vaf < 0){FALSE}  else if (pbinom(0,as.integer(x[dp2i]),vaf) < 0.05){TRUE} else {FALSE}}
  # else {TRUE}
  #}, mafa1i=mafa1Index,mafa2i=mafa2Index,dp1i=dp1i,dp2i=dp2i))
  
  
  
  #check depth power to reject a presence of a mutation
  depthPowerKeep <- as.vector(
    apply( sampAB, 1, 
           function(x,mafa1i,mafa2i,dp1i,dp2i) {
             if( (as.numeric(x[mafa1i]) == 0) & (as.numeric(x[mafa2i]) == 0) ){
               return(FALSE)
             } 
             
    if(as.numeric(x[mafa1i]) == 0){
      vaf = as.numeric(x[mafa2i])
      if (vaf > 1 | vaf < 0){
        FALSE
      }else if(pbinom(0,as.integer(x[dp1i]),vaf) < 0.05){
        TRUE
      }else{
        FALSE
      }
    }else if (as.numeric(x[mafa2i]) == 0){
      vaf = as.numeric(x[mafa1i])
      if (vaf > 1 | vaf < 0){
        FALSE
      }else if(pbinom(0,as.integer(x[dp2i]),vaf) < 0.05){
        TRUE
      }else{
        FALSE
      }
    }else{
      TRUE
    }       
             
    }, mafa1i=mafa1Index,mafa2i=mafa2Index,dp1i=dp1i,dp2i=dp2i))

  sampAB = sampAB[depthPowerKeep,]
  
  subMuts = subclonalMut(sampAB, sn1s, sn2s, minAF, ssAF=ssAF, pob=pob)       #subclonal mutations
  subMuts$rAUC = AUCstuff$rAUC
  subMuts$weightAF = AUCstuff$weightAF
  subMuts$subArow = AUCstuff$subArow
  subMuts$subBrow = AUCstuff$subBrow
  
  subMuts$sn1s = sn1s
  subMuts$sn2s = sn2s
  
  #allA_Rows = which(sampAB[,maf1Index] > minAF & sampAB[,mafa1Index] > minAF & sampAB[,maf1Index] <= 1)
  allA_Rows = intersect(union(subMuts$subAi, subMuts$pubTi), which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))
  #message(paste(subMuts$subArow[match(allA_Rows[which(sampAB[allA_Rows, maf1Index] <= 0.05)], subMuts$subAi)], collapse="  "))
  subA_Rows = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))
  ssA_Rows  = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1 & sampAB[,maf2Index] <= ssAF))
  BinWidthA = round(dpih(sampAB[allA_Rows, maf1Index]/ratio),2)
  if ((BinWidthA < 0.02 & length(allA_Rows) < 3000) | BinWidthA == 0) { BinWidthA = 0.02 }
  
  #allB_Rows = which(sampAB[,maf2Index] > minAF & sampAB[,mafa2Index] > minAF & sampAB[,maf2Index] <= 1)
  allB_Rows = intersect(union(subMuts$subBi, subMuts$pubTi), which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))
  allAB_Rows = union(allA_Rows,allB_Rows)
  subB_Rows = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))
  ssB_Rows  = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1 & sampAB[,maf1Index] <= ssAF))
  
  BinWidthB = round(dpih(sampAB[allB_Rows, maf2Index]/ratio),2)
  if ((BinWidthB < 0.02 & length(allB_Rows) < 3000) | BinWidthB == 0) { BinWidthB = 0.02 }
  BinWidth = min(c(BinWidthA, BinWidthB, 0.1))
  if ( binw != 0 ) {
    BinWidth = binw
  }
  message(paste("bin width: ", BinWidthA, BinWidthB, BinWidth, sep=" "))
  nbreaksA = round((max(sampAB[allA_Rows, maf1Index]/ratio)-min(sampAB[allA_Rows, maf1Index]/ratio)+0.01)/BinWidth)
  nbreaksB = round((max(sampAB[allB_Rows, maf2Index]/ratio)-min(sampAB[allA_Rows, maf2Index]/ratio)+0.01)/BinWidth)
  nbreaks = ceiling(diff(range(minAF,1))/BinWidth)
  message(paste("nbreaks: ", nbreaksA, nbreaksB, nbreaks, sep=" "))
  breaksA = seq(minAF,1,length.out=nbreaks)
  breaksB = seq(minAF,1,length.out=nbreaks)
  
  
  sampAh = hist(sampAB[allA_Rows, maf1Index]/ratio, breaks=breaksA,  plot=F)
  sampAhsub = hist(sampAB[subA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  sampAhss = hist(sampAB[ssA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  ylimup = max(sampAh$count)
  
  sampBh = hist(sampAB[allB_Rows, maf2Index]/ratio, breaks=breaksB, plot=F)
  sampBh$counts = sampBh$counts*(-1)
  sampBhsub = hist(sampAB[subB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhsub$counts = sampBhsub$count*(-1)
  sampBhss = hist(sampAB[ssB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhss$counts = sampBhss$counts*(-1)
  ylimdown = min(sampBh$count)
  
  #set ylim up and down equal to the same scale
  if (abs(ylimup) >= abs(ylimdown)) {
    ylimdown = (-1)*ylimup
  } else {
    ylimup = (-1)*ylimdown
  }
  
  
  ssfr = seq(0.05,0.24,by=0.01)
  ssfs = as.vector(sapply(ssfr, function(x, sAB, maf1I, maf2I){ length(which((sAB[,maf1I] >= x & sAB[, maf1I] < x+0.01 & sAB[, maf2I] == 0) | (sAB[,maf2I] >= x & sAB[, maf2I] < x+0.01 & sAB[, maf1I] == 0)))/
      length(which((sAB[,maf1I] >= x & sAB[, maf1I] < x+0.01) | (sAB[,maf2I] >= x & sAB[, maf2I] < x+0.01)))}, sAB = sampAB, maf1I=maf1Index, maf2I=maf2Index))
  
  if (plotAF == TRUE) {
    
    if (pdf == TRUE) {
      pdf(file = paste(outputfolder, sampName, ".hist.pdf", sep=""), width = 8+widthadj, height = 8+heightadj, useDingbats=FALSE)
    }
    par(mar=c(4.5,5,4.5,0))
    
    plot( sampAh, col=rgb(0,0,0,1/4), xlim=c(0, 1), ylim=c(ylimdown,ylimup), border=F, ylab="# of Mutations", xlab="Allele Frequency", axes = F, main = main, cex.lab = 2.3, cex.main = 2.3)  # first histogram
    plot( sampAhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green set border
    plot( sampAhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
    
    plot( sampBh,col=rgb(0,0,0,1/4), border=F, add=T )  # second histogram
    plot( sampBhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green
    plot( sampBhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
    sampAhss$counts = 0                                  #make black line
    sampBhss$counts = 0                                  #make black line
    plot( sampAhss, col="black", add=T, border=F)        #make black line
    plot( sampBhss, col="black", add=T, border=F)        #make black line
    
    #sn1s = sn1n
    if (sn1s == ""){
      sn1s = gsub(".+(Core\\d+)","\\1",sn1, perl=T, ignore.case=T)
      sn1s = gsub(".+(Cor\\d+)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Sec(\\d+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Primary\\_(\\d+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Primary\\_(\\w+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Met\\_(\\d+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Met\\_(\\w+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub("Primary","Pri",sn1s, perl=T, ignore.case=T)
      sn1s = gsub("mafc","",sn1s)
      sn1s = gsub("mafa","",sn1s)
      sn1s = gsub("CRCTumor","",sn1s)
      sn1s = gsub("HCT116_","",sn1s)
    }
    #sn2s = sn2n
    if (sn2s == ""){
      sn2s = gsub(".+(Core\\d+)","\\1",sn2, perl=T, ignore.case=T)
      sn2s = gsub(".+(Cor\\d+)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Sec(\\d+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Primary\\_(\\d+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Primary\\_(\\w+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Met\\_(\\d+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Met\\_(\\w+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub("Primary","Pri",sn2s, perl=T, ignore.case=T)
      sn2s = gsub("mafc","",sn2s)
      sn2s = gsub("mafa","",sn2s)
      sn2s = gsub("CRCTumor","",sn2s)
      sn2s = gsub("HCT116_","",sn2s)
    }
    message(sn1s)
    message(sn2s)
    
    axis(side=1,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=1.7)
    axis(side=2,at=decideTickAt(ylimdown, ylimup),labels=allAbs(decideTickAt(ylimdown, ylimup)),cex.axis=1.7)
    text(x=0.5,y=ylimdown,labels=paste(sn2s,",",sum(subMuts$pubTn,subMuts$sharedBn,subMuts$ssBn), "SSNVs", sep = " "), cex=2.2)
    text(x=0.5,y=ylimup,labels=paste(sn1s,",",sum(subMuts$pubTn,subMuts$sharedAn,subMuts$ssAn), "SSNVs", sep = " "), cex=2.2)
    
    #stats starting from here!
    fHsub = round(mean(c(subMuts$ratioHighSubA,subMuts$ratioHighSubB)),3)
    fst = round(subMuts$FST,3)
    ksd = round(subMuts$KSD,3)
    text(x=0.8,y=3*ylimup/4, labels=bquote(paste("fH"["sub"], " = ", .(fHsub))),cex=2.2)
    text(x=0.8,y=(3/4-0.167)*ylimup, labels=bquote(paste("FST", " = ", .(fst))),cex=2.2)
    text(x=0.8,y=(3/4-0.334)*ylimup, labels=bquote(paste("KSD", " = ", .(ksd))),cex=2.2)
    
    npub = subMuts$pubTn
    legendText = c(paste("Public ","(",npub,")",sep=""),"Pvt-Shared","Pvt-Rgn Specific")
    legendXpos = 0.57
    legendYpos = ylimdown/3
    lengedCol = c(rgb(0,0,0,1/4),rgb(178/255,223/255,138/255,1),rgb(31/255,120/255,180/255,1))
    if (nohistlegend == TRUE){
      legendText = c(paste("",npub,sep=""))
      legendXpos = 0.75
      legendYpos = ylimdown/2
      lengedCol = c(rgb(0,0,0,1/4))
    }
    legend(legendXpos,ylimdown/3, legend=legendText, col=lengedCol, pch=15, bty="n", cex=1.7)
    if (pdf == TRUE) {
      dev.off()
    }
  }
  
  return(subMuts)
}


plt.ccf.pairs = function(sampAB, sn1 = "", sn2 = "", pob="pubOrSubadj"){
  sampAB %>%
    ggplot(aes_string( x = sn1,  y = sn2, col = pob ) ) + theme_classic() +
    geom_jitter(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
    #xlim(0, 1) +
    #ylim(0, 1) +
    scale_color_manual( values = c(rgb(31/255,120/255,180/255,1),
                                   "#A2A0C8",
                                   rgb(0,0,0,1/4),
                                   rgb(178/255,223/255,138/255,1)
                                   )  ) +
    theme(legend.position = "none")
}


#' rNeutral
#' Giving a neutral test. The input is same with rAUC

rNeutral <- function(data, lower, mafs, depths) {
  
  if (length(mafs) <= 8) {
    lower = round((0.08/length(mafs)),2)
  } else {
    lower = max(0.005, round((0.08/length(mafs)),3))
  }
  #message(paste("lower: ", lower,sep=""))
  weightAF = weightAFs(data, mafs, depths)
  #message(paste(weightAF, collapse=" "))
  weightAF = weightAF[which(weightAF >= lower & weightAF <= 0.25)]
  #message(paste("weightedMuts: ", length(weightAF), sep=""))
  vafs = seq(lower,0.25,by=0.01)
  nstep = length(vafs)
  counts = vector()
  ncounts = ((1/vafs)-(1/0.25))/((1/lower)-(1/0.25))
  
  for ( i in 1:length(vafs) ) {
    counts = append(counts, length(which(weightAF > vafs[i]))-length(which(weightAF > vafs[nstep])))
  }
  counts = counts/counts[1]
  
  #bei = bezierCurve(vafs, counts, length(vafs))
  #AUC = trapz(bei$x,bei$y)
  
  #modify
  AUC = trapz(vafs,counts)
  
  #AUC = trapz(mafs,counts)
  nAUC = trapz(vafs,ncounts)
  rAUC = AUC/nAUC
  rAUC = round(rAUC,8)
  
  #neutral test using neutralitytest
  neuTout = neutralitytestr::neutralitytest(weightAF, fmin = lower, fmax = 0.25)
  #plot_all(neuTout) %>% print()
  
  summary(neuTout)
  
  p1 = data.frame(vafs = rep(vafs, 2), ncounts = c(ncounts, counts), type = rep(c("T","R"), each = length(vafs)) ) %>%
    ggplot(aes(x = vafs, y = ncounts, col = type, size = type)) + theme_classic() + 
    geom_point() +
    geom_line() +
    annotate( geom = "text",  x = 0.20, y =0.85, label = sprintf("rAUC: %.2f", rAUC) , size = 5 ) +
    scale_color_manual( values = c("red3","grey80") ) +
    scale_size_manual(values = c(0.8, 2)) +
    labs(x = "f (pooled VAF)", y = "Fraction SSNVs in [f, fmax]") +
    theme( legend.position = "none" )
  
  
  return(list(neuTout = neuTout, p1 = p1))
  
}


subclonalMut <- function(sampAB, snA, snB, minAF=0.08, statsAF=0.08, highAF=0.2, ratio=1, crv=2.58, ssAF=0, pob="pubOrSub")  {                   #determinine subclonal mutations
  ccfAi = match(paste(snA, "ccf", sep=""), colnames(sampAB))
  ccfBi = match(paste(snB, "ccf", sep=""), colnames(sampAB))
  ccfsdAi = match(paste(snA, "ccfSD", sep=""), colnames(sampAB))
  ccfsdBi = match(paste(snB, "ccfSD", sep=""), colnames(sampAB))
  mafaAi = match(paste(snA, "mafa", sep=""), colnames(sampAB))
  mafaBi = match(paste(snB, "mafa", sep=""), colnames(sampAB))
  nbAi = match(paste(snA, "nb", sep=""), colnames(sampAB))
  nbBi = match(paste(snB, "nb", sep=""), colnames(sampAB))
  depthAi = match(paste(snA, "d", sep=""), colnames(sampAB))
  depthBi = match(paste(snB, "d", sep=""), colnames(sampAB))
  pobi = match(pob, colnames(sampAB))
  
  gNIndex = match("geneName",colnames(sampAB))
  gLIndex = match("geneLoc",colnames(sampAB))
  fCIndex = match("functionalClass",colnames(sampAB))
  
  rownames(sampAB) = 1:nrow(sampAB)
  
  # subclonal mutations
  
  #subAi: subclonal muttaions of A ( including private = A or private = A,B )
  subAi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaAi] > minAF & ((sampAB[,mafaBi] == 0 & (sampAB[,nbBi] != 0 | sampAB[,nbAi] == 0)) | sampAB[,mafaBi] != 0) )
  subArow = rownames(sampAB)[subAi]
  mutsA = sampAB[subAi,mafaAi]/ratio #mafa: adjust maf
  #ssAi: subclonal specific of A (limit to private = A), Region-specific subclonal SSNVs.
  ssAi  = intersect(subAi, which( sampAB[,mafaAi] > minAF & sampAB[,mafaBi] <= ssAF ))
  
  
  subBi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaBi] > minAF & ((sampAB[,mafaAi] == 0 & (sampAB[,nbAi] != 0 | sampAB[,nbBi] == 0)) | sampAB[,mafaAi] != 0) )
  subBrow = rownames(sampAB)[subBi]
  mutsB = sampAB[subBi,mafaBi]/ratio
  ssBi  = intersect(subBi, which( sampAB[,mafaBi] > minAF & sampAB[,mafaAi] <= ssAF ))
  
  KSD = as.numeric(ks.test( mutsA[which(mutsA > statsAF)], mutsB[which(mutsB > statsAF)] )$statistic)
  
  #for rAUC
  allSubRows = union(subAi,subBi) #Union of private=A or private=B or private=A,B
  mafs = paste(c(snA,snB), "mafa", sep="")
  depths = paste(c(snA,snB), "d", sep="")
  rAUCout = rAUC(sampAB[allSubRows,], 0.04, mafs, depths)     #need to add for lower bound
  rAUC = round(rAUCout$rAUC,8)
  weightAF = rAUCout$weightAF
  
  rNeutralout = rNeutral(sampAB[allSubRows,], 0.04, mafs, depths)
  
  
  # for PUB
  pubTi = which( sampAB[,pobi] == "public" &
                   !((sampAB[,mafaAi] < 0.15 & sampAB[,nbAi] == 0) | (sampAB[,mafaBi] < 0.15 & sampAB[,nbBi] == 0)) )
  pubTrow = rownames(sampAB)[pubTi]
  
  # for FST
  mutsSub = sampAB[allSubRows,]
  mutsSub = data.frame( maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi] )
  FST = fst.hudson(mutsSub, minAF=statsAF)
  
  
  # for other stats
  #subclona A mutations. Including A specific subclonal and A-B shared subclonal mutations
  mutsA2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF )), mafaAi]
  mutsAh2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF )), mafaAi]
  mutsASr2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] > 0.02)), mafaAi]   #shared between A and B
  
  #removed shared-B subclonal mutations. Regional-specific muttaions
  mutsASp2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] == 0)), mafaAi]
  mutsASph2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF & sampAB[,mafaBi] == 0)), mafaAi]
  
  mutsB2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF )), mafaBi]
  mutsBh2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF )), mafaBi]
  mutsBSr2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] > 0.02)), mafaBi]  #shared
  
  mutsBSp2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] == 0)), mafaBi]
  mutsBSph2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF & sampAB[,mafaAi] == 0)), mafaBi]
  
  
  # mutational counts
  pubTn = length(pubTi) #public mutations
  
  lenSubA=length(mutsA2) # total subclonal mutations
  lenSubB=length(mutsB2) #
  
  ssAn = length(ssAi)   #Region-specific subclonal mutations.
  ssBn = length(ssBi)
  sharedAn = length(mutsASr2) #shared sub-clonal mutations.
  sharedBn = length(mutsBSr2)
  #sharedn = length(union(subAi, subBi))-length(union(ssAi,ssBi))
  
  # list for output    
  muts = list(A=mutsA,B=mutsB,subAi=subAi,subBi=subBi, ssAi=ssAi, ssBi=ssBi, pubTi=pubTi, subArow=subArow, subBrow=subBrow, pubTrow = pubTrow,
              rNeutralout = rNeutralout,
              lenSubA=length(mutsA2),lenSubAh=length(mutsAh2),ratioHighSubA=length(mutsAh2)/length(mutsA2),
              lenSubB=length(mutsB2),lenSubBh=length(mutsBh2),ratioHighSubB=length(mutsBh2)/length(mutsB2),
              lenSsA=length(mutsASp2),lenHighSsA=length(mutsASph2),ratioHighSsA=length(mutsASph2)/length(mutsASp2),pSsA=length(mutsASp2)/length(mutsA2),
              lenSsB=length(mutsBSp2),lenHighSsB=length(mutsBSph2),ratioHighSsB=length(mutsBSph2)/length(mutsBSp2),pSsB=length(mutsBSp2)/length(mutsB2),
              
              lenSharedA=length(mutsASr2), lenSharedB=length(mutsBSr2), 
              
              ratioSharedA=length(mutsASr2)/length(mutsA2), 
              ratioSharedB=length(mutsBSr2)/length(mutsB2),
              
              FST=FST, KSD=KSD, rAUC=rAUC, weightAF=weightAF, pubTn=pubTn, sharedAn=sharedAn, sharedBn=sharedBn, ssAn=ssAn, ssBn=ssBn)
  return(muts)
}


subSSColor <- function(allindex, subindex, ssindex) {
  subsscolor = as.vector(sapply(allindex, function(x, subindex, ssindex){
    if (x %in% ssindex) {
      rgb(31/255,120/255,180/255,1)
    } else if (x %in% subindex) {
      rgb(178/255,223/255,138/255,1)
    } else {
      rgb(0,0,0,1/4)
    }},subindex=subindex, ssindex=ssindex))
  return(subsscolor)
}


nearestDecimal <- function(x) {
  r = x %% 10
  nearD = 0
  if (r > 5) {
    nearD = x + (10-r)
  } else {
    nearD = x - r
  }
  return(nearD)
}


decideTickAt <- function(ylimdown, ylimup) {  #assuming they are equal in abs
  
  if (abs(ylimdown) != abs(ylimup)){
    stop("ylimdown and up should be equal in abs!")
  }
  
  realmax = nearestDecimal(ylimup)
  incre = nearestDecimal(realmax/5)
  if (incre == 0) {
    incre = 5
  }
  njump = floor(realmax/incre)
  tailgap = realmax %% incre
  
  ylimup = realmax
  ylimdown = (-1)*realmax
  tickAt = c(ylimdown, seq(ylimdown + tailgap, 0, by=incre), seq(0, ylimup - tailgap, by=incre), ylimup)
  return(tickAt)
  
}


allAbs <- function(x) {
  for (i in 1:length(x)) {
    x[i] = abs(x[i])
  }
  return(x)
}


# this is a modified function modified from Clonal Heterogeneity Analysis Tool (CHAT)
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0473-4
#'@param A: Alt number of reads
#'@param S:
#'@param pa: P(CNA), cellular  prevalence of the CNA.
#'@param nt: number of total CNAs;
#'@param nb: number of minor CNAs.
#'
computeCCF <- function(f, A, S, pu, pa, sAGP, nt, nb, prior="unknown", overadj=1.6, sigTh=0.90) {
  
  ccf = 0
  ccf2 = 0
  sd = 0
  cc <- seq(0.02, 1, by = 0.01)
  evoType = "A1/A2/B/C"
  N = A + S
  nc = nt * pa + 2 * (1 - pa)
  nc2 = nt * sAGP + 2 * (1 - sAGP)
  
  if (nb == 1 & nt == 2) {   #normal diploid
    ccf = 2*(f/pu)
    ff = pu*cc/2
    Ms = computeSD(N, S, ff)
    ccf2 <- Ms$M1
    sd <- Ms$SD
  }
  else if (nt == 1) {        #het deletion
    ccf = (f/pu)*nc
    ff.C <- pu*cc/nc                       #dbinom
    Ms.C <- computeSD(N, S, ff.C)          #dbinom
    ccf2 <- Ms.C$M1                        #dbinom
    sd <- Ms.C$SD
    fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
    fl.ea <- (sAGP * (nt - nb))/nc2
    fh.t <- sAGP/nc2
    fh.e <- (1 - sAGP)/nc2
    pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
    pLate <- pbeta(fh.t, S+1, A+1)
    pEuploid <- pbeta(fh.e, S+1, A+1)
    Ptot <- pEarly.a + pLate + pEuploid
    cp.A <- pEarly.a/Ptot
    cp.CD <- 1 - cp.A
    cp.C <- pLate/Ptot
    cp.D <- pEuploid/Ptot
    cp.AC <- 1 - cp.D
    cp.AD <- 1 - cp.C
    if (cp.A >= sigTh) {
      evoType <- "A1"
    } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh) {
      evoType <- "B/C"
    } else if (cp.C >= sigTh) {
      evoType <- "B"
    } else if (cp.D >= sigTh) {
      evoType <- "C"
    } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh) {
      evoType <- "A1/B"
    } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh) {
      evoType <- "A1/C"
    }
  } else if (nb == 0 & nt == 0) {        #homozygous deletion
    evoType <- "C"
    ccf = (f*nc2)/sAGP
    ff.C = cc[cc<=(1-sAGP)]/nc2
    Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=(1-sAGP)])
    ccf2 <- (Ms.C$M1)/pu                    #dbinom
    sd <- Ms.C$SD
  } else if (nb == 0 | nt == 2 * nb) {   #NLOH or other balanced CNAs
    fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
    fl.ea <- (sAGP * (nt - nb))/nc2
    fh.t <- sAGP/nc2
    fh.e <- (1 - sAGP)/nc2
    pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
    pLate <- pbeta(fh.t, S+1, A+1)
    pEuploid <- pbeta(fh.e, S+1, A+1)
    Ptot <- pEarly.a + pLate + pEuploid
    cpEarly.a <- pEarly.a/Ptot
    cpLate.eup <- 1 - cpEarly.a
    cpLate <- pLate/Ptot
    cpEup <- pEuploid/Ptot
    if (Ptot > 0) {
      if (cpEarly.a >= sigTh){
        evoType <- "A1"
      } else if (cpLate.eup >= sigTh){
        evoType <- "B/C"
      } else if (cpLate >= sigTh){
        evoType <- "B"
      } else if (cpEup >= sigTh){
        evoType <- "C"
      }
    }
    allprobs = c(pEarly.a, pLate, pEuploid)
    names(allprobs) = c("pEarly.a", "pLate", "pEuploid")
    maxType = names(allprobs[match(max(allprobs),allprobs)])
    if (evoType == "A1" & prior != "late") {
      ccf = (f/pu)*nc - (nt - nb - 1)*pa
      ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
      Ms.A <- computeSD(N, S, ff.A)               #dbinom
      ccf2 <- Ms.A$M1                             #dbinom
      sd <- Ms.A$SD
    } else {
      ccf = (f/pu)*nc
      ff.C <- pu*cc/nc                        #dbinom
      Ms.C <- computeSD(N, S, ff.C)           #dbinom
      ccf2 <- Ms.C$M1                         #dbinom
      sd <- Ms.C$SD
    }
  } else if (nb >= 1 & nt > 2) {
    fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
    fl.ea <- (sAGP * (nt - nb))/nc2
    fh.eb <- (nb * sAGP + 1 - sAGP)/nc2
    fl.eb <- nb * sAGP/nc2
    fh.t <- sAGP/nc2
    fh.e <- (1 - sAGP)/nc2
    pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
    pEarly.b <- pbeta(fh.eb, S+1, A+1) - pbeta(fl.eb, S+1, A+1)
    pLate <- pbeta(fh.t, S+1, A+1)
    pEuploid <- pbeta(fh.e, S+1, A+1)
    Ptot <- pEarly.a + pEarly.b + pLate + pEuploid
    cp.A <- pEarly.a/Ptot
    cp.B <- pEarly.b/Ptot
    cp.C <- pLate/Ptot
    cp.D <- pEuploid/Ptot
    cp.AB <- 1 - cp.C - cp.D
    cp.AC <- 1 - cp.B - cp.D
    cp.AD <- 1 - cp.B - cp.D
    cp.BC <- 1 - cp.A - cp.D
    cp.BD <- 1 - cp.A - cp.C
    cp.CD <- 1 - cp.A - cp.B
    cp.ABC <- 1 - cp.D
    cp.ABD <- 1 - cp.C
    cp.ACD <- 1 - cp.B
    cp.BCD <- 1 - cp.A
    if (Ptot > 0) {
      if (cp.A >= sigTh) {                   # earl A
        evoType = "A1"
      } else if (cp.B >= sigTh){
        evoType <- "A2"
      } else if (cp.C >= sigTh){
        evoType <- "B"
      } else if (cp.D >= sigTh){
        evoType <- "C"
      } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh){
        evoType <- "B/C"
      } else if (cp.AB >= sigTh & cp.A < sigTh & cp.B < sigTh){
        evoType <- "A1/A2"
      } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh){
        evoType <- "A1/B"
      } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh){
        evoType <- "A1/C"
      } else if (cp.BC >= sigTh & cp.B < sigTh & cp.C < sigTh){
        evoType <- "A2/B"
      } else if (cp.BD >= sigTh & cp.B < sigTh & cp.D < sigTh){
        evoType <- "A2/C"
      } else if (cp.BCD >= sigTh & cp.BC < sigTh & cp.BD < sigTh & cp.CD < sigTh & cp.B < sigTh & cp.C < sigTh & cp.D < sigTh){
        evoType <- "A2/B/C"
      } else if (cp.ABC >= sigTh & cp.BC < sigTh & cp.AB < sigTh & cp.AC < sigTh & cp.B < sigTh & cp.C < sigTh & cp.A < sigTh){
        evoType <- "A1/A2/B"
      } else if (cp.ABD >= sigTh & cp.AB < sigTh & cp.AD < sigTh & cp.BD < sigTh & cp.B < sigTh & cp.D < sigTh & cp.A < sigTh){
        evoType <- "A1/A2/C"
      } else if (cp.ACD >= sigTh & cp.AC < sigTh & cp.AD < sigTh & cp.CD < sigTh & cp.A < sigTh & cp.D < sigTh & cp.C < sigTh){
        evoType <- "A1/B/C"
      }
    }
    allprobs = c(pEarly.a, pEarly.b, pLate, pEuploid)
    names(allprobs) = c("pEarly.a", "pEarly.b", "pLate", "pEuploid")
    maxType = names(allprobs[match(max(allprobs),allprobs)])
    if (evoType == "A1" & prior != "late") {
      ccf = (f/pu)*nc - (nt - nb - 1)*pa          #early A1
      ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
      #ff.A <- (cc - sAGP + (nt - nb) * sAGP)/nc2
      Ms.A <- computeSD(N, S, ff.A)               #dbinom
      ccf2 <- Ms.A$M1                             #dbinom
      #ccf2 = ccf2/pu
      sd <- Ms.A$SD
    } else if (evoType == "A2" & prior != "late") {    
      ccf = (f/pu)*nc - (nb - 1)* pa         #early A2
      ff.B <- pu*(cc - pa + nb * pa)/nc      #dbinom
      Ms.B <- computeSD(N, S, ff.B)          #dbinom
      ccf2 <- Ms.B$M1                        #dbinom
      sd <- Ms.B$SD
    } else {
      ccf = (f/pu)*nc                        #other
      ff.C <- pu*cc/nc                       #dbinom
      Ms.C <- computeSD(N, S, ff.C)          #dbinom
      ccf2 <- Ms.C$M1                        #dbinom
      sd <- Ms.C$SD
    }
  }
  if ( f > 0.1 & ccf >= overadj ) {          #over-adjustment
    if (evoType != "A1") {
      if ( (nt-nb) >= 3 ) {
        ccf = (f/pu)*2
      } else if ( nt >= 2 & (nt-nb) < 3 ) {
        ccf = (f/pu)*nc - (nt - nb - 1)*pa
      } else {
        ccf = (f/pu)*nc
      }
    } else {
      ccf = f*2
    }
  }
  if ( f > 0.7 & ccf2 < 0.05 & is.nan(sd)) {  #un-identified LOH
    ccf = f*2
    ccf2 = 1
    sd = 0.01
  }
  return(c(ccf, evoType, ccf2, sd))
}

computeSD <- function(N, S, f, cc=seq(0.02, 1, by = 0.01)) {
  M1list <- c()
  M2list <- c()
  MLElist <- c()
  for (ii in 1:length(N)) {
    PF <- sum(dbinom(S[ii], N[ii], f), na.rm = TRUE)
    M1 <- sum(dbinom(S[ii], N[ii], f) * cc, na.rm = TRUE)/PF
    M2 <- sum(dbinom(S[ii], N[ii], f) * cc^2, na.rm = TRUE)/PF
    M1list <- c(M1list, M1)
    M2list <- c(M2list, M2)
    MLElist <- c(MLElist, cc[which.max(dbinom(S[ii], N[ii], f))])
  }
  return(list(M1 = MLElist, SD = sqrt(M2list - M1list^2)))
}


partialRound <- function(x) {
  r = x
  for (i in 1:length(x)) {
    r[i] = round(x[i])
    if (abs(r[i]-x[i]) > 0.3 & x[i] > 1) {
      r[i] = round(x[i],1)
    }
    if (r[i] < 0) {
      r[i] = 0
    }
  }
  return(r)
}


effectivePurity <- function(p, Nt) {
  pu = Nt*p/(Nt*p + 2*(1-p))
  return(pu)
}


# classify mutations into private or public for a mutation table object
pubOrSub <- function(sampAB, samples, minAF=0.05, minDepTotal=5*length(samples), groupName = "") {
  
  originalColNames = colnames(sampAB)
  
  maxPa = vector()
  maxsAGP = vector()
  for (i in 1:length(samples)) {           #each sample
    sn = samples[i]
    pai = match(paste(sn, "pa", sep=""), originalColNames)
    sAGPi = match(paste(sn, "sAGP", sep=""), originalColNames)
    
    maxPa = c(maxPa, max(sampAB[,pai]))
    maxsAGP = c(maxsAGP, max(sampAB[,sAGPi]))
  }
  message(paste(maxPa,collapse="\t"))
  message(paste(maxsAGP,collapse="\t"))
  
  pstype = as.vector(apply(sampAB, 1, function(x, coln, samples, maxPa, maxsAGP, minAF, minDepTotal) {
    mutVector = as.vector(x)
    cppres = pubOrSub.Calc(x, coln, samples, maxPa, maxsAGP)
    #if (as.numeric(mutVector[2]) == 1390624 & mutVector[1] == "chr1") {
    #    message(paste(cppres, collapse="\t"))
    #}
    cpstype = "unknown"
    totalDepth = sum(as.numeric(mutVector[match(paste(samples,"d",sep=""),coln)]))
    if (totalDepth < minDepTotal) {
      cpstype = "unknown"
    } else if (cppres[1] >= 0.05 | cppres[3] == 1) {   #accept public
      cpstype = "public"
    } else if (cppres[2] >= 0.05) {                    #accept absent
      cpstype = "absent"
    } else {                                           #subclonal
      foundSamples = ""
      for (i in 1:length(samples)) {                 #get subclonal info for each sample
        sn = samples[i]
        mafai = match(paste(sn, "mafa", sep=""), coln)
        nbi = match(paste(sn, "nb", sep=""), coln)
        depthi = match(paste(sn, "d", sep=""), coln)
        if (as.numeric(mutVector[mafai]) > minAF) {
          foundSamples = paste(foundSamples, sn, ",", sep="")
        }
      }
      if (foundSamples != ""){
        foundSamples = gsub(",$","", foundSamples)
        cpstype = paste("private", foundSamples, sep="=")
      }
    }
    cpstype
  }, coln=originalColNames, samples=samples, maxPa = maxPa, maxsAGP = maxsAGP, minAF=minAF, minDepTotal=minDepTotal))
  
  if (groupName == "") {
    if ("pubOrSub" %in% originalColNames) {
      sampAB$pubOrSub = pstype
    } else {
      sampAB = data.frame(sampAB, pubOrSub=pstype)
      colnames(sampAB) = c(originalColNames,"pubOrSub")
    }
  } else {
    if (groupName %in% originalColNames) {
      sampAB[,groupName] = pstype
    } else {
      sampAB = data.frame(sampAB, pubOrSubTmp=pstype)
      colnames(sampAB) = c(originalColNames,groupName)
    }
  }
  return(sampAB)
  
}


# classify a particular mutation into private or public
pubOrSub.Calc <- function(mutVector, originalNames, samples, maxPa = vector(), maxsAGP = vector()) {
  
  CCFaboveOne = all((as.numeric(mutVector[match(paste(samples, "ccf", sep=""), originalNames)]) +
                       2.58*as.numeric(mutVector[match(paste(samples, "ccfSD", sep=""), originalNames)])) >= 1)
  aboveContri = length(which((as.numeric(mutVector[match(paste(samples, "ccf", sep=""), originalNames)]) +
                                2.58*as.numeric(mutVector[match(paste(samples, "ccfSD", sep=""), originalNames)])) >= 1))
  VAFaboveQua = all(as.numeric(mutVector[match(paste(samples, "mafa", sep=""), originalNames)]) >= 0.25)
  
  cpop = sapply(samples, function(x, samples, originalNames, mutVector, aboveContri, maxPa, maxsAGP) {
    sn = x
    si = match(x, samples)
    ccfi = match(paste(sn, "ccf", sep=""), originalNames)
    ccfsdi = match(paste(sn, "ccfSD", sep=""), originalNames)
    mafai = match(paste(sn, "mafa", sep=""), originalNames)
    refci = match(paste(sn, "refc", sep=""), originalNames)
    altci = match(paste(sn, "altc", sep=""), originalNames)
    pui = match(paste(sn, "pu", sep=""), originalNames)
    pai = match(paste(sn, "pa", sep=""), originalNames)
    sAGPi = match(paste(sn, "sAGP", sep=""), originalNames)
    nti = match(paste(sn, "nt", sep=""), originalNames)
    nbi = match(paste(sn, "nb", sep=""), originalNames)
    depthi = match(paste(sn, "d", sep=""), originalNames)
    
    cPa = as.numeric(mutVector[pai])
    prior = "unknown"
    if (cPa == maxPa[si] | cPa == 0) {
      cPa = 1
    }
    if (aboveContri >= 2 | (aboveContri >= 1 & length(samples) == 2)) {                    #looks like public, change prior
      prior = "pub"
    }
    cpop.probs = pubOrSub.prob.binom(as.integer(mutVector[refci]), as.integer(mutVector[altci]), as.numeric(maxsAGP[si]), as.numeric(mutVector[sAGPi]), as.integer(mutVector[nti]), as.integer(mutVector[nbi]), pa=cPa, prior=prior)
    cpop.probs
  }, originalNames=originalNames, samples=samples, mutVector=mutVector, aboveContri=aboveContri, maxPa=maxPa, maxsAGP=maxsAGP)
  
  cpp = prod(cpop[1,])
  cpa = prod(cpop[2,])
  
  return(c(cpp, cpa, as.numeric(CCFaboveOne | VAFaboveQua)))
  
}

# classify mutations into private or public for a mutation table object, from simulated virtual tumor
pubOrSub.simu <- function(sampAB, samples, minAF=0.05, minDepTotal=5*length(samples), groupName = "", pAF=0.25) {    #for simulated data
  
  originalColNames = colnames(sampAB)
  
  pstype = as.vector(apply(sampAB, 1, function(x, coln, samples, pAF, minAF, minDepTotal) {
    mutVector = x
    names(mutVector) = coln
    cppres = pubOrSub.Calc.sim(mutVector, samples, rep(1,length(samples)), rep(1,length(samples)), pAF=pAF)
    cpstype = "unknown"
    totalDepth = sum(as.numeric(mutVector[paste("depth",samples,sep="")]))
    if (totalDepth < minDepTotal) {
      cpstype = "unknown"
    } else if (cppres[1] >= 0.05 | cppres[3] == 1) {   #accept public
      cpstype = "public"
    } else if (cppres[2] >= 0.05) {  #accept absent
      cpstype = "absent"
    } else {     #subclonal
      foundSamples = ""
      for (i in 1:length(samples)) {                                     #get subclonal ones 2nd round: Refine
        sn = samples[i]
        mafai = match(paste("maf", sn, sep=""), coln)
        depthi = match(paste("depth", sn, sep=""), coln)
        if (as.numeric(mutVector[mafai]) > minAF) {
          foundSamples = paste(foundSamples, sn, ",", sep="")
        }
      }
      if (foundSamples != "") {
        foundSamples = gsub(",$","", foundSamples)
        cpstype = paste("private", foundSamples, sep="=")
      }
    }
    cpstype
  }, coln=originalColNames, samples=samples, pAF=pAF, minAF=minAF, minDepTotal=minDepTotal))
  
  if (groupName == "") {
    if ("pubOrSub" %in% originalColNames) {
      sampAB$pubOrSub = pstype
    } else {
      sampAB = data.frame(sampAB, pubOrSub=pstype)
      colnames(sampAB) = c(originalColNames,"pubOrSub")
    }
  } else {
    if (groupName %in% originalColNames) {
      sampAB[,groupName] = pstype
    } else {
      sampAB = data.frame(sampAB, pubOrSubTmp=pstype)
      colnames(sampAB) = c(originalColNames,groupName)
    }
  }
  return(sampAB)
}

# classify a particular mutation into private or public, for simulated virtual tumor in pure "diploid" stats
pubOrSub.Calc.sim <- function(mutVector, samples, maxPa = vector(), maxsAGP = vector(), pAF=0.25) {
  
  originalNames = names(mutVector)
  
  aboveContri = 0
  VAFaboveQua = all(as.numeric(mutVector[match(paste("maf", samples, sep=""), originalNames)]) >= pAF)
  
  cpop = sapply(samples, function(x, originalNames, mutVector) {
    sn = x
    mafai = match(paste("maf", sn, sep=""), originalNames)
    depthi = match(paste("depth", sn, sep=""), originalNames)
    depth = as.numeric(mutVector[depthi])
    maf = as.numeric(mutVector[mafai])
    refc = round(depth * (1-maf))
    altc = round(depth * maf)
    
    cPa = 1
    prior = "unknown"
    
    cpop.probs = pubOrSub.prob.binom(refc, altc, 1, 1, 2, 1, pa=1, prior=prior)
    cpop.probs
  }, originalNames=originalNames, mutVector=mutVector)
  
  cpp = prod(cpop[1,])
  cpa = prod(cpop[2,])
  
  return(c(cpp, cpa, as.numeric(VAFaboveQua)))
  
}


# probability of being public or being absent
pubOrSub.prob.binom <- function(A, S, pu, sAGP, nt, nb, pa=1, prior="unknown") {      #using binomial test
  
  N = A + S
  
  nc = nt * sAGP + 2 * (1 - sAGP)
  f.abs <- 0.02
  f.pub <- (pu * (nt - nb))/nc
  
  f.pub = ifelse(nb == 0 & nt == 0, 0.5*pu,
                 ifelse(nb == 0 & prior == "pub", max((pu - sAGP)/nc, 0),
                        ifelse(nt >= 2, pu/nc, (pu * (nt - nb))/nc)))
  
  f.pub = min(1,f.pub)
  #p.pub = pbeta(fh.pub,S+1,A+1) - pbeta(fl.pub,S+1,A+1)
  p.pub = pbinom(S, N, f.pub)       #p for pub (under pub vaf, the prob reaching S)
  #p.abs = 1 - pbinom(S, N, f.abs)   #p for abs (under abs, the prob reaching S)
  p.abs = pbeta(f.abs, S+1, A+1)
  return(c(p.pub, p.abs))
  
}


#Fst hudson estimate
#af = data.frame( maf1, depth1, maf2)
fst.hudson <- function(af, minAF=0.08) {  #combine for multi SSNV sites using ratio of averages
  mafis = which(grepl("maf", colnames(af)))
  keep = as.vector(apply(af, 1, function(x, mafis) {
    maxmaf = max(as.numeric(x[mafis]))
    if (maxmaf > minAF) {
      TRUE
    } else {
      FALSE
    }
  }, mafis=mafis))     #filter data
  af = af[keep,]
  
  Ns = c()
  Ds = c()
  for(k in 1:nrow(af)) {
    n1 = af$depth1[k]
    n2 = af$depth2[k]
    p1 = af$maf1[k]
    p2 = af$maf2[k]
    N = (p1-p2)^2-(p1*(1-p1))/(n1-1)-(p2*(1-p2))/(n2-1)     # covariance
    D = p1*(1-p2)+p2*(1-p1)                                 # standard deviations
    Ns = c(Ns, N)
    Ds = c(Ds, D)
  }
  Fst.h = mean(Ns)/mean(Ds)
  return(Fst.h)
}



#allSubRows = union(subAi,subBi)
#mafs = paste(c(snA,snB), "mafa", sep="")
#depths = paste(c(snA,snB), "d", sep="")
#rAUC(sampAB[allSubRows,], 0.04, mafs, depths)

#'rAUC: calculation of rAUC
#' contrast the theoretical sf with neutral ones. VAFs within 0.04 - 0.25.

#calculation of rAUC
rAUC <- function(data, lower, mafs, depths) {
  
  if (length(mafs) <= 8) {
    lower = round((0.08/length(mafs)),2)
  } else {
    lower = max(0.005, round((0.08/length(mafs)),3))
  }
  #message(paste("lower: ", lower,sep=""))
  weightAF = weightAFs(data, mafs, depths)
  #message(paste(weightAF, collapse=" "))
  weightAF = weightAF[which(weightAF >= lower & weightAF <= 0.25)]
  #message(paste("weightedMuts: ", length(weightAF), sep=""))
  vafs = seq(lower,0.25,by=0.01)
  nstep = length(vafs)
  counts = vector()
  ncounts = ((1/vafs)-(1/0.25))/((1/lower)-(1/0.25))
  
  for ( i in 1:length(vafs) ) {
    counts = append(counts, length(which(weightAF > vafs[i]))-length(which(weightAF > vafs[nstep])))
  }
  counts = counts/counts[1]
  
  #bei = bezierCurve(vafs, counts, length(vafs))
  #AUC = trapz(bei$x,bei$y)
  
  #modify
  AUC = trapz(vafs,counts)
  
  #AUC = trapz(mafs,counts)
  nAUC = trapz(vafs,ncounts)
  rAUC = AUC/nAUC
  rAUC = round(rAUC,8)
  return(list(rAUC=rAUC,weightAF=weightAF))
  
}

# weight AF for multi samples
weightAFs <- function(af.data, mafs, depths, minAF=0) {
  
  keep = vector()
  for (i in 1:length(mafs)) {
    cmaf = mafs[i]
    cmafi = match( cmaf, colnames(af.data) )
    if (length(keep) == 0) {
      keep = af.data[,cmafi] > minAF
    } else {
      keep = keep | (af.data[,cmafi] > minAF)
    }
  }
  af.data = af.data[keep,]
  #message(paste("keepRows4AUC: ", dim(af.data)[1], sep=""))
  
  tdepth = rep(0,dim(af.data)[1])
  talt = rep(0,dim(af.data)[1])
  for (i in 1:length(mafs)) {
    cmaf = mafs[i]
    cmafi = match(cmaf,colnames(af.data))
    cdep = depths[i]
    cdepi = match(cdep,colnames(af.data))
    tdepth = tdepth + af.data[,cdepi]
    #message(paste(tdepth, collapse=" "))
    talt = talt + af.data[,cdepi]*af.data[,cmafi]
  }
  tmaf = talt/tdepth
  return(tmaf)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plotting simVAF data and generate stats for virtual tumors


#for plotting simulated vaf data
plotRes.simVAF.matrix.pdf <- function(sampAB, samples, depths, pdfsize = 16, plotType = "AF", snr="", sns=vector()) {
  combinations = combn(length(samples),2)
  nplots = dim(combinations)[2]
  ndims = length(samples)-1
  
  resStats = list()
  outfile = paste(snr,"multi_hist.pdf",sep="")
  if (plotType == "Scatter") {
    outfile = paste(snr,"multi_scatter.pdf",sep="")
  }
  if (plotType == "Density") {
    outfile = paste(snr,"multi_density.pdf",sep="")
  }
  pdf(file = outfile, width=pdfsize, height=pdfsize)
  layout(matrix(seq(1,ndims^2), ndims, ndims, byrow = FALSE))
  #par(mfcol=c(ndims,ndims))
  pindex = 0
  for (ci in 1:ndims) {                           #for each column
    noplotRow = vector()
    if (ci > 1) {
      noplotRow = (1:ndims)[1:(ci-1)]         #which row do not plot?
    }
    for (ri in 1:ndims) {                       #for each row
      if (ri %in% noplotRow) {                #no plot
        plot(NULL,NULL,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
      } else {                                #plot
        pindex = pindex + 1
        pair = combinations[,pindex]
        sn1 = samples[pair[1]]              #samplename1
        sn1mafa = sn1
        sn1d = depths[pair[1]]
        sn2 = samples[pair[2]]              #samplename2
        sn2mafa = sn2
        sn2d = depths[pair[2]]
        main.title = paste(snr, sns[pair[1]], "vs", sns[pair[2]], sep=" ")
        statName = paste(snr, sns[pair[1]], sns[pair[2]], "stats", sep="_")
        if (plotType == "AF") {
          resStats[[statName]] = plotRes.simVAF.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], dp1=depths[pair[1]], dp2=depths[pair[2]],
                                                    sn1=sn1mafa, sn2=sn2mafa, plotDensity=F, plotScatter=F, pdf=F)
        } else if (plotType == "Scatter") {
          resStats[[statName]] = plotRes.simVAF.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], dp1=depths[pair[1]], dp2=depths[pair[2]],
                                                    sn1=sn1mafa, sn2=sn2mafa, plotDensity=F, plotAF=F, pdf=F)
        } else if (plotType == "Density") {
          resStats[[statName]] = plotRes.simVAF.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], dp1=depths[pair[1]], dp2=depths[pair[2]],
                                                    sn1=sn1mafa, sn2=sn2mafa, plotScatter=F, plotAF=F, pdf=F)
        }
      }
    }
  }
  dev.off()
  return(1)
}



plotRes.simVAF.pdf <- function(sampAB, sampName, main=sampName, sn1n="A", sn2n="B", sn1="maf1", sn2="maf2", dp1="depth1", dp2="depth2",minAF=0.05, ratio=1, plotAF=TRUE, pdf=TRUE, alpha=1, binw=0, reportSSR=FALSE) {
  
  sn1s = sn1n     #name of sample 1
  sn2s = sn2n     #name of sample 2
  
  dp1i = match(dp1,colnames(sampAB))
  dp2i = match(dp2,colnames(sampAB))
  
  maf1Index = match(sn1, colnames(sampAB))
  maf2Index = match(sn2, colnames(sampAB))
  
  #check depth power to reject a presence of a mutation
  depthPowerKeep <- as.vector(apply(sampAB, 1, function(x,maf1i,maf2i,dp1i,dp2i) {
    if(as.numeric(x[maf1i]) == 0){vaf = as.numeric(x[maf2i])
    if (pbinom(0,as.numeric(x[dp1i]),vaf) < 0.05){TRUE} else {FALSE}}
    else if (as.numeric(x[maf2i]) == 0){vaf = as.numeric(x[maf1i])
    if (pbinom(0,as.numeric(x[dp2i]),vaf) < 0.05){TRUE} else {FALSE}}
    else {TRUE}
  }, maf1i=maf1Index,maf2i=maf2Index,dp1i=dp1i,dp2i=dp2i))
  sampAB = sampAB[depthPowerKeep,]
  
  subMuts = subclonalMutSim(sampAB, sn1, sn2, dp1, dp2, minAF=minAF)
  
  allA_Rows = which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1)                                                       #all sample A rows
  subA_Rows = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))                             #sample A sub rows
  ssA_Rows  = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1 & sampAB[,maf2Index] == 0))   #sample A specific rows
  
  allB_Rows = which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1)                                                       #all sample B rows
  subB_Rows = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))                             #sample B sub rows
  ssB_Rows  = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1 & sampAB[,maf1Index] == 0))   #sample B specific rows
  allAB_Rows = union(allA_Rows, allB_Rows)
  
  ssR = length(union(ssA_Rows,ssB_Rows))/length(union(subA_Rows,subB_Rows))
  if (reportSSR == TRUE) {
    return(ssR)
  }
  
  BinWidthA = round(dpih(sampAB[allA_Rows, maf1Index]/ratio),2)
  if (BinWidthA == 0) { BinWidthA = 0.02 }
  BinWidthB = round(dpih(sampAB[allB_Rows, maf2Index]/ratio),2)
  if (BinWidthB == 0) { BinWidthB = 0.02 }
  BinWidth = min(c(BinWidthA, BinWidthB, 0.1))
  if ( binw != 0 ) {
    BinWidth = binw
  }
  #message(paste("bin width: ", BinWidthA, BinWidthB, BinWidth, sep=" "))
  nbreaksA = round((max(sampAB[allA_Rows, maf1Index]/ratio)-min(sampAB[allA_Rows, maf1Index]/ratio)+0.01)/BinWidth)
  nbreaksB = round((max(sampAB[allB_Rows, maf2Index]/ratio)-min(sampAB[allA_Rows, maf2Index]/ratio)+0.01)/BinWidth)
  nbreaks = ceiling(diff(range(minAF,1))/BinWidth)
  #message(paste("nbreaks: ", nbreaksA, nbreaksB, nbreaks, sep=" "))
  breaksA = seq(minAF,1,length.out=nbreaks)
  breaksB = seq(minAF,1,length.out=nbreaks)
  
  
  sampAh = hist(sampAB[allA_Rows, maf1Index]/ratio, breaks=breaksA,  plot=F)
  sampAhsub = hist(sampAB[subA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  sampAhss = hist(sampAB[ssA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  ylimup = max(sampAh$count)
  
  
  sampBh = hist(sampAB[allB_Rows, maf2Index]/ratio, breaks=breaksB, plot=F)
  sampBh$counts = sampBh$counts*(-1)
  sampBhsub = hist(sampAB[subB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhsub$counts = sampBhsub$count*(-1)
  sampBhss = hist(sampAB[ssB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhss$counts = sampBhss$counts*(-1)
  ylimdown = min(sampBh$count)
  
  if (abs(ylimup) >= abs(ylimdown)) {
    ylimdown = (-1)*ylimup
  } else {
    ylimup = (-1)*ylimdown
  }
  
  if (plotAF == TRUE) {
    
    if (pdf == TRUE) {
      pdf(file = paste(sampName, "hist.pdf", sep="_"), width = 8, height = 8, useDingbats=FALSE)
    }
    par(mar=c(4.5,5,4.5,0))
    
    plot( sampAh, col=rgb(0,0,0,1/4), xlim=c(0, 1), ylim=c(ylimdown,ylimup), border=F, ylab="# of Mutations", xlab="Allele Frequency", axes = F, main = main,cex.lab = 2.9, cex.main = 2.9)    # first histogram
    plot( sampAhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green no border
    plot( sampAhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
    
    
    plot( sampBh,col=rgb(0,0,0,1/4), border=F, add=T )  # second histogram
    plot( sampBhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green
    plot( sampBhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
    sampAhss$counts = 0                                  #make black line
    sampBhss$counts = 0                                  #make black line
    plot( sampAhss, col="black", add=T, border=F)        #make black line
    plot( sampBhss, col="black", add=T, border=F)        #make black line
    
    axis(side=1,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=2.4)       #x-axis
    axis(side=2,at=decideTickAt(ylimdown, ylimup),labels=allAbs(decideTickAt(ylimdown, ylimup)),cex.axis=2.4)
    text(x=0.5,y=ylimdown,labels=paste(sn2s,",",length(which(sampAB[,maf2Index] > minAF)), "SSNVs", sep = " "), cex=2.8)
    text(x=0.5,y=ylimup,labels=paste(sn1s,",",length(which(sampAB[,maf1Index] > minAF)), "SSNVs", sep = " "), cex=2.8)
    
    #stats starting from here!
    fHsub = round(mean(c(subMuts$ratioHighSubA,subMuts$ratioHighSubB)),3)
    fst = round(subMuts$FST,3)
    ksd = round(subMuts$KSD,3)
    text(x=0.8,y=3*ylimup/4, labels=bquote(paste("fH"["sub"], " = ", .(fHsub))),cex=2.6)
    text(x=0.8,y=(3/4-0.167)*ylimup, labels=bquote(paste("FST", " = ", .(fst))),cex=2.6)
    text(x=0.8,y=(3/4-0.334)*ylimup, labels=bquote(paste("KSD", " = ", .(ksd))),cex=2.6)
    
    npub = subMuts$pubTn
    legend(0.5,ylimdown/4, legend=c(paste("Public ","(",npub,")",sep=""),"Pvt-Shared","Rgn Specific"),
           col=c(rgb(0,0,0,1/4),rgb(178/255,223/255,138/255,1),rgb(31/255,120/255,180/255,1)), pch=15, bty="n", cex=2.4)
    if (pdf == TRUE) {
      dev.off()
    }
  }
  
  return(subMuts)
}



subclonalMutSim <- function(sampAB, snA, snB, dpA, dpB, minAF=0.05, statsAF=0.08, highAF=0.2, ratio=1, pob="pubOrSub")   {                  #determinine subclonal mutations
  mafaAi = match(snA, colnames(sampAB))
  mafaBi = match(snB, colnames(sampAB))
  depthAi = match(dpA, colnames(sampAB))
  depthBi = match(dpB, colnames(sampAB))
  pobi = match(pob, colnames(sampAB))
  
  # for KSD
  subAi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaAi] > minAF )    
  ssAi  = intersect(subAi, which( sampAB[,mafaAi] > minAF & sampAB[,mafaBi] == 0 ))
  mutsA = sampAB[subAi,mafaAi]/ratio
  
  subBi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaBi] > minAF )
  ssBi  = intersect(subBi, which( sampAB[,mafaBi] > minAF & sampAB[,mafaAi] == 0 ))
  mutsB = sampAB[subBi,mafaBi]/ratio
  
  KSD = as.numeric(ks.test( mutsA[which(mutsA > statsAF)], mutsB[which(mutsB > statsAF)] )$statistic)
  
  # pub
  pubTi = which( sampAB[,pobi] == "public" )
  
  # FST
  allSubRows = union(subAi,subBi)
  mutsSub = sampAB[allSubRows,]
  mutsSub = data.frame( maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi] )
  FST = fst.hudson(mutsSub, minAF=statsAF)
  
  # for other stats
  mutsA2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF )), mafaAi]
  mutsAh2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF )), mafaAi]
  mutsASr2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] > 0.02)), mafaAi]   #shared
  
  mutsASp2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] == 0)), mafaAi]
  mutsASph2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF & sampAB[,mafaBi] == 0)), mafaAi]
  
  mutsB2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF )), mafaBi]
  mutsBh2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF )), mafaBi]
  mutsBSr2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] > 0.02)), mafaBi]  #shared
  
  mutsBSp2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] == 0)), mafaBi]
  mutsBSph2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF & sampAB[,mafaAi] == 0)), mafaBi]
  
  # mutational counts
  pubTn = length(pubTi)
  ssAn = length(ssAi)
  ssBn = length(ssBi)
  sharedAn = length(mutsASr2)
  sharedBn = length(mutsBSr2)
  
  # list for output    
  muts = list(A=mutsA,B=mutsB,subAi=subAi,subBi=subBi, ssAi=ssAi, ssBi=ssBi, pubTi=pubTi, #subArow=subArow, subBrow=subBrow, pubTrow = pubTrow,
              lenSubA=length(mutsA2),lenSubAh=length(mutsAh2),ratioHighSubA=length(mutsAh2)/length(mutsA2),
              lenSubB=length(mutsB2),lenSubBh=length(mutsBh2),ratioHighSubB=length(mutsBh2)/length(mutsB2),
              lenSsA=length(mutsASp2),lenHighSsA=length(mutsASph2),ratioHighSsA=length(mutsASph2)/length(mutsASp2),pSsA=length(mutsASp2)/length(mutsA2),
              lenSsB=length(mutsBSp2),lenHighSsB=length(mutsBSph2),ratioHighSsB=length(mutsBSph2)/length(mutsBSp2),pSsB=length(mutsBSp2)/length(mutsB2),
              lenSharedA=length(mutsASr2), lenSharedB=length(mutsBSr2), ratioSharedA=length(mutsASr2)/length(mutsA2), ratioSharedB=length(mutsBSr2)/length(mutsB2),
              FST=FST, KSD=KSD, pubTn=pubTn, sharedAn=sharedAn, sharedBn=sharedBn, ssAn=ssAn, ssBn=ssBn)
  
  return(muts)
}



###### SVM training based on virtual tumors #################################
trainSVM <- function(data, lent, featureCols=2:5, classCol=1, ncores=2, trainY="", subSample=FALSE, seed=1943) {
  registerDoMC(cores = ncores)
  
  x = apply(data[,featureCols], 2, as.numeric)
  x = apply(x, 2, function(x) {
    mean.pool = mean(x)
    sd.pool = sd(x)
    (x-mean.pool)/sd.pool
  })
  
  trainX = x[(lent+1):dim(x)[1],]
  testX = x[1:lent,]
  if (trainY == "") {
    y = data[,classCol]
    trainY = y[(lent+1):length(y)]
    trainY = sapply(trainY, function(x){if (x == "s=5" | x == "s=10" | x == "s=2" | x == "s=3")      # group into eneutral and selection modes
    {"selection"} else {"eneutral"}})
    trainY = as.factor(trainY)
  } else {
    trainY = trainY
  }
  testX = trainX
  testY = trainY
  if (subSample == TRUE) {
    set.seed(seed)
    tSize = length(trainY)
    sSize = round(tSize/5)
    message(paste("testSize and trainSize:", tSize,sSize,sep=" "))
    testI = sample(tSize, sSize)
    keepI = setdiff(1:tSize, testI)
    testX = trainX[testI,]
    testY = trainY[testI]
    trainX = trainX[keepI,]
    trainY = trainY[keepI]
    message(paste("trainSize:", length(trainY), sep=" "))
  }
  
  ## SVM start
  # First pass
  set.seed(seed)
  # Setup for cross validation
  ctrl <- trainControl(method="repeatedcv",                   # 10fold cross validation
                       repeats=5,                             # do 5 repititions of cv
                       summaryFunction=twoClassSummary,       # Use AUC to pick the best model
                       classProbs=TRUE)
  
  #Train and Tune the SVM
  message("first round training")
  
  svm.tune <- train(x=trainX,
                    y=trainY,
                    method = "svmRadial",                     # Radial kernel
                    tuneLength = 15,                          # 15 values of the cost function
                    #tuneGrid = grid,
                    metric="ROC",
                    trControl=ctrl,
                    scaled = FALSE)
  
  sigma1 = as.numeric(svm.tune$bestTune["sigma"])
  message(sigma1)
  s_incre = round(sigma1/20, 2)
  message(s_incre)
  sigmaTestRange = seq(sigma1-2*s_incre, sigma1+2*s_incre, by=s_incre)
  C1 = as.numeric(svm.tune$bestTune["C"])
  message(C1)
  c1_incre = round(C1/20, 2)
  message(c1_incre)
  CTestRange = seq(C1-2*c1_incre, C1+2*c1_incre, by=c1_incre)
  
  # Second pass
  set.seed(seed)
  # Use the expand.grid to specify the search space   
  grid <- expand.grid(sigma = sigmaTestRange, C = CTestRange)
  
  #Train and Tune the SVM
  message("second round training: refining sigma and C")
  svm.tune <- train(x=trainX,
                    y= trainY,
                    method = "svmRadial",
                    #preProc = c("center","scale"),
                    metric="ROC",
                    tuneGrid = grid,
                    trControl = ctrl,
                    scaled = FALSE)
  roc = ""
  if (subSample == TRUE) {
    pred = predict.train(svm.tune, testX)
    roc = roc(testY, as.numeric(pred))
    svm.tune = list(svm.tune=svm.tune, roc = roc)
  }
  return(svm.tune)
}

