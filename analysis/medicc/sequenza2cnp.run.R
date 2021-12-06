#object:
#to prepare the input formate of medicc, we need to transform the seuquenza into copy number profiles(cnps).

library(GenomicRanges)
library(magrittr)


#write into fasta
write.fasta = function(merge_A, major = "major", out.dir = "data", tumor = "tumor"){
  system(sprintf("mkdir %s", out.dir))
  chrs = unique(merge_A$seqnames)[!( unique(merge_A$seqnames) %in% c("chrX", "chrY")) ]
  num = ncol(merge_A)
  
  for(i in chrs){
    merge_B = subset(merge_A, seqnames == i)
    text = rep("0", 2*(num - 3))
    text[1] = ">diploid"
    text[2] = paste0( rep(1, nrow(merge_B)), collapse = "")
    
    for(j in 5:num){
      text[ 2*(j-3) -1 ] = sprintf(">%s", colnames(merge_B)[j] )
      text[ 2*(j-3)] = merge_B[,j] %>% stringr::str_c(collapse = "")
    }
    write.table(text, file = sprintf("%s/%s_%s_%s.fasta", out.dir, tumor ,major, i), quote = F, col.names = F, row.names = F )
  }
  
  #write desc
  chrinfo = rep("0", length(chrs))
  for(i in 1:length(chrs)){
    chrinfo[i] = sprintf("%s %s %s", chrs[i], 
                         sprintf("%s_%s_%s.fasta", tumor ,"major", chrs[i]), 
                         sprintf("%s_%s_%s.fasta", tumor ,"minor", chrs[i]) )
  }
  write.table(chrinfo, file = sprintf("%s/%s.descr.txt", out.dir, tumor), quote = F, col.names = F, row.names = F )
  
}

#combind data

#' @param segfiles: segment files position.
#' @param sampleid: the corresponding sample ids.
#' @param out.dir: output dir
#' @param tumor: the output prefix
#' @param N.baf: quality control for the sequenza output
#' @param cnv_min_length: quality control for the sequenza output
#' @param max_CNt: quality control for the sequenza output
#' @param minLength: the min length of CNVs to output.
#' @param maxCNV: the max CNV to output. When the raw CNV greater than maxCNV, then its value was set to maxCNV.


red.seg =  function(segfiles, sampleid, out.dir = "data",  tumor = "tumor", N.baf = 30, cnv_min_length = 1e5, max_CNt = 10 , minLength = 1e5, maxCNV = 4 ){
  
  #Parameters for filtering before downstream analysis.
  #N.baf
  #cnv_min_length
  #max_CNt
  
  #Parameters for finial output 
  #minLength = 1e5, , maxCNV = 4
  
  
  seglist = list()
  #read segs
  for(i in 1:length(sampleid)){
    
    seg = read.delim(file = segfiles[i], header = T, stringsAsFactors = F) %>%
      #dplyr::filter(chromosome == "chr1") %>%
      dplyr::mutate(sample =sampleid[i] )
    
    seglist[[i]] = seg
  }
  
  #removing low-confidence regions
  #seglist = base::Reduce(rbind, seglist) %>%
  #  dplyr::filter(N.BAF >= N.baf &  (end.pos - start.pos) >=  cnv_min_length & CNt <= max_CNt)
  
  #change low-confidence regions into 2,1,1
  seglist = base::Reduce(rbind, seglist) %>%
    dplyr::mutate( keep = ifelse( N.BAF >= N.baf &  
                                    (end.pos - start.pos) >=  cnv_min_length & 
                                    CNt <= max_CNt, TRUE, FALSE)) %>%
    dplyr::mutate(CNt = ifelse(keep, CNt, 2),
                  A = ifelse(keep, A, 1),
                  B = ifelse(keep, B, 1)
    )
  
  
  
  gseg = GRanges(seqnames = seglist$chromosome,
                 ranges = IRanges(seglist$start.pos, seglist$end.pos),
                 strand = "+"
                 
  )
  
  #metadata columns can be added to a GRanges object
  mcols(gseg) = seglist
  
  
  #split regions into small regions.
  segdis =  GenomicRanges::disjoin(gseg)
  #add region infor to data
  mcols(segdis) = data.frame(segdis)
  
  #set the minimal site of each segs.
  #minLength = 1e5
  
  #get overlaps regions
  overlaps =  findOverlaps(segdis, gseg )
  
  #combined information
  merge = cbind( mcols(segdis[queryHits(overlaps) , ]) , mcols(gseg[subjectHits(overlaps) ,] ) )  %>%
    as.data.frame() %>%
    #filter seg length
    dplyr::filter( width >= minLength) %>%
    dplyr::select(seqnames, start, end, width,  CNt, A, B, sample)
  
  #removing neutral region. CNt = 2 and A = 1 and B =1.
  
  merge_summary = merge %>%
    group_by(seqnames, start, end, width, CNt, A, B) %>%
    summarise( num = n()) %>%
    filter(num == length(sampleid) )
  
  merge = left_join(merge, merge_summary) %>%
    filter(is.na(num) ) %>%
    mutate(num = NULL)
  
  #major info: A
  merge_A = merge %>%
    dplyr::select(seqnames, start, end, width, A, sample) %>%
    tidyr::spread(key = sample, value = A, fill = 1 )
  #change cnv number when it >=10
  tmp = merge_A[, 5:ncol(merge_A) ]
  tmp[tmp>= maxCNV ] = maxCNV
  merge_A = cbind(merge_A[,1:4], tmp)
  
  #minor info: B
  merge_B = merge %>%
    dplyr::select(seqnames, start, end, width, B, sample) %>%
    tidyr::spread(key = sample, value = B, fill = 1 )
  
  tmp = merge_B[, 5:ncol(merge_B) ]
  tmp[tmp>= maxCNV ] = maxCNV
  merge_B = cbind(merge_B[,1:4], tmp)
  
  if(!file.exists(out.dir)){
    system( paste0("mkdir ", out.dir))
  }
  
  file = sprintf("%s/%s.descr.txt", out.dir, tumor)
  write.table(merge_A, file = sprintf("%s/%s.major.txt", out.dir, tumor), quote = F, row.names = F, sep = "\t" )
  write.table(merge_B, file = sprintf("%s/%s.minor.txt", out.dir, tumor), quote = F, row.names = F, sep = "\t" )
  
  write.fasta(merge_A = merge_A, major = "major", out.dir = out.dir, tumor =  tumor )
  write.fasta(merge_A = merge_B, major = "minor", out.dir = out.dir, tumor =  tumor )
  
  return(list(major = merge_A, minor=merge_B))
}



#write.seg
write.seg = function(segfiles, sampleid, out.dir = "data", N.baf = 30, cnv_min_length = 1e5, max_CNt = 10){
  #for Lung and UterusM, we change the cnvs.
  
  if(!file.exists(out.dir)){
    system( paste0("mkdir ", out.dir))
  }
  
  seglist = list()
  #read segs
  for(i in 1:length(sampleid)){
    
    seg = read.delim(file = segfiles[i], header = T, stringsAsFactors = F) %>%
      dplyr::mutate(sample =sampleid[i] ) %>%
      dplyr::mutate( keep = ifelse( N.BAF >= N.baf &  
                                      (end.pos - start.pos) >=  cnv_min_length & 
                                      CNt <= max_CNt, TRUE, FALSE)) %>%
      dplyr::mutate(CNt = ifelse(keep, CNt, 2),
                    A = ifelse(keep, A, 1),
                    B = ifelse(keep, B, 1)
      )
    file = sprintf("%s/%s_segments.txt", out.dir, sampleid[i])
    message( sprintf("Writing: %s", file ))
    write.table(seg %>% mutate(sample = NULL, keep = NULL), 
                file = file, quote = F, row.names = F, sep = "\t" )
    
    seglist[[i]] = seg
  }
  
  seg
}



# running -----------------------------------------------------------------

#input
#folder = "/data1/qingjian/Rproject/Three/medicc/Seg.new/"

#segfiles = list.files( folder, pattern = ".txt", full.names = T)
#sampleid = list.files( folder, pattern = ".txt") %>% stringr::str_remove("_segments.txt")

#for each cancer types.
#tumors = c("Breast","Coad","Lung","OveryLM","OveryRM","UterusM")

#for(i in tumors){
#  data = red.seg(segfiles = segfiles[grepl(segfiles, pattern = i)], 
#          sampleid = sampleid[grepl(sampleid, pattern = i)], out.dir = i , tumor = i, minLength #= 1e6)
#}



#for panCancer, combined all cancer into one types.
#data = red.seg(segfiles = segfiles, 
#        sampleid = sampleid , out.dir = "fasta/pan1", tumor = "pan1", N.baf = 30,
#        cnv_min_length = 5e6, max_CNt = 10 , 
#        minLength = 5e6, maxCNV = 4
#        )



#major = data$major
#minor = data$minor





################################################################################

#New run on 2020-10-27

#for Lung and UterusM, their cnv changes are rare, we set a more strict criteria to filtering the cnv, removing the potential false positive cnvs.

#change Lung and UterusM

#i = "Lung"
#data = write.seg(segfiles = segfiles[grepl(segfiles, pattern = i)],
#                  sampleid = sampleid[grepl(sampleid, pattern = i)], 
#                  out.dir = "Seg.new",
#                  N.baf = 30, max_CNt = 4 , cnv_min_length = 1e6)

#i = "UterusM"
#data = write.seg(segfiles = segfiles[grepl(segfiles, pattern = i)],
#                 sampleid = paste("UterusM", c(2,4:7), sep = "_"), 
#                 out.dir = "Seg.new",
#                 N.baf = 30, max_CNt = 4 , cnv_min_length = 1e6)

#keep U1 and U3 un-changes.
#system("cp Seg/UterusM_1_segments.txt Seg.new/UterusM_1_segments.txt")
#system("cp Seg/UterusM_3_segments.txt Seg.new/UterusM_3_segments.txt")


#folder = "/data1/qingjian/Rproject/Three/medicc/Seg.new"
#segfiles = list.files( folder, pattern = ".txt", full.names = T)
#sampleid = list.files( folder, pattern = ".txt") %>% stringr::str_remove("_segments.txt")

#for panCancer, combined all cancer into one types.
#data = red.seg(segfiles = segfiles, 
#               sampleid = sampleid , out.dir = "panCancer-New", tumor = "panCancer",
#               N.baf = 30, cnv_min_length = 5e5, max_CNt = 8, #quality control for cnv calling
#               minLength = 1e6, maxCNV = 4 #output controls.
#)

#We can construct Breast-Lung-UterusM firstly, then merge them with COAD.



folder = "/data1/qingjian/Rproject/Three/medicc/Seg.new/"
segfiles = list.files( folder, pattern = ".txt", full.names = T)
sampleid = list.files( folder, pattern = ".txt") %>% stringr::str_remove("_segments.txt")

#for panCancer, combined all cancer into one types.
# Met1

sampleid =c(
  paste0("Coad_", 1:5), paste0("OveryLM_", 1:5), paste0("OveryRM_", 1:6), paste0("UterusM_", c(1, 3) )
)

data = red.seg(segfiles = segfiles, 
               sampleid = sampleid , 
               out.dir = "panCancer", tumor = "Pan",
               N.baf = 30, cnv_min_length = 5e5, max_CNt = 8, #quality control for cnv calling
               minLength = 1e6, maxCNV = 4 #output controls.
)



# for the cnp2cpn analysis ------------------------------------------------

#Example of fasta file: 
#  > my_first_cnp
#1,2,0,3,2,1
#> my_second_cnp
#1,2,3,0,1,1
# 
# #combind data
# red.seg =  function(segfiles, sampleid, out.dir = "data",  tumor = "tumor", minLength = 1e5,  N.baf =30   ){
#   
#   seglist = list()
#   #read segs
#   for(i in 1:length(sampleid[1:3])){
#     
#     seg = read.delim(file = segfiles[i], header = T, stringsAsFactors = F) %>%
#       #dplyr::filter(chromosome == "chr1") %>%
#       dplyr::mutate(sample =sampleid[i] )
#     
#     seglist[[i]] = seg
#   }
#   
#   seglist = base::Reduce(rbind, seglist) %>%
#     dplyr::filter(N.BAF >= N.baf)
#   
#   gseg = GRanges(seqnames = seglist$chromosome,
#                  ranges = IRanges(seglist$start.pos, seglist$end.pos),
#                  strand = "+"
#                  
#   )
#   
#   #metadata columns can be added to a GRanges object
#   mcols(gseg) = seglist
#   
#   
#   #split regions into small regions.
#   segdis =  GenomicRanges::disjoin(gseg)
#   #add region infor to data
#   mcols(segdis) = data.frame(segdis)
#   
#   #set the minimal site of each segs.
#   #minLength = 1e5
#   
#   #get overlaps regions
#   overlaps =  findOverlaps(segdis, gseg )
#   
#   #combined information
#   merge = cbind( mcols(segdis[queryHits(overlaps) , ]) , mcols(gseg[subjectHits(overlaps) ,] ) )  %>%
#     as.data.frame() %>%
#     #filter seg length
#     dplyr::filter( width >= minLength)
#   
#   
#   #CNt info: 
#   merge_CNt = merge %>%
#     dplyr::select(seqnames, start, end, width, CNt, sample) %>%
#     tidyr::spread(key = sample, value = CNt, fill = 2 )
# 
#   return( list(data = merge, merge_CNt = merge_CNt))
# 
# }
# 
# 
# #write into fasta
# write.fasta = function( merge_CNt , out.dir = "data", tumor = "tumor"){
#   system(sprintf("mkdir %s", out.dir))
#   
#   #filter chrX and chrY.  
#   merge_CNt = dplyr::filter( !(seqnames %in% c("chrX", "chrY")))
#   
#   
#     num = ncol(merge_CNt)
#     
#     text = rep("0", 2*(num - 3))
#     text[1] = ">diploid"
#     text[2] = paste0( rep(1, nrow(merge_B)), collapse = "")
#     
#     for(j in 5:num){
#       text[ 2*(j-3) -1 ] = sprintf(">%s", colnames(merge_B)[j] )
#       text[ 2*(j-3)] = merge_B[,j] %>% stringr::str_c(collapse = "")
#     }
#   write.table(text, file = sprintf("%s/%s_%s_%s.fasta", out.dir, tumor ,major, i), quote = F, col.names = F, row.names = F )
#  
# }

