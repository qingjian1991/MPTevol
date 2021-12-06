
#estimate.dist.gene
#' Estimate gene's distances.
#' @param method: methods to estimate genes. the distance measure to be used. This must be one of "pairwise", "percentage" or "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.  

estimate.dist.gene = function(mut_dat, method = "pairwise" ){
  
  if(method %in% c("pairwise","percentage") ){
    dist.gene(mut_dat, method = method)
  }else if(method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
    dist(m, method= method)
  }
}


#MP 
byMP <- function(mut_dat){
  #matTree <- nj(dist.gene(mut_dat))
  estimate_tr <- function(m) nj(dist(m, method="manhattan"))
  matTree <- estimate_tr(mut_dat)
  tree_dat <- phangorn::as.phyDat(mut_dat, type="USER", levels = c(0, 1))
  tree_pars <- suppressMessages(phangorn::optim.parsimony(matTree, tree_dat, method = "sankoff", trace = FALSE)) 
  matTree <- phangorn::acctran(tree_pars, tree_dat)
  return(matTree)
}

#ML
byML <- function(mut_dat){
  matTree <- nj(dist.gene(mut_dat))
  tree_dat <- phangorn::as.phyDat(mut_dat, type="USER", levels = c(0, 1))
  fitJC <- phangorn::pml(matTree, tree_dat)
  fitJC <- try(phangorn::optim.pml(fitJC,control = phangorn::pml.control(trace = FALSE))) 
  matTree <- fitJC$tree
  return(matTree)
}

getTrees = function(binary.matrix, method = "NJ", bootstrap.rep.num = 100){
  
  binary.matrix$NORMAL = 0
  
  mut_dat <- t(binary.matrix)
  
  if(method == "NJ"){
    matTree <- nj(dist.gene(mut_dat))
    root_num <- which(matTree$tip.label == "NORMAL")
    matTree <- root(matTree, root_num)
    bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e){nj(dist.gene(e))}, B = bootstrap.rep.num,quiet = TRUE, rooted = TRUE)/(bootstrap.rep.num)*100
  }else if(method == "MP"){
    matTree <- byMP(mut_dat)
    root_num <- which(matTree$tip.label == "NORMAL")
    matTree <- root(matTree, root_num)
    bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e){byMP(e)}, B = bootstrap.rep.num, quiet = TRUE, rooted = TRUE, mc.cores = 10)/(bootstrap.rep.num)*100 
  }else if(method == "ML"){
    matTree <- byML(mut_dat)
    root_num <- which(matTree$tip.label == "NORMAL")
    matTree <- root(matTree, root_num)
    bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e)byML(e), B = bootstrap.rep.num, quiet = TRUE, rooted = TRUE,  mc.cores = 10)/(bootstrap.rep.num)*100
  }
  
  list(tree = matTree,
       bootstrap.value = bootstrap.value,
       method = method,
       binary.matrix = binary.matrix
       )
}




# set certain colors
colorScale <- c("#3C5488FF", "#00A087FF", "#F39B7fFF",
                "#8491B4FF","#E64B35FF","#4DBBD5FF",
                "#E41A1C", "#377EB8", "#7F0000",
                "#35978f", "#FC8D62", "#2166ac",
                "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#8DD3C7", "#6E016B" ,
                "#BEBADA", "#e08214", "#80B1D3",
                "#d6604d", "#ffff99", "#FCCDE5",
                "#FF6A5A", "#BC80BD", "#CCEBC5" ,
                "#fb9a99", "#B6646A", "#9F994E", 
                "#7570B3" , "#c51b7d" ,"#66A61E" ,
                "#E6AB02" , "#003c30", "#666666")


#Plot trees from clonevol
plot.tree = function( output, mutdata, vaf.col.names = NULL, ccf.col.names =NULL ,sample.groups = NULL, 
                      founding.cluster = 1, ignore.clusters = NULL, cluster.col.name = "cluster", 
                      weighted = FALSE, depth.col.names = NULL, plt.pairwise = F,
                      scale.monoclonal.cell.frac = TRUE,
                      subclonal.test.model = "non-parametric", cancer.initiation.model = "monoclonal",
                      sum.p = 0.05, alpha = 0.05,
                      highlight.note.col.name = "gene", highlight = "is.driver", highlight.CCF = FALSE ){
  
  if(plt.pairwise == TRUE){
    plot.pairwise(mutdata, col.names = vaf.col.names,
                  group.col.name = cluster.col.name,
                  onePage = FALSE,
                  multiPages = TRUE,
                  out.prefix = sprintf("%s.pair", output),
                  yMaxSmall = 100, xMaxSmall = 120,
                  colors = clone.colors)
  }
  message("infer.clonal.models")
  
  y.B = infer.clonal.models(variants = mutdata,
                            cluster.col.name = cluster.col.name,
                            ccf.col.names	= ccf.col.names,
                            vaf.col.names	= vaf.col.names, 
                            sample.groups = sample.groups,
                            cancer.initiation.model= cancer.initiation.model,
                            subclonal.test = "bootstrap",
                            subclonal.test.model = subclonal.test.model,
                            num.boots = 1000,
                            founding.cluster = founding.cluster,
                            cluster.center = "median",
                            ignore.clusters = ignore.clusters ,
                            clone.colors = clone.colors,
                            min.cluster.vaf = 0.01,
                            merge.similar.samples = F,
                            weighted = weighted,
                            depth.col.names = depth.col.names,
                            # min probability that CCF(clone) is non-negative
                            sum.p = sum.p,
                            # alpha level in confidence interval estimate for CCF(clone)
                            alpha = alpha)
  
  message("convert.consensus.tree.clone.to.branch")
  #Converting node-based trees to branch-based trees
  y.B <- convert.consensus.tree.clone.to.branch(y.B ,  cluster.col = cluster.col.name, branch.scale = "sqrt")
  
  message("plot.clonal.models")
  plot.clonal.models(y.B,
                     # box plot parameters
                     box.plot = TRUE,
                     fancy.boxplot = TRUE,
                     fancy.variant.boxplot.highlight = "is.driver",
                     fancy.variant.boxplot.highlight.size = 2.5,
                     fancy.variant.boxplot.highlight.shape = 21,
                     fancy.variant.boxplot.highlight.color = "blue",
                     fancy.variant.boxplot.highlight.fill.color = "green",
                     fancy.variant.boxplot.highlight.note.col.name = highlight.note.col.name,
                     fancy.variant.boxplot.highlight.note.color = "blue",
                     fancy.variant.boxplot.highlight.note.size = 2,
                     fancy.variant.boxplot.jitter.alpha = 1,
                     fancy.variant.boxplot.jitter.center.color = "grey50",
                     fancy.variant.boxplot.base_size = 12,
                     fancy.variant.boxplot.plot.margin = 1,
                     fancy.variant.boxplot.vaf.suffix = "",
                     fancy.variant.boxplot.founding.cluster = founding.cluster,
                     fancy.variant.boxplot.ccf = highlight.CCF,
                     
                     # bell plot parameters
                     clone.shape = "bell",
                     bell.event = TRUE,
                     bell.event.label.color = "blue",
                     bell.event.label.angle = 60,
                     clone.time.step.scale = 1,
                     bell.curve.step = 2,
                     # node-based consensus tree parameters
                     merged.tree.plot = TRUE,
                     tree.node.label.split.character = NULL,
                     tree.node.shape = "circle",
                     tree.node.size = 30,
                     tree.node.text.size = 0.5,
                     merged.tree.node.size.scale = 1.25,
                     merged.tree.node.text.size.scale = 2.5,
                     merged.tree.cell.frac.ci = FALSE,
                     # branch-based consensus tree parameters
                     merged.tree.clone.as.branch = TRUE,
                     mtcab.event.sep.char = ",",
                     mtcab.branch.text.size = 1,
                     mtcab.branch.width = 0.3,
                     mtcab.node.size = 3,
                     mtcab.node.label.size = 1,
                     mtcab.node.text.size = 1.5,
                     # cellular population parameters
                     cell.plot = TRUE,
                     num.cells = 100,
                     cell.border.size = 0.25,
                     cell.border.color = "black",
                     clone.grouping = "horizontal",
                     
                     #meta-parameters
                     scale.monoclonal.cell.frac = scale.monoclonal.cell.frac,
                     show.score = FALSE,
                     cell.frac.ci = TRUE,
                     disable.cell.frac = FALSE,
                     # output figure parameters
                     out.dir = sprintf("%s", output),
                     out.format = "pdf",
                     overwrite.output = TRUE,
                     width = 15,
                     #height = 4,
                     # vector of width scales for each panel from left to right
                     panel.widths = c(3,4,2,4,4))
  
  
  pdf(file = sprintf("%s/%s.trees.pdf", output, output), width = 6, height = 6)
  plot.all.trees.clone.as.branch(y.B, branch.width = 0.5,
                                 node.size = 2, node.label.size = 0.5,
                                 tree.rotation	=180
  )
  dev.off()
  
  return(y.B)
}







