
#######################
# Auxiliary functions #
#######################

plot.tree = function( output, mutdata, vaf.col.names, sample.groups = NULL, 
                      founding.cluster = 1, ignore.clusters = NULL, cluster.col.name = "cluster", 
                      weighted = FALSE, depth.col.names = NULL, plt.pairwise = F,
                      subclonal.test.model = "non-parametric", sum.p = 0.05, alpha = 0.05,
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
  
  y.B = infer.clonal.models(variants = mutdata,
                            cluster.col.name = cluster.col.name,
                            #ccf.col.names	= vaf.col.names,
                            vaf.col.names	= vaf.col.names, 
                            sample.groups = sample.groups,
                            cancer.initiation.model= "monoclonal",
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

  #Converting node-based trees to branch-based trees
  y.B <- convert.consensus.tree.clone.to.branch(y.B ,  cluster.col = cluster.col.name, branch.scale = "sqrt")
  
  
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
                     fancy.variant.boxplot.vaf.suffix = ".VAF",
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
                     scale.monoclonal.cell.frac = TRUE,
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



#########################################################################################
##test condes
if(FALSE){

vaf.col.names = vaf.col.names.sub[c(4,7)]
sample.groups = sample.groups.sub[c(4,7)]
founding.cluster = 1 
ignore.clusters = c(6, 7)
cluster.col.name = "cluster"
plt.pairwise = F
weighted = T
depth.col.names = str_c(vaf.col.names, ".depth")



y.B = infer.clonal.models(variants = mutdata,
                          cluster.col.name = cluster.col.name,
                          #ccf.col.names	= vaf.col.names,
                          vaf.col.names	= vaf.col.names, 
                          sample.groups = sample.groups,
                          cancer.initiation.model= "monoclonal",
                          subclonal.test = "bootstrap",
                          subclonal.test.model = "normal-truncated",
                          num.boots = 1000,
                          founding.cluster = founding.cluster,
                          cluster.center = "median",
                          ignore.clusters = ignore.clusters ,
                          clone.colors = clone.colors,
                          min.cluster.vaf = 0.01,
                          merge.similar.samples = F,
                          weighted = weighted,
                          depth.col.names = depth.col.names,
                          score.model.by = "metap",
                          # min probability that CCF(clone) is non-negative
                          sum.p = 0.05,
                          # alpha level in confidence interval estimate for CCF(clone)
                          alpha = 0.05)




#only plot bell shapes and clones cells

plot.clonal.models(y.merge,
                   
                   # box plot parameters
                   box.plot = FALSE,
                   fancy.boxplot = FALSE,
                   fancy.variant.boxplot.highlight = "is.driver",
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.color = "black",
                   fancy.variant.boxplot.highlight.note.col.name = "gene_site",
                   fancy.variant.boxplot.highlight.note.color = "blue",
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = "grey50",
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = ".VAF",
                   
                   # bell plot parameters
                   clone.shape = "bell", #polygon or bell
                   bell.event = FALSE,
                   bell.event.label.color = "blue",
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,

                   # node-based consensus tree parameters
                   merged.tree.plot = FALSE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = "circle",
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = FALSE,
                   mtcab.event.sep.char = ",",
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
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
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   
                   # output figure parameters
                   max.num.models.to.plot = 2,
                   out.dir = sprintf("%s/%s", output, "simple"),
                   out.format = "pdf",
                   overwrite.output = TRUE,
                   width = 6,
                   #height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(4,2),
                   
                   color.node.by.sample.group = T
                   
                   )



}



# Map calls: mutations to CNA segments with CNAqc
map_calls = function(CNA_calls, mutation_calls, samples, purities){
  mp = function(sample_id)
  {
    sample = samples[sample_id]
    purity = purities[sample_id]
  
    # Diploid segments with at least 500 SNVs
    CNA_calls = CNA_calls %>% select(chr, from, to, starts_with(sample))
    colnames(CNA_calls)[4:5] = c('minor', 'Major')
    
    SNV_calls = mutation_calls %>% select(chr, from, to, ref, alt, starts_with(sample), -ends_with('_N.VAF'))
    colnames(SNV_calls)[6:8] = c('DP', 'NV', 'VAF')
    
    # Use CNAqc to map mutations to segments
    init(snvs = SNV_calls, cna = CNA_calls, purity = purity)
  }
  
  CNAqc_objects = lapply(
    seq(samples),
    mp)
  names(CNAqc_objects) = samples
  
  return(CNAqc_objects)
}



# Plot CNA segments and VAF distributions 
plot_calls = function(CNAqc_objects){
  # Comparative CNA plot
  plot_cna = CNAqc::plot_multisample_CNA(CNAqc_objects)
  
  plot_snvs = lapply(names(CNAqc_objects),
                     function(x) {
                       CNAqc::plot_data_histogram(CNAqc_objects[[x]]) +
                         xlim(0.05, 1) +
                         labs(title = x)
                     }
)
  
  # Assemble figure: CNA (all) and diploid mutations (right)
  ggarrange(plot_cna,
            ggarrange(
              plotlist = plot_snvs,
              ncol = ceiling(length(plot_snvs)/2),
              nrow = 2
            ),
            ncol = 2)
}

# Run MOBSTER on a list of samples, using mutation data
fit_mobsters = function(mutations, samples){
  # Set the seed
  set.seed(123456)
  
  mobster_fits = lapply(
    samples,
    function(x)
    {
      SNV_calls = mutations %>%
        select(chr , from , to , ref, alt , ends_with(!!x)) 
      
      colnames(SNV_calls) = gsub(pattern = x, replacement = '', x = colnames(SNV_calls))
      colnames(SNV_calls) = gsub(pattern = '\\.', replacement = '', x = colnames(SNV_calls))
      
      if( any(SNV_calls$VAF >1)  ){
        stop("SNV_calls VAF greater than 1") 
      }
        
      mobster_fit(
        SNV_calls  %>% filter(VAF > 0.05 & VAF <1) ,
        description = x,
        parallel = FALSE,
        K = 1:2,
        init = 'peaks',
        model.selection = 'ICL',
        samples = 3,
        maxIter = 150,
        epsilon = 1e-6
      )$best
    })
  names(mobster_fits) = samples
  
  return(mobster_fits)
}

# Get ids for tail mutations
get_nontail_mutations = function(mutations, mobster_fits, filter = "any"){
  # Clusters, retain non-tail mutations
  non_tail_mutations = lapply(names(mobster_fits),
                              function(x)
                                mobster::Clusters(mobster_fits[[x]]) %>%
                                mutate(id = paste(chr, from, to, ref, alt, sep = ':')) %>%
                                select(id, cluster) %>%
                                mutate(sample = x))
  
  non_tail_mutations = Reduce(bind_rows, non_tail_mutations) %>%
    spread(id, cluster) %>%
    select(starts_with('chr'))
  
  # Use the mutation id
  ids = colnames(non_tail_mutations)
  if(filter == "any"){
    message("Filtering is ", filter)
    ids_left = ids[apply(non_tail_mutations, 2, function(x) all(x != 'Tail', na.rm = T))]
  }else{
    message("Filtering is ", filter)
    ids_left = ids[apply(non_tail_mutations, 2, function(x) !(all(x == 'Tail', na.rm = T)))]
  }
  
  mutations %>%
    mutate(id = paste(chr, from, to, ref, alt, sep = ':')) %>%
    filter(id %in% ids)
}

# Return a set of colors for VIBER clusters, using the wesanderson palettes
get_cluster_colors = function(viber_fit, W = 'M')
{
  # Get 5 nice colours from the 
  # colors = sapply(palettes, wesanderson::wes_palette) %>% as.vector
  
  require(ggsci)
  require("scales")
  
  if(W == "S")
  {
    # nej = pal_nejm("default")(4)
    # jco = (pal_jco("default")(3))[-1]
    # lzo = (pal_locuszoom("default")(4))[-c(1,2,3)]
    # igv = (pal_igv("default")(4))[-2]
    # 
    # # show_col(c(nej, jco, lzo, igv))
    # colors = c(nej, jco, lzo, igv)
    jco = (pal_jco("default")(10))
    colors = c(jco)
  }
  
  if(W == "M")
  {
    # str = pal_d3()(7)
    str = wesanderson::wes_palette("Darjeeling1", 5)

    # show_col(c(str))
    colors = c(str)
  }
    
    order_by_size = order(viber_fit$pi_k * viber_fit$N, decreasing = F) 
    colors = colors[order_by_size]
    
    non_zero_clusters = which((viber_fit$pi_k * viber_fit$N) %>% round > 0) %>% names
    names(colors) = non_zero_clusters
  
  colors
}

# Squared complex plot
squareplot = function(mobster_fits, viber_fit_bottom, viber_fit_top, samples_list, colors_bottom, colors_top){
  row_plots = NULL
  for (s in seq(samples_list))
  {
    sn = samples_list[s]
    mb = list(plot(mobster_fits[[sn]]) + labs(title = sn) )
    
    idx_pre = 1:s
    idx_post = s:length(samples_list)
    
    pl_r = pl_l = NULL
    
    if (length(idx_pre) > 1)
      pl_r = lapply(setdiff(idx_pre, s), function(x) {
        VIBER::plot_2D(viber_fit_bottom, d1 = sn, d2 = samples_list[x], colors = colors_bottom)
      })
    
    if (length(idx_post) > 1)
      pl_l = lapply(setdiff(idx_post, s), function(x) {
        VIBER::plot_2D(viber_fit_top, d1 = sn, d2 = samples_list[x], colors = colors_top)
      })
    
    row_plot = cowplot::plot_grid(
      plotlist = append(append(pl_r, mb), pl_l),
      nrow = 1,
      ncol = length(pl_r) + length(pl_l) + 1,
      align = 'h',
      axis = 'x'
    )
    
    row_plots = append(row_plots, list(row_plot))
  }
  
  cowplot::plot_grid(
    plotlist = row_plots,
    ncol = 1,
    nrow = length(row_plots),
    align = 'v'
  )
}



plot_sample_occurrences = function(x){
  clusters = x$labels %>% unique %>% unlist
  
  N = x$x %>% select(-cluster.Binomial)
  D = x$y %>% select(-cluster.Binomial)
  V = N/D
  
  V$counts = apply(V, 1, function(x) sum(x > 0))
  V$cluster = x$x$cluster.Binomial
  
  N = V %>% group_by(cluster) %>% summarise(N = n())
  N = pio:::nmfy(N$cluster, N$N)
  
  Np = N/sum(N)
  Np = round(Np, 2) * 100
  
  caption = paste0('n = ', N)
  caption = paste0(names(N), ', ', caption)
  caption = paste0(caption, ' [', Np, '%]')
  caption = paste0(caption, collapse = '; ')
  
  colors = RColorBrewer::brewer.pal(n = 9, name = 'Purples')
  colors[1] = 'steelblue'
  
  Vp = V %>% 
    group_by(cluster, counts) %>% 
    summarise(O = n()) %>%
    ungroup() %>%
    mutate(perc = O/N[cluster])
  
  ggplot(Vp,
         aes(x = cluster, y = perc, fill = paste(counts))
         )+
    geom_bar(stat = 'identity') +
    VIBER:::my_ggplot_theme() +
    scale_fill_manual(values = colors) +
    labs(
      title = "Mapping of cluster points to biopsies",
      subtitle = caption,
      x = "Cluster",
      y = "Percentage of points"
    ) +
    guides(fill = guide_legend("Number of biopsies with VAF > 0", nrow = 1))
  
  
}
  
# Plot cluster mapping among 2 anlyses
plot_mapping = function(mutations, mobster_fits, with_mobster, without_mobster){
  # Add ID to each mutation
  mutations = mutations %>%
    mutate(id = paste(chr, from, to, ref, alt, sep = ':'))
  
  # Standard analysis
  standard_analysis_assignments = 
    pio:::nmfy(
      mutations$id,
      without_mobster$labels %>% unlist()
    )
  
  non_tail_mutations = get_nontail_mutations(mutations, mobster_fits)
  
  with_mobster_analysis_assignments = 
    pio:::nmfy(
      non_tail_mutations$id,
      with_mobster$labels %>% unlist()
    )
  
  # All together
  all_assignments = data.frame(
    id = mutations$id,
    stringsAsFactors = FALSE
  )
  
  all_assignments$standard = standard_analysis_assignments[all_assignments$id]
  all_assignments$mobster = with_mobster_analysis_assignments[all_assignments$id]
  all_assignments$mobster[is.na(all_assignments$mobster)] = "Missing"
  
  all_assignments = all_assignments %>%
    group_by(standard, mobster) %>%
    summarise(N = n())
  
  colors_m = c(
    get_cluster_colors(W = "M", with_mobster),
    `Missing` = 'gray'
  )
  
  colors_s = c(
    get_cluster_colors(W = "S", without_mobster),
    `Missing` = 'gray'
  )
  
  p1 = ggplot(all_assignments, aes(x = mobster, y = N, fill = mobster)) +
    geom_bar(stat = 'identity') +
    mobster:::my_ggplot_theme() +
    facet_wrap(~standard, scales = 'free_y', nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors_m)
  
  p2 = ggplot(all_assignments, aes(x = standard, y = N, fill = standard)) +
    geom_bar(stat = 'identity') +
    mobster:::my_ggplot_theme() +
    facet_wrap(~mobster, scales = 'free_y', nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors_s)
  
  ggarrange(p1, p2, nrow = 2)
  
}

plot_dnds = function(x, mode, gene_list, dndscv_plot,colors, mask_colors = FALSE){
    x$globaldnds %>%
      dplyr::filter(name == 'wall') %>%
      ggplot2::ggplot(ggplot2::aes(
        x = name,
        y = mle,
        ymin = cilow,
        ymax = cihigh
      )) +
      mobster:::my_ggplot_theme() +
      facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
      ggplot2::xlab("") +
      ggplot2::ylab("dN/dS") +
      labs(
        title = paste0("dN/dS values via dndscv")
      ) +
      ggplot2::geom_hline(yintercept = 1.0,
                          lty = 2,
                          size = .3) +
      guides(color = F, fill = F) + 
      geom_pointrange(color = 'black')
    
    
  }

