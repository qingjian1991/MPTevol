#' set.colors
#'
#' Some colors could be used.
#'
#' @param n number of colors to select. A total of 36 colors are saved.
#' @param rev reverse output the colors.
#'
#' @export


set.colors = function( n = 36, rev = FALSE){

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


  if(n <= 36){
     if(rev){colors = colorScale[36:(36-n+1)]}else{colors = colorScale[1:n]}
  }else{
     if(rev){colors = rev(colorScale)}else{colors = colorScale}
  }

  colors

}








