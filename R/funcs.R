
#' getGseaRes
#' for heatmapPathway app
#' @param input shiny input
#'
#' @return GseaRes
#' @export
#'
getGseaRes=function(input){
  phenotype=input$phenotype
  gseaRes=gseaRess[[phenotype]]
  gseaRes=gseaRes[gseaRes[[input$pthwySig]]<0.05,]
  return(gseaRes)
}

#' getPthwy
#' for heatmapPathway app
#' @param gseaRes global variable gseaRes
#'
#' @return Pthwy
#' @export
#'
getPthwy=function(gseaRes){
  pthwy=pathway[gseaRes$pathway]
  pthwy=lapply(pthwy, function(x) x[x %in% rownames(log2_expressions)])
  return(pthwy)
}

#' getGenes
#' for heatmapPathway app
#'
#' @param pthwy global variable pthwy (res of getPthwy)
#' @param input
#'
#' @return genes
#' @export
#'
getGenes=function(pthwy,input){
  genes=pthwy[[input$pathwyName]]
  genes=genes[genes%in%rownames(log2_expressions)]
  return(genes)
}







