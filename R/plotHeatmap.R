#' Plot a heatmap for selected genes
#'
#' This function acts as a wrapper function for the heatmap.2 function of the gplots package and normalizes the
#' subjected expression matrix to the log2 of the mean expression of the gene across all time-points.
#'
#' @import RColorBrewer
#' @importFrom gplots heatmap.2
#'
#' @param peakdet A list returned by the peakDetection function.
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
#' @param clustermembers An optional character vector defining genes to be selected.

#' @return This function does not return any value but generates a heatmap plot.
#'
#' @author David Lauenstein
#'
#' @examples
#' # Example based on the heat-shock dataset
#' data(heat)
#' heat = as.matrix(heat)
#' # Define series
#' series <- c(37,40,41,42,43)
#' # Run the peak detection algorithm
#' peakdet <- peakDetection(heat, series, type ='rnaseq', actstrength = 1.5,
#' prominence = 1.3, minexpr = 5000)
#' # cluster exploration using kmeans with a maximum of 4 clusters to be assigned
#' clusters <- findClusters(peakdet, heat, maxclusters = 4, method = 'kmeans')
#' # Plot the heatmap for one of the clusters returned by the findClusters function
#' heatmap <- plotHeatmap(peakdet, heat, clustermembers = clusters$clustermembers[[1]])
#'
#' @export

plotHeatmap <- function(peakdet, exprmat, clustermembers = c()){
  if (length(clustermembers>0)){
    expr <- exprmat[which(rownames(exprmat) %in% clustermembers),]
  }
  else {
    expr <- exprmat
  }
  for (i in 1:dim(expr)[1]){ # normalize to log2
    meanexpr <- mean(expr[i,])
    for (j in 1:dim(expr)[2]){
      if (expr[i,j]>0){
        expr[i,j] <- log(expr[i,j]/meanexpr,2)
      }
      else{
        expr[i,j] <- log((expr[i,j]+.Machine$double.eps)/meanexpr,2)
      }
    }
  }
  heatmap.2(expr,col=rev(brewer.pal(8,"RdYlBu")),breaks=seq(-1,1,0.25),trace='none', Rowv = T)
}

##updates:

expr = exprmat
for (i in 1:dim(expr)[1]){ # normalize to log2
  meanexpr <- mean(expr[i,])
  for (j in 1:dim(expr)[2]){
    if (expr[i,j]>0){
      expr[i,j] <- log(expr[i,j]/meanexpr,2)
    }
    else{
      expr[i,j] <- log((expr[i,j]+.Machine$double.eps)/meanexpr,2)
    }
  }
}


library(pheatmap)

data(heat)
heat = as.matrix(heat)
series <- c(37,40,41,42,43)
peakdet <- peakDetection(heat, series, type ='rnaseq', actstrength = 1.5,
                         prominence = 1.3, minexpr = 1000)


genes = peakdet$peakgenes[which(peakdet$peakloc == 40 | peakdet$peakloc == 43)]

genes = sub("//.*", "", genes)

cormat = cor(t(heat[which(rownames(heat) %in% genes),]))

exprmat_sel = heat[which(rownames(heat) %in% genes),]

anno = data.frame(peakdet$peakloc[which(peakdet$peakloc == 40 |
                                          peakdet$peakloc == 43)])

rownames(anno) = genes
colnames(anno) = "Peak Loc: Temp"
anno$`Peak Location: Temp` = as.factor(anno$`Peak Location: Temp`)

pheatmap(cormat, annotation_row = anno, border_color = F)

pheatmap(scale(exprmat_sel),annotation_row = anno, border_color = T)


