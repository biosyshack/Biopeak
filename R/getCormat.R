#' Identification of co-expressing genes
#'
#' The getCormat function calculates a pair-wise correlation matrix and plots a bi-clustered heatmap.
#'
#' @import Pheatmap
#'
#' @param peakdet A list returned by the peakDetection function.
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
#' @param genes A character vector with genes to be selected. If set to NULL, the function will select all genes present in the original expression matrix.
#' @param legendvar A character string defining the legend title for the annotation column.
#' @param method A character string defining the correlation algorithm. Options are: c('pearson', 'kendall', 'spearman').
#'
#' @return Returns both the heatmap object and the re-ordered correlation matrix:
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
#' # calculate and plot correlation matrix
#' corobjects <- getCormat(peakdet, heat,
#' genes = peakdet$peakgenes[which(peakdet$peakloc == 43)],
#' legendvar = "Peak Loc: Temp", method = 'spearman')
#'
#' @export

getCormat <- function(peakdet, exprmat, genes = NULL, legendvar = "Peak Loc", method = 'spearman'){
  exprmat_act <- t(exprmat[!is.na(peakdet$activation),])
  if (is.null(genes)){
    genes = rownames(exprmat)
  }
  cormat <- cor(t(exprmat[which(rownames(exprmat) %in% genes),]), method=method)

  anno <- data.frame(peakdet$peakloc[which(peakdet$peakgenes %in% genes)])
  rownames(anno) <- genes[which(genes %in% peakdet$peakgenes)]
  colnames(anno) <- legendvar
  anno[,1] <- as.factor(anno[,1])

  hm <- pheatmap(cormat, annotation_row = anno, border_color = F)

  # reproduce re-ordered matrix from heatmap
  rowclust <- hclust(dist(cormat))
  reordered <- cormat[rowclust$order,]
  colclust <- hclust(dist(t(cormat)))
  hm_cormat <- reordered[, colclust$order]

  return(list(hm = hm, hm_cormat = hm_cormat))
}
