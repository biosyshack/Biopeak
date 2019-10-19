#' Identification of co-expressing genes
#'
#' The getCormat function calculates a pair-wise correlation matrix and plots a bi-clustered heatmap.
#'
#' @import RColorBrewer
#' @import graphics
#' @import stats
#' @importFrom gplots heatmap.2
#'
#' @param peakdet A list returned by the peakDetection function.
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
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
#' corobjects <- getCormat(peakdet, heat, method = 'spearman')
#'
#' @export

getCormat <- function(peakdet, exprmat, method = 'spearman'){
  exprmat_act <- t(exprmat[!is.na(peakdet$activation),])
  cormat <- cor(exprmat_act, method=method)
  hm <- heatmap.2(cormat, col = rev(brewer.pal(9,"RdYlBu")), trace = 'none')
  hm_cormat <- cormat[rev(hm$rowInd), hm$colInd]
  return(list(hm = hm, hm_cormat = hm_cormat))
}
