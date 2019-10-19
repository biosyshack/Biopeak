#' Plot the expression signal of individual genes
#'
#' This function plots the expression signal of a defined gene and marks the main peak location with a dashed line.
#'
#' @import graphics
#'
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
#' @param gene A character string (not case-sensitive) defining the gene to be plotted.
#' @param series A numeric vector defining the experimental series (e.g. time-points of sample acquisition).
#' @param peakdet A list returned by the peakDetection function.
#'
#' @return This function does not return any value but generates a plot.
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
#' # Plot the expression signal of the gene CXCL5
#' plotExpression(heat, 'CXCL5', series, peakdet)
#'
#' @export

plotExpression <- function(exprmat, gene, series, peakdet){
  plot(series, exprmat[which(tolower(rownames(exprmat)) == tolower(gene)),], type = 'l',
       main = gene, xlab = 'Time', ylab = 'Expression')
  clip(0,max(series),0,max(exprmat[which(tolower(rownames(exprmat)) == tolower(gene)),]))
  abline(v =  peakdet$peakloc[which(tolower(peakdet$peakgenes) == tolower(gene))], lty = 2)
}
