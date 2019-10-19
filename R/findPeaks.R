#' Identification of peaks in an expression signal
#'
#' This helper function identifies peaks in an expression signal by treating the gene expression as a signal that
#' propagates along an experimental axis. A peak is defined as a local maximum in the expression signal satisfying:
#' y(t) > y(t+1) and y(t) > y(t-1), where y(t) represents the gene expression as a function of series condition t.
#'
#' @param expr A numeric vector with gene expression values
#'
#' @return Returns a list comprising of a numeric vector with the location of each peak (peakloc), a numeric vector
#' with the absolute height of each peak (peakheight) and a character vector of gene symbols for which at least one
#' peak has been identified (peakgenes).
#'
#' @author David Lauenstein
#'
#' @examples
#' # Example based on the heat-shock dataset
#' data(heat)
#' heat = as.matrix(heat)
#' # Run the findPeaks function for the first gene in the expression matrix
#' peaks <- findPeaks(heat[1,])
#'
#' @export

findPeaks <- function(expr){
  peakloc <- c()
  peakheight <- c()
  peakcount <- 1
  for (i in 2:length(expr)){
    if (i<length(expr)){
      if (expr[i] > expr[i-1] & expr[i] > expr[i+1]){ # find local maxima
        peakloc[peakcount] <- i
        peakheight[peakcount] <- expr[i]
        peakcount <- peakcount+1
      }
    }
  }
  return(list(peakloc, peakheight, peakcount))
}
