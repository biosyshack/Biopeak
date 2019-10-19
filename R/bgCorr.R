#' Background noise correction of expression data
#'
#' This helper function performs a background noise correction before subjecting the corrected matrix to the peakDetection
#' function. Genes with an expression lower than the 5 % quantile of the entire expression matrix are discarded.
#'
#' @import stats
#'
#' @param exprmat A numeric matrix with time-series expression data with variables as rownames.
#'
#' @return Returns background noise corrected expression matrix.
#'
#' @author David Lauenstein
#'
#' @examples
#' # Example based on the heat-shock dataset
#' data(heat)
#' heat = as.matrix(heat)
#' # Execute the bgCorr function
#' exprmat_corrected <- bgCorr(heat)
#'
#' @export

bgCorr <- function(exprmat){ # helper function performing background correction
  estbg <- quantile(as.matrix(exprmat),probs = seq(0, 1, .05))[[2]] # estimate background noise as 5% quantile of expression
  exprmat[exprmat < estbg] <- NA
  exprmat <- na.omit(exprmat)
  return(exprmat)
}
