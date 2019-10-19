#' Pre-processing of microarray datasets
#'
#' This helper function pre-processes microarray datasets by performing an exponentiation with number 2 as the base
#' on the expression values.
#'
#' @param expr A numeric vector with gene expression values
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
#'
#' @return Returns a numeric vector with the exponentiated expression values.
#'
#' @author David Lauenstein
#'
#' @examples
#' # Example based on the heat-shock dataset
#' data(heat)
#' heat = as.matrix(heat)
#' # Run the findPeaks function for the first gene in the expression matrix
#' peaks <- maProcessing(heat[1,], heat)
#'
#' @export

maProcessing <- function(expr, exprmat){
  normexprmat <- matrix(nrow = dim(exprmat)[1], ncol = dim(exprmat)[2])
  for (gene in 1:dim(exprmat)[1]){
    for (timepoint in 1:length(expr)){
      expr[timepoint] <- 2^(expr[timepoint])
    }
  }
  return(expr)
}
