#' Save the peak detection output to a text file
#'
#' This function saves the output of the peakDetetection funtion (peakgenes, peaklocation and peakheight)
#' to a text file.
#'
#' @import utils
#'
#' @param peakdet A list returned by the peakDetection function.
#' @param filename A character string defining the output file.
#'
#' @return This function does not return any value but saves data to a text file.
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
#' # Save the peak detection output to a text file
#' saveOutput(peakdet,' file.path(tempdir(),'heat_out.txt'))

#'
#' @export

saveOutput <- function(peakdet, filename){
  outmat <- cbind(peakdet$peakgenes, peakdet$peakloc, peakdet$peakheight)
  colnames(outmat) <- c('Gene Name', 'Peak Location', 'Peak Height')
  write.table(outmat, filename, sep = '\t', quote = F, row.names = F, col.names = T)
}
