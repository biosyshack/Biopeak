#' Identification of biomarkers specific to distinct phases of the underlying biological process
#'
#' The peakDetection function facilitates the identification of impulse-like gene expression changes based on user-defined
#' selection criteria. This function calls the helper functions: bgCorr(), maProcessing() and findPeaks().
#'
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
#' @param series A numeric vector defining the experimental series (e.g. time-points of sample acquisition).
#' @param actstrength Threshold for minimal activation relative to the mean expression across all time-points.
#' @param prominence Threshold for minimal peak prominence relative to the second highest peak.
#' @param type A character string defining the sequencing platform. Possible values are c('microarray', 'rnaseq').
#' @param minexpr An optional threshold for minimal mean expression across all time-points for a given gene.
#' @param peakwidth An optional definitino of the minimal number of time-points that a peak spans (based on sustact threshold).
#' @param sustact An optional threshold for minimal peakheight relative to the main peak to be considered as sustained activation.
#' @param bgcorr An optional logical constant (TRUE or FALSE) defining if a background noise correction is performed or not.

#' @return Returns a list comprising of multiple vectors and matrices. A numeric vector with the location of each peak (peakloc),
#' a numeric vector with the absolute height of each peak (peakheight), a character vector of gene symbols for which at least
#' one peak has been identified (peakgenes), a numeric matrix containing time-points with sustained activation, the logical
#' vector defining which gene index has been selected and the numeric input vector defining the time-series.
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
#'
#' @export


peakDetection <- function(exprmat, series, actstrength = 1.3, prominence = 1.3, type = 'rnaseq',
                          minexpr = 0, peakwidth = 0, sustact = 0.6, bgcorr = T){

  if (bgcorr == T){
    exprmat <- bgCorr(exprmat)
  }

  activation <- c()
  neighbors <- matrix(0,nrow = dim(exprmat)[1], ncol = dim(exprmat)[2])
  peakheight_out <- c()
  actcount <- 1

  for (gene in 1:dim(exprmat)[1]){
    expr <- exprmat[gene,]
    meanexpr <- mean(expr)

    if (type == 'microarray'){
      expr <- maProcessing(expr, exprmat)
    }

    rawexpr <- expr
    expr <- c(0, expr, 0) # set baseline of expression signal to zero

    peakdata <- findPeaks(expr)
    peakloc <- peakdata[[1]]
    peakheight <- peakdata[[2]]
    peakcount <- peakdata[[3]]

    if (length(peakheight) > 1){ # sort peaks by height
      peakheight_sort <- sort(peakheight, decreasing = T)
    }
    else{
      peakheight_sort <- peakheight
    }

    if (length(peakheight) > 0 & meanexpr > minexpr){# filter criteria: at least one peak and a minimal expression
      if (peakheight_sort[1] > mean(rawexpr)*actstrength){# exclude small peaks
        if (length(peakheight)>1){# in case of more than 1 peak
          if (peakheight_sort[1] > peakheight_sort[2]*prominence){# select for single large peaks
            activation[actcount] <- peakloc[which.max(peakheight)]
            peakheight_out[actcount] <- peakheight_sort[1]
            neighbors[gene,peakloc[which.max(peakheight)]-1] <- 1
            actcount <- actcount+1
            for (width in 1:length(series)){
              if (peakloc[which.max(peakheight)]-1-width>0){
                if (rawexpr[peakloc[which.max(peakheight)]-1-width]>peakheight_sort[1]*sustact){
                  #check for sustained peak activation in neighbourhood
                  neighbors[gene,peakloc[which.max(peakheight)]-1-width] <- width+1
                }
              }
              if  (peakloc[which.max(peakheight)]-1+width<=length(rawexpr)){
                if (rawexpr[peakloc[which.max(peakheight)]-1+width]>peakheight_sort[1]*sustact){
                  #check for sustained peak activation in neighbourhood
                  neighbors[gene,peakloc[which.max(peakheight)]-1+width] <- width+1
                }
              }
            }
          }
          else { # no single large peak
            activation[actcount] <- NA
            peakheight_out[actcount] <- NA
            actcount <- actcount+1
          }
        }
        else { # only 1 peak
          activation[actcount] <- peakloc[which.max(peakheight)]
          peakheight_out[actcount] <- peakheight_sort[1]
          neighbors[gene,peakloc[which.max(peakheight)]-1] <- 1
          actcount <- actcount+1
          for (width in 1:length(series)){
            if (peakloc[which.max(peakheight)]-1-width>0){
              if (rawexpr[peakloc[which.max(peakheight)]-1-width]>peakheight_sort[1]*sustact){
                #check for sustained peak activation in neighbourhood
                neighbors[gene,peakloc[which.max(peakheight)]-1-width] <- width+1
              }
            }
            if  (peakloc[which.max(peakheight)]-1+width<=length(rawexpr)){
              if (rawexpr[peakloc[which.max(peakheight)]-1+width]>peakheight_sort[1]*sustact){
                #check for sustained peak activation in neighbourhood
                neighbors[gene,peakloc[which.max(peakheight)]-1+width] <- width+1
              }
            }
          }
        }
      }
      else {# activation not met
        activation[actcount] <- NA
        peakheight_out[actcount] <- NA
        actcount <- actcount+1
      }
    }
    else { # meanexpr not met or no peak
      activation[actcount] <- NA
      peakheight_out[actcount] <- NA
      actcount <- actcount+1
    }
  }

  peakgenes <- rownames(exprmat)[which(is.na(activation) == FALSE)]
  peakgenes <- as.matrix(peakgenes)

  peakloc <- activation[which(is.na(activation) == FALSE)]
  peakloc <- series[peakloc-1]
  peakheight_out <- peakheight_out[which(is.na(peakheight_out) == FALSE)]

  rownames(neighbors) <- rownames(exprmat)
  neighbors <- neighbors[which(is.na(activation) == FALSE),]

  peakgenes_filtered <- c()
  peakloc_filtered <-  c()
  peakheight_filtered <- c()

  if (peakwidth > 0){ # filter for genes with sustained peak activation in neighborhood
    neighborgenes <- c()
    neighborfilter <- c()
    for (row in 1:dim(neighbors)[1]){
      mainpeak = which(neighbors[row,] == 1)
      neighborcount = 0
      for (width in 1:peakwidth){
        if (mainpeak+width <= length(series)){
          if (neighbors[row,mainpeak+width] == width+1){
            neighborcount <- neighborcount+1
          }
        }
        if (mainpeak-width > 0){
          if (neighbors[row,mainpeak-width] == width+1){
            neighborcount <- neighborcount+1
          }
        }
      }
      if (neighborcount >= peakwidth){
        neighborgenes <- append(neighborgenes,rownames(neighbors)[row])
        neighborfilter <- append(neighborfilter, row)
      }
    }
    peakgenes_filtered <- neighborgenes
    peakloc_filtered <- peakloc[neighborfilter]
    peakheight_filtered <- peakheight[neighborfilter]
    neighbors_filtered <- neighbors[neighborfilter,]
    activation_filtered <- activation[neighborfilter]
  }
  else {
    peakgenes_filtered <- peakgenes
    peakloc_filtered <- peakloc
    peakheight_filtered <- peakheight_out
    neighbors_filtered <- neighbors
    activation_filtered <- activation
  }

  results <- list(peakgenes = peakgenes_filtered,peakloc = peakloc_filtered, peakheight=peakheight_filtered,
                  neighbors=neighbors_filtered, activation = activation_filtered, series = series)
  return(results)
}
