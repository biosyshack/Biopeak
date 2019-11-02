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

plotExpression(peakdet, heat, genes = c("CXCL1","CXCL5","CXCL3"),
               series, xlab = "Temperature", type = "area")

plotExpression <- function(peakdet, exprmat, genes, series, type = "area",
                           colpal = "Blues", xlabel = "Series",
                           ylabel = "Expression",linewidth = 0.8,
                           alpha = 0.8){
  exprmat_sel = c()
  for (gene in genes){
    geneseries <- data.frame(cbind(cbind(series,
                                        exprmat[which(tolower(rownames(exprmat))
                                                      == tolower(gene)),])))
    geneseries= cbind(geneseries,rep(gene,dim(exprmat)[2]))
    colnames(geneseries) <- c("Series", "Expression", "Gene")
    exprmat_sel <- rbind(exprmat_sel,geneseries)
  }

  if (type == "area"){
    ggplot(data=exprmat_sel, aes(x = Series, y = Expression)) +
      geom_ribbon(aes(ymin = 0, ymax = Expression,
                    group = Gene, fill = Gene), alpha = alpha)+
      theme_minimal()+
      scale_color_brewer(palette=colpal)+
      scale_fill_brewer(palette=colpal)+
      labs(x = xlabel, y = ylabel)
  }

  else if (type == "line"){
    ggplot(data=exprmat_sel, aes(x = Series, y = Expression)) +
      geom_point(aes(col = Gene))+
      geom_line(aes(col = Gene), size = linewidth)+
      theme_minimal()+
      scale_color_brewer(palette=colpal)+
      labs(x = xlabel, y = ylabel)
  }
}
