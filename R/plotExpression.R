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

colpal = "Blues"

exprmat = heat
peakdet <- peakDetection(heat, series, type ='rnaseq', actstrength = 1.5, prominence = 1.3, minexpr = 400)


genes = peakdet$peakgenes[which(peakdet$peakloc == 40)]

plotExpression2 <- function(exprmat, genes, series, peakdet, colpal = "Dark2"){

}

exprmat_sel = c()
for (gene in genes){
  geneseries = data.frame(cbind(cbind(series,exprmat[which(tolower(rownames(exprmat)) == tolower(gene)),])))
  geneseries= cbind(geneseries,rep(gene,dim(exprmat)[2]))
  colnames(geneseries) = c("Series", "Expression", "Gene")
  exprmat_sel = rbind(exprmat_sel,geneseries)
}


ggplot(data=exprmat_sel, aes(x = Series, y = Expression)) +
  #geom_point()+
  geom_area(aes(color = Gene, fill = Gene), alpha = 0.8)+
  #geom_line(aes(y = Expression), position = "stack")+
  theme_minimal()+
  scale_color_brewer(palette=colpal)+
  scale_fill_brewer(palette=colpal)


# ggplot(data=exprmat_sel, aes(x = Series, y = Expression)) +
#   geom_point(aes(col = Gene))+
#   geom_line(aes(col = Gene))+
#   theme_minimal()+
#   scale_color_brewer(palette=colpal)

