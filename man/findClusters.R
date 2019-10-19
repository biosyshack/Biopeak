#' Identification of clusters with similar temporal regulation
#'
#' The findClusters function estimates the number of genes with similar temporal regulation and supports three different
#' clustering algorithms: kmeans, dbscan and hierarchical clustering. Clustering is based on a PCA projection of the input
#' data.
#'
#' @import cluster
#' @import factoextra
#' @import dbscan
#'
#' @param peakdet A list returned by the peakDetection function.
#' @param exprmat A numeric matrix with expression series data with variables as rownames.
#' @param maxclusters Maximal number of clusters used for kmeans cluster estimation.
#' @param eps Epsilon value used by the dbscan algorithm.
#' @param clusters Number of clusters used for the cutree function of the hierarchical clustering.
#' @param method A character string defining the clustering algorithm with options: c('kmeans', 'dbscan', 'hclust').
#'
#' @return Returns a cluster assignment of each variable and the number of identified clusters.
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
#' # Cluster exploration using kmeans with a maximum of 4 clusters to be assigned
#' clusters <- findClusters(peakdet, heat, maxclusters = 4, method = 'kmeans')
#'
#' @export

findClusters <- function(peakdet, exprmat, maxclusters = 3, eps = 0.02, clusters = 3,method = 'kmeans'){
  pca <- prcomp(t(exprmat[!is.na(peakdet$activation),]), center=T, scale=T)

  if (method == 'kmeans'){
    gap <- clusGap(pca$rotation[,1:2], FUNcluster = hcut,K.max = maxclusters, B = 2)
    k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "Tibs2001SEmax")
    cl <- kmeans(pca$rotation[,1:2],centers=k,nstart = 20)
    clusterAS <- cl$cluster
  }

  else if(method == 'dbscan'){
    cl <- dbscan(pca$rotation[,1:2], eps = eps)
    clusterAS <- cl$cluster
  }

  else if(method == 'hclust'){
    hclust <- hclust(dist(exprmat[!is.na(peakdet$activation),]))
    clusterAS <- cutree(hclust, clusters)
  }

  pci <- data.frame(pca$rotation[,1:2],cluster = clusterAS)
  plot(pci$PC1, pci$PC2, col = pci$cluster, main = paste('Biomarker clusters (',method,')', sep=''), xlab = 'PC1', ylab = 'PC2')
  k <- max(clusterAS)

  if (k > 0){
    clusterlist = list()
    for (i in 1:k){
      clusterlist[[i]] <- peakdet$peakgenes[which(clusterAS == i)]
    }

    results <- list(clusters = clusterAS, k = k, clustermembers = clusterlist)
  }
  else{
    results <- list(clusters = clusterAS, k = k)
  }
  return(results)
}
