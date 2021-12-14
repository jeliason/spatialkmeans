#' Calculates the spatial distance between cells and each cluster. Distance to cluster is determined by metric.
#' @param cell.dists An n x n matrix of pairwise physical distances between all cells
#' @param clust A vector length n indiacting cluster membership for each cell
#' @param metric A column-wise quantile function to determine distance to cluster (eg. colMins, colMedians, colMaxs from matrixStats)
#' @return An nclust X n matrix with the physical distance to each cluster for each cell
#' @export
spatial.distance = function(cell.dists, clust, metric=matrixStats::colMins) {
  n = length(clust)
  grouped = sapply(1:max(clust), function(i) {
    mat = cell.dists[which(clust == i),]
    if(is.matrix(mat)) {
      return(list(metric(mat)))
    } else {
      return(list(mat))
    }
  })
  d = do.call(rbind,grouped)
  d
}
#' Clusters spatial omics data with spatial regularization
#' @param X An n x p matrix with n cells and p features
#' @param cell.dists An n x n matrix with pairwise distances between all cells
#' @param nclust A scalar indicating the number of clusters
#' @param lambda The weight that spatial penalty takes
#' @param metric A column-wise quantile function to determine distance to cluster (eg. colMins, colMedians, colMaxs from matrixStats)
#' @param max.iter Maximum number of iterations before stopping clustering
#' @returns a list containing mus (the centroid of each cluster), est (cluster membership for each cell), num_it (number of iterations taken to converge)
#' @export
k.means.spatial = function(X, cell.dists, nclust=3, lambda=0.1, metric=matrixStats::colMins, max.iter=100) {
  it = 0
  N = dim(X)[1]
  assign.changed = TRUE
  mus = t(X[sample(1:N,nclust),])
  last.clust = rep(1:nclust,length.out = N)
  clust = last.clust
  scale = max(colMeans(abs(X)))
  while(it < max.iter & assign.changed) {
    dists = outer(1:nclust,1:nrow(X), Vectorize(function(i,j) {
      sum((mus[,i] - X[j,])^2)
    }))
    phys.dists = spatial.distance(cell.dists, clust)
    phys.dists = lambda*scale*phys.dists
    dists = dists + phys.dists
    clust = apply(dists, 2, function(d) {
      which.min(d)
    })
    mus = sapply(1:nclust, function(i) {
      X.clust = X[clust == i,]
      if(is.vector(X.clust)) {
        return(X.clust)
      } else {
        return(colMeans(X.clust))
      }
    })
    if(sum(last.clust == clust) == nrow(X)) {
      assign.changed = FALSE
    }
    last.clust = clust
    it = it + 1
    if(it %% 10 == 0) {
      print(it)
    }
  }
  return(list(mus=mus,est=factor(clust),num_it=it))
}
