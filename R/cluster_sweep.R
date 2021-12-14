#' Performs a hyperparameter sweep across lambda, nclust and quantile function for a dataset
#' @param X An n x p matrix with n cells and p features
#' @param coord An n x 2 matrix with X and Y coordinates for each cell
#' @param nclust A vector of numbers of clusters to sweep
#' @param lambda A vector of lambda values to sweep
#' @param metrics A named list of quantile functions to sweep
#' @param max.iter Maximum iterations per sweep
#' @returns A 'long' dataframe with cluster membership for each cell along with sweep information
#' @export
cluster.sweep = function(X, coord, nclust = c(3,6,10),
                         lambda = c(1),
                         metrics = list(colMins=matrixStats::colMins, colMeds=matrixStats::colMedians, colMaxs=matrixStats::colMaxs),
                         max.iter = 100) {
  N = dim(X)[1]
  grid = expand.grid(nclust,lambda,metrics)
  num.exp = dim(grid)[1]
  df = matrix(0,nrow = num.exp*N,ncol = 9)
  cell.dists = as.matrix(stats::dist(coord))
  # make all cell dists on [0,1]
  cell.dists = cell.dists / max(cell.dists)
  # normalize features
  X = scale(X)
  for(j in 1:num.exp) {
    i = j - 1
    par = grid[j,]
    cat(par[[1]],par[[2]],names(par[[3]]),'\n')
    d = k.means.spatial(X,cell.dists,par[[1]],par[[2]],par[[3]],max.iter)
    d_normal = kmeans(X,par[[1]],iter.max = 50)$cluster
    d_hier = stats::dist(X)
    fit<-hclust(d_hier, method="ward.D2") # Ward method
    hier_res = cutree(fit, par[[1]])
    r = (i*N + 1):((i+1)*N)
    df[r,3] = as.vector(d$est)
    df[r,1:2] = as.matrix(coord)
    df[r,4] = par[[1]]
    df[r,5] = par[[2]]
    df[r,6] = names(par[[3]])
    df[r,7] = d$num_it
    df[r,8] = d_normal
    df[r,9] = hier_res
  }
  df = as.data.frame(df)
  colnames(df) = c("x","y","cluster","nclust","lambda","quantile_fun","num_it","cluster_kmeans","cluster_hier")
  df$cluster = factor(df$cluster)
  df$cluster_kmeans = factor(df$cluster_kmeans)
  df$cluster_hier = factor(df$cluster_hier)
  df$x = as.numeric(df$x)
  df$y = as.numeric(df$y)
  df$y = as.numeric(df$y)
  df$nclust = as.numeric(df$nclust)
  df$lambda = as.numeric(df$lambda)
  df$quantile_fun = factor(df$quantile_fun)
  df$num_it = as.numeric(df$num_it)
  return(df)
}

