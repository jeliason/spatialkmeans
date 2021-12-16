# Implementation of Spatial k-means and Comparison with Existing Clustering Algorithms
This R package consists of implementation of the spatial k-means algorithm and its comparison with different algorithms. The implementation of spatial k-means incorporates the idea of "physical distance" into the distance metric - in our case, physical distance is incorporated as a 'spatialyl-regularizing' penalty term. This will generate a different clustering result depending on the penalty parameter $\lambda$. 

The comparison of different algorithms focuses on the different results obtained when changing the algorithms, the parameters each algorithm uses, and the metrics used for benchmarking - see the file "Comparison_with_visualizations.Rmd" for examples.

After installing from Github, the spatial k-means algorithm can be used similarly to this example:

    library(spatialkmeans)

    X = read.csv("../X_pca.csv", header = F)
    coord = read.csv("../10x_spatial_positions.csv", header = F)

    cell.dists = as.matrix(stats::dist(coord))
    cell.dists = cell.dists / max(cell.dists)

    res = k.means.spatial(X,cell.dists)

This library also contains a function `cluster_sweep()`, which will automatically perform a parameter sweep over the main parameters of `lambda`,`nclust` and `metric` (a quantile function for distance calculation).
