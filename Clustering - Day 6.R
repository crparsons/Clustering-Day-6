# QMB 2016 - Clustering
# Nils Gehlenborg  
# August 22, 2016  

#We need the following R packages installed and loaded for this lecture.

```
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("dendextend")
install.packages("cluster")
install.packages("fpc")
install.packages("e1071")

library(gplots)
library(RColorBrewer)
library(dendextend)
library(cluster)
library(fpc)
library(e1071)
```

## Recap from Day 2

#- Clustering is an _unsupervised_ classification approach. _Unsupervised_ because the classes are not provided. Clustering is used to find groups of related objects within a set of objects. 
#- Clustering requires the following:
#  1. Set of objects to cluster.
#1. Method to compute similarity of (= distance between) two objects.
#1. Algorithm to determine cluster assignments based on pairwise similarities of all objects.
#- Common distance measures include Euclidean distance, Manhattan distance, or correlation measures.
#- Hierarchical clustering is a type of clustering that creates a hierarchy of clusters. It can either be _agglomerative_ (bottom up) or _divisive_ (top down).

## Goals for Today
#
#- Learn about different types of clustering
#- hierarchical, k-means, and probabilistic/fuzzy clustering
#- Evaluate clustering results
#- computational and visual methods
#- Pointers to more advanced topics
#- model-based clustering

## Hierarchical Clustering

#First, we retrieve a gene expression matrix for glioblastoma patients_.

library(curl)
gbm_data <- read.table(curl("https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/TCGA_unified_CORE_ClaNC840.txt"), sep="\t", row.names=1, header=1)
dim(gbm_data)

## [1] 841 174

#  This dataframe has been filtered to a subset of genes and samples but contains some meta data for the subtype assignments in row 1 and column 1.



#  Let's extract these into separate variables, trim them from the dataframe, and prepare the data for clustering:

gbm_subtype_genes <- factor(gbm_data[-1,1])
gbm_subtypes <- gbm_data[1,-1]

# trim the subtypes and genes/subtype information from the dataframe
gbm_data <- gbm_data[-1,-1]
```

#  Now we can use this matrix. First, we scale and center the columns of the matrix, i.e. we subtract the _mean_ and divide by the _standard deviation_ of each column. Then we compute the Euclidean distances between all patient samples, i.e. the columns of the matrix.


```r
# convert into numeric data frame
gbm_data <- apply(gbm_data,2,as.numeric)

gbm_data_scaled <- scale(gbm_data)
ed <- dist(t(gbm_data_scaled))

#  Then we can compute a hierarchical clustering of the patients and plot the dendrogram.


```r
hclusters <- hclust(ed)
plot(hclusters)
```

# ![](clustering_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

Color leaves of the dendrogram.


```r
subtype_colors <- c(Mesenchymal="red", Neural="green", Proneural="blue", Classical="orange");

dendrogram <- as.dendrogram(hclusters)
labels_colors(dendrogram) <- subtype_colors[as.matrix(gbm_subtypes)][order.dendrogram(dendrogram)]
plot(dendrogram)
```


par(mfrow=c(1,1))
dendrogram <- dendrogram %>% set("leaves_pch", 19) %>% set("leaves_cex", 1) %>% set("leaves_col", labels_colors(dendrogram))
plot(dendrogram)
```

#We will remove the sample identifier labels because they are not relevant for our purposes.

dendrogram <- dendrogram %>% set("labels",  rep( "", dim(gbm_data_scaled)[2] ))
plot(dendrogram)
```

# Let's compare different distance measures and linkage methods. First, we write a function to handle all the steps described so far and add parameters for the input matrix and class assignments as well as for distance and linkage methods.

cluster_comparison <- function( matrix, class.colors, linkage, distance ) {
  distances <- dist(t(matrix), method=distance)
  hclusters <- hclust(distances, method=linkage)
  subtype_colors <- class.colors;
  dendrogram <- as.dendrogram(hclusters)
  labels_colors(dendrogram) <- subtype_colors[as.matrix(gbm_subtypes)][order.dendrogram(dendrogram)]
  dendrogram <- dendrogram %>% set("leaves_pch", 19) %>% set("leaves_cex", 1) %>% set("leaves_col", labels_colors(dendrogram))
  dendrogram <- dendrogram %>% set("labels", rep( "", dim(matrix)[2] ))
  
  # reset the plot
  par(mfrow=c(1,1))
  plot(dendrogram)  
  title(main=paste( linkage, "-", distance ))
  
  return( dendrogram )
}
```


```r
cluster_comparison( gbm_data_scaled, c(Mesenchymal="red", Neural="green", Proneural="blue", Classical="orange"), "complete", "euclidean" )
```

![](clustering_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
  
  ```
## 'dendrogram' with 2 branches and 173 members total, at height 53.08395
```


```r
cluster_comparison( gbm_data_scaled, c(Mesenchymal="red", Neural="green", Proneural="blue", Classical="orange"), "single", "manhattan" )
```

![](clustering_files/figure-html/unnamed-chunk-11-1.png)<!-- -->
  
  ```
## 'dendrogram' with 2 branches and 173 members total, at height 701.1403
```


```r
cluster_comparison( gbm_data_scaled, c(Mesenchymal="red", Neural="green", Proneural="blue", Classical="orange"), "average", "euclidean" )
```

![](clustering_files/figure-html/unnamed-chunk-12-1.png)<!-- -->
  
  ```
## 'dendrogram' with 2 branches and 173 members total, at height 44.34999
```

#  Let's god back to our orignal clustering using Euclidean distance and complete linkage (the default values for `dist` and `hclust`). We can extract the four subtrees corresponding to `cutree(dendrogram, 4)` using the list structure of the dendrogram and plot each subtree in a separate plot.


```r
# create 2x2 grid for subplots
par(mfrow=c(2,2))

plot(dendrogram[[1]][[1]])
plot(dendrogram[[1]][[2]])
plot(dendrogram[[2]][[2]])
plot(dendrogram[[2]][[1]])
```

```r
# reset the plot
par(mfrow=c(1,1))
```

#  We can also plot a heatmap showing the expression values that we used to compute the clustering.


```r
heatmap.2(gbm_data_scaled, col=rev(brewer.pal(11,"RdBu")), keysize=0.5, margins=c(1,1), dendrogram="column", Colv=dendrogram, labCol=F, labRow=F, key=F, trace="none")
```

![](clustering_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

The subtype colors are currently associated with the dendrogram but the heatmap also supports plotting of this kind of meta data:


```r
heatmap.2(gbm_data_scaled, col=rev(brewer.pal(11,"RdBu")), keysize=0.5, margins=c(1,1), dendrogram="column", Colv=dendrogram, labCol=F, labRow=F, key=F, trace="none", ColSideColors=subtype_colors[as.matrix(gbm_subtypes)])
```

![](clustering_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

#  It is kind of hard to see where the "cluster" boundaries are, let's break up the heatmap a bit to make that more obvious.


```r
heatmap.2(gbm_data_scaled, col=rev(brewer.pal(11,"RdBu")), keysize=0.5, margins=c(1,1), dendrogram="column", Colv=dendrogram, labCol=F, labRow=F, key=F, trace="none", ColSideColors=subtype_colors[as.matrix(gbm_subtypes)], colsep=cumsum(unlist(lapply(lapply(dendrogram,lapply,unlist),lapply,length))), sepwidth=c(0.5, 0.5))
```

![](clustering_files/figure-html/unnamed-chunk-16-1.png)<!-- -->
  
  Here is a more detailed view of one of the clusters. It is very clearly visible that a small number of patients in this cluster have somewhat different expression patterns.


```r
heatmap.2(gbm_data_scaled, col=rev(brewer.pal(11,"RdBu")), keysize=0.5, margins=c(1,1), dendrogram="column", Colv=dendrogram[[1]][[2]], labCol=F, labRow=F, key=F, trace="none")
```

![](clustering_files/figure-html/unnamed-chunk-17-1.png)<!-- -->
  
  ## k-means Clustering
  
#- k-Means is an extremely popular technique for clustering of high-dimensional data
#- it is fast and generally converges
#- requires number of cluster _k_ a priori
#- it is limited to Euclidean distance
#- can get stuck in local minima (suboptimal solutions)
#- highly dependent on the initialization

### k-means Algorithm

# - *Inputs*: matrix, number of clusters _k_ and maximum number of iterations 
# - *Objective*: find assignment of elements to clusters so that overall distance between cluster centroids and elements is minimal

# - *Initialize*: Select _k_ cluster centroids (see below)
# - *Iterate*:
#  1. Assign elements to cluster based on distance to cluster centroid.
# 1. Compute new location of cluster centroid.
# 1. Repeat.
# - *Terminate*: Once maximum number of iterations is reached or no more change in cluster centroids.



### Initialization of Centroids

- *random*
  - pick _k_ data points from input
- pro: initial centroinds will be near data
- con: they might be close to each other
- *farthest*
  - pick a random data point from input as first centroid
- then select next centroid as the data point farthest away from first centroid
- pro: inital centroids will be distribted
- con: outliers might get selected
- *random + farthest*
  
  ### k-means in Practice
  
  Generate a 2D data set with cluster structure.


```r
random_clusters <- function(points, clusters, seed) {
  set.seed(seed)
  d <- data.frame(x = unlist(lapply(1:clusters, function(i) rnorm(points/clusters, runif(1)*i^2))), 
                  y = unlist(lapply(1:clusters, function(i) rnorm(points/clusters, runif(1)*i^2))))
}

sample_data <- random_clusters( 1000, 5, 9 )
```

Run k-means with default settings.


```r
par(pch=19)
kc <- kmeans(sample_data,5);
plot(sample_data,col=kc$cluster)
```

![](clustering_files/figure-html/unnamed-chunk-19-1.png)<!-- -->
  
#  Let's see what happens if we let the algorithm choose random initializations from the data. 


```r
par(mfcol=c(3,3))
for (i in 1:9) {
kc <- kmeans(sample_data, 5);
plot(sample_data, col=kc$cluster, main=kc$tot.withinss);
}
```

![](clustering_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
par(mfcol=c(1,1))
```

#  In the previous plot we used up to 10 iterations, now we will use just 1. Of course that is a bad idea in general, but this illustrates quite nicely what happens when we end the optimization before reaching an optimum.


```r
par(mfcol=c(3,3))
for (i in 1:9) {
kc <- kmeans(sample_data, 999);
plot(sample_data, col=kc$cluster, main=kc$tot.withinss);
}
```

```
## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration

## Warning: did not converge in 1 iteration
```

![](clustering_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
par(mfcol=c(1,1))
```

#### Determining the Number of Clusters

One way to determine the number of clusters is to plot the total within cluster variance or total within-cluster sum of squares. We will look for the "elbow" in the curve, which is a heuristic to determine the number of clusters. 


```r
tot.withinss <- c();
for (i in 1:20) {
kc <- kmeans(sample_data, i, nstart=100);
tot.withinss <- c( tot.withinss, kc$tot.withinss )
}

plot( 1:20, tot.withinss)
```

#  ![](clustering_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

#  Note that due to the optimization criterion of k-means (distance to centroid), more clusters will always yield a better fit, as illustrated by the curve approaching 0. Therefore, it is advisable to take into account the *model complexity*, which is the number of clusters in this case. There are so-called "information criteria" that can be used for these kind of model selection tasks. One example is the "Bayesian Information Criteriton" (BIC).

#  We can also compute the silhouette width for each element to get another measure of how well the data is clustered. Silhouette width describes how similar an element is to the other members of its cluster relative to how similar it is to the elements in the nearest cluster. Silhouette width ranges from _1_ to _-1_, where _1_ indicates that the element is well clustered, _0_ indicates that the element lies between clusters, and _-1_ indicates that the element should have been assigned to another cluster.


```r
kc <- kmeans(sample_data, 5, nstart=1, iter.max = 1);
si <- silhouette(kc$cluster, dist(sample_data));
par(mfcol=c(1,2));
plot(sample_data,col=kmeans(sample_data,5)$cluster);
plot(si,col=1:5, border=NA);
```

```r
par(mfcol=c(1,1));
```

Now we will apply k-means with `k = 4` to our glioblastoma data from above. We use the Silhouette width to evaluate the clustering.


```r
kc <- kmeans(t(gbm_data_scaled), 4, nstart=100);
si <- silhouette(kc$cluster, dist(t(gbm_data_scaled)));
par(mfcol=c(1,1));
plot(si,col=1:4, border=NA);
```

![](clustering_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
par(mfcol=c(1,1));
```

Compare hierarchical clustering to k-means clustering by plotting cluster assignments along dendrogram.


```r
heatmap.2(gbm_data_scaled, col=rev(brewer.pal(11,"RdBu")), keysize=0.5, margins=c(1,1), dendrogram="column", Colv=dendrogram, labCol=F, labRow=F, key=F, trace="none", ColSideColors=as.character(kc$cluster))
```

![](clustering_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

# How stable are these clusters? What if there is noise in the data or we don't have the same set of patients? We can use bootstrapping to obtain a measure for the stability of the clustering. 


```r
kc.boot <- clusterboot(t(gbm_data_scaled), B=20, bootmethod="boot", clustermethod=kmeansCBI, k=4, seed=1, count=F)
plot(kc.boot)
#### There is a thing called Jacar something (coefficient) - there are sets which are mutually exclusive, but they partition (each cluster is a subset)
# - with the Jacar (or Drakar) coefficent, we look how much they overlap

```

![](clustering_files/figure-html/unnamed-chunk-26-1.png)<!-- -->
  
  The Jaccard similarity coefficient (or "Jaccard index") measures how much "overlap" two sets (here: clusters) have. It is defined as the size of the intersection of two sets divided by the size of their union. Perfect overlap is described by a Jaccard index of 1. The theoretical minimium is 0.

### Fuzzy c-Means Clustering

# Fuzzy clustering is also called soft clustering. The most popular algorithm in the _fuzzy c-means_ algorithm. It works like k-means but has a "fuzziness factor" _m_.


```r
# m = 2 is the default
cc <- cmeans(sample_data, 5, m=2)

par(mfcol=c(2,3))

for( i in 1:5) {
  plot(sample_data,col=rev(heat.colors(100))[ceiling(100 * cc$membership[,i])] )  
  title( paste( "Cluster ", i, " (m = 2)" ))
}
```

![](clustering_files/figure-html/unnamed-chunk-27-1.png)<!-- -->
  
  Less fuzzy:
  
  
  ```r
# m = 2 is the default
cc <- cmeans(sample_data, 5, m=1.5)

par(mfcol=c(2,3))

for( i in 1:5) {
  plot(sample_data,col=rev(heat.colors(100))[ceiling(100 * cc$membership[,i])] )  
  title( paste( "Cluster ", i, " (m = 1.5)" ))
}
```

![](clustering_files/figure-html/unnamed-chunk-28-1.png)<!-- -->
  
  Fuzzier:
  
  
  ```r
# m = 2 is the default
cc <- cmeans(sample_data, 5, m=4)

par(mfcol=c(2,3))

for( i in 1:5) {
  plot(sample_data,col=rev(heat.colors(100))[ceiling(100 * cc$membership[,i])] )  
  title( paste( "Cluster ", i, " (m = 4)" ))
}
```

![](clustering_files/figure-html/unnamed-chunk-29-1.png)<!-- -->
  
  ## Advanced Topics 
  
  
  ### Model-based Clustering
  
  - k-means assumes that all dimensions have the same weight
- that is very limiting because generally clusters are not spheres
- R package `mclust` offers methods for model-based clustering to deal with mixtures of Gaussians
- model-based clustering uses Expectation-Maximization algorithm (EM) to learn models
