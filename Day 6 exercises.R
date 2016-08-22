#  Day 6 exercises
#  Aug 22, 2016
#  Chase Parsons 

rm(list=ls()) # clean the workspace
setwd("/Users/Chase/Documents/Biomedical Informatics/Boot Camp") #replace location with local path_to_file
prad_data <- read.table("PRAD_norm_counts.txt",sep="\t",header=T)
res_deseq2 <- read.table("RES_DESeq2_PRAD.txt",sep="\t",header=T)
gbm_data <- read.table("GBM_norm_counts_filt.txt",sep="\t",header=T)


### Data preparation ###

prad_data <-log2(as.matrix(prad_data)+1) # log2-transformed matrix with offset =1

# **A.** Prior to clustering data, it is preferable to rescale variables (i.e. genes) for comparability. NOTE: the function "scale" scales the columns of the matrix, therefore the matrix has to be transposed using t() command.

head(prad_data_scaled[,1:5])

##          TCGA_HC_8260_11A TCGA_HC_8259_11A TCGA_EJ_7123_11A TCGA_G9_6496_01A TCGA_EJ_7781_01A
## TSPAN6        -0.23407687      -0.61373633       -2.7864599        0.6515002      -0.30269024
## TNMD          -0.07024842       0.44020751        0.2243651       -0.1078316      -0.55956643
## DPM1           0.27614379       0.98379431       -0.9089804       -1.6804555      -0.65729243
## SCYL3         -0.47652270       0.43371233       -1.2617604       -1.7531606       0.86863311
## C1orf112       0.07150625       0.09609684       -0.8374033       -1.9955209      -0.04052485
## FGR           -0.42508305       0.35207128        1.1097773        0.7020206       0.07906976

### Partitioning clustering: K-means  ###

# K-means clustering is the most popular partitioning method. It requires the user to specify the number of clusters. First look at help(kmeans) to become familiar with the function.

#  **A.** Try to perform k-means clustering for all the samples. Set the number of clusters to k=2. How the samples are distributed across the two clusters with respect to the classification "Primary Tumor" vs "Solid Tissue Normal"?
#  NOTE: k-means uses a random starting solution and thus different runs can lead to different results. Therefore, trying several random starts is often recommended. Set the parameter "nstart" equal to 100.

#  samp_info<-read.table("PRAD_sample_info.txt",sep="\t",header=T)

##   condition
##    Primary Tumor Solid Tissue Normal
##  1            42                  10
##  2             1                  33

#  **B.** Try to extract the sub-matrix keeping the top 100 genes with the lowest p-value and repeat the clustering (NOTE: use the command "order"). What happens to the clustering results with respect to the classification "Primary Tumor" vs "Solid Tissue Normal"? (NOTE: need to re-scale the sub-matrix)

##   condition
##    Primary Tumor Solid Tissue Normal
##  1            43                   3
##  2             0                  40

#  **C.** Try to repeat the clustering by genes using k=2 and plot the results using the command clusplot of the package "cluster". NOTE: use "cex=0.7".

install.packages("cluster")
library(cluster)

#  The two clusters seems well separated. How much well? We can use the silhouette method to evaluate the clustering.

#  For each gene i, let a(i) be the average dissimilarity of gene i (e.g. euclidean distance) with all the other genes within the same cluster. We can interpret a(i) as how well the gene i is assigned to its cluster (the smaller the value, the better the assignment). We then define the average dissimilarity of the gene i to a cluster c as the average of the distance from the gene i to all the genes belonging to c.

#  Let b(i) be the lowest average dissimilarity of gene i to any other cluster, of which gene i is not a member. The cluster with this lowest average dissimilarity is said to be the "neighbouring cluster" of gene i because it is the next best fit cluster for the gene i. We now define a silhouette:

#  As a(i) is a measure of how dissimilar the gene i is to its own cluster, a small value means it is well matched. Furthermore, a large  b(i) implies that the gene i is badly matched to its neighbouring cluster. Thus an s(i) close to one means that the gene is appropriately clustered.

#  **D.** Calculate the silhouette for the clustering results by gene using k=2. Plot the silhouette and repeat the same procedure using k=3 and look at the differences.

install.packages("fpc")
library(fpc)

cl.genes3 <- kmeans(prad_data_scaled[ind,], centers=3, nstart=100)
sil3 <- silhouette(cl.genes3$cluster, dd)
plot(sil3, main="Silhouette plot k=3")

### Hierarchical clustering  ###
#  **A.** Let's repeat the clustering analysis by genes and samples using the function hclust, with euclidean distance and average linkage, using k=2. Use functions "hclust" and "cutree". Compare the results with the k-means.


hcl.samples
1  2
1  0 46
2 40  0
```

```
hcl.genes
1  2
1 43  0
2  0 57
```

#  **B.** Let's visualize these results by a heatmap plot. The pheatmap package provides very nice figures which can be easily scaled according to the dimensions of our data. Install and load the pheatmap package and use the pheatmap function to realize the heatmap plot.
Use "fontsize"=4, "clustering_method"="average", "cutree_rows" = 2, "cutree_cols" = 2

install.packages("pheatmap")

#  Now, let's consider the GBM data. In this case all the samples are primary tumors, but we can indentify sub-groups of samples that might correspond to different cancer sub-types. A plot of the within groups sum of squares by number of clusters extracted can help to determine the appropriate number of clusters. The function cluster.stats() in the "fpc" package provides a series of several measures and criteria for comparing different clustering solutions.

#  **C.** Perform the hierarchical clustering per sample, using euclidean distance and single, average and complete linkage. Use the function cluster.stats and plot the within sum of squares (within.cluster.ss) and the average silhouette (avg.silwidth) as shown in the pictures, for k=2:20. Use funtions plot, lines and legend.


wss_s<-wss_c<-wss_av<-(nrow(gbm_data_scaled)-1)*sum(apply(gbm_data_scaled,2,var)) #wss for k=1
for(i in 2:20)
{hclres.samples <- hclust(d.samples,method="single")
hcl.samples <- cutree(hclres.samples, k=i)
wss_s[i] <- cluster.stats(d.samples, hcl.samples)$within.cluster.ss
...
}
...
```

![](Images_Day6/WSS_HC.png)

```r
sil_s<-sil_c<-sil_av<-NA
for(i in 2:20)
{hclres.samples <- hclust(d.samples,method="single")
hcl.samples <- cutree(hclres.samples, k=i)
sil_s[i] <- cluster.stats(d.samples, hcl.samples)$avg.silwidth
...
}
...
```

![](Images_Day6/SIL_HC.png)

The pvclust( ) function in the pvclust package provides p-values for hierarchical clustering based on multiscale bootstrap resampling. Clusters that are highly supported by the data will have large p values. Interpretation details are provided here: http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust/. Be aware that pvclust clusters by columns, not rows. 

**C.** Install and download pvclust. Use functions pvclust and pvrect to perform the analysis, using average linkage and euclidean distance.

```r
install.packages("pvclust")
library(pvclust)
...
``` 

![](Images_Day6/PVclust.png)

pvclust provides two types of p-values: AU (Approximately Unbiased) p-value and BP (Bootstrap Probability) value. AU p-value, which is computed by multiscale bootstrap resampling, is a better approximation to unbiased p-value than BP value computed by normal bootstrap resampling.
pvclust performs hierarchical cluster analysis via function hclust and automatically computes p-values for all clusters contained in the clustering of original data. It also provides graphical tools such as plot function or useful pvrect function which highlights clusters with relatively high/low p-values.
For a cluster with AU p-value > 0.95, the hypothesis that "the cluster does not exist" is rejected with significance level 0.05.