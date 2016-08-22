# Day 4 R 
# Chase R Parsons
# August 18

biocLite("curl") #run if necessary
library(curl)
my_data <- read.table(curl("https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/unifiedScaledFiltered.txt"),sep="\t",header=1)

# A. Before we get started, we should examine the range of our gene expression vectors. Which genes have the highest and lowest mean expression? What about the variances of the genes? What does this tell you about the data?
max(apply(my_data,1 , max))
which.max(apply(my_data, 2, max))
which.max(apply(my_data, 1, max))

min(apply(my_data,1 , min))
which.min(apply(my_data, 2, min))
which.min(apply(my_data, 1, min))

var(lapply(my_data, var))

# B. Perform a PCA on the GMB expression data, paying special attention to centering and scaling data (hint: what are the default options for the R function “prcomp”?) How many PCs do you find?
p <- princomp(my_data)
str(p)

## I THINK I find 202 different PCs

# C. What proportion of total variance is explained by the first PC? What proportion of total variance is explained by the first 10 PCs?
p$vectors / 
summary(p)
***

# D. Plot the first fifteen principal components by variance. How many should be retained?
p11 <- princomp(my_data, npcs = 15)
plot(p11, npcs = 15, type = "lines")

# E. What genes have the highest absolute loading values for the first principal component? For the second component? Would you expect to find the same genes with the highest loadings for the first and second PCs?



# F. Create a scatterplot of the values for the first and second PCs. Label the points with the sample names, and note which samples stand out. Does this change if you plot a different set of components? (Hint: Use plot() to create a scatterplot and text() to add labels to those points)

