# Linear Algebra and Dimensionality Reduction
# Peter Kharchenko  
# August 18, 2016  




## Working with matrices


```r
# constructing a matrix
a <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
a <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow=T)

print(a)

# multiplying each rows and columns
a*2
```

```
##      [,1] [,2]
## [1,]    2    4
## [2,]    6    8
## [3,]   10   12
```

```r
# multiply each row by a different number
a * c(1,2,3)
```

```
##      [,1] [,2]
## [1,]    1    2
## [2,]    6    8
## [3,]   15   18
```

```r
# multiply each column by a different number
t( t(a) * c(1,2) )
```

```
##      [,1] [,2]
## [1,]    1    4
## [2,]    3    8
## [3,]    5   12
```

```r
# multiply matrix by a vector
a %*% c(1,2)
```

```
##      [,1]
## [1,]    5
## [2,]   11
## [3,]   17
```

```r
# multiply two matrices
b <- matrix(c(4,2,-1,0), nrow=2, ncol=2)
a %*% b
print(b)
```

```
##      [,1] [,2]
## [1,]    8   -1
## [2,]   20   -3
## [3,]   32   -5
```

```r
t(b) %*% t(a)
```

```
##      [,1] [,2] [,3]
## [1,]    8   20   32
## [2,]   -1   -3   -5
```

```r
# dot product
c(1,7,2,5) %*% c(1,2,3,4)
```

```
##      [,1]
## [1,]   41
```

```r
# outer product
c(1,7,2,5) %*% t(c(1,2,3,4))
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    1    2    3    4
## [2,]    7   14   21   28
## [3,]    2    4    6    8
## [4,]    5   10   15   20
```


## Solving systems of linear equations
Here we try to solve a simple linear system $Ax=b$ for x

```r
a <- matrix(c(1,0.4,-2,1), nrow=2)
b <- c(-2,1)

x <- solve(a,b)
x
```

```
## [1] 0 1
```

```r
a %*% x
```

```
##      [,1]
## [1,]   -2
## [2,]    1
```

```r
# through matrix inverse
solve(a) %*% b
```

```
##      [,1]
## [1,]    0
## [2,]    1
```

## Matrix eigen decomposition

```r
a <- matrix(c(2,2,3,1),nrow=2)
e <- eigen(a)
e
```

```
## $values
## [1]  4 -1
## 
## $vectors
##           [,1]       [,2]
## [1,] 0.8320503 -0.7071068
## [2,] 0.5547002  0.7071068
```

```r
(a %*% e$vectors[,1]) / e$values[1]
```

```
##           [,1]
## [1,] 0.8320503
## [2,] 0.5547002
```

## Singular Value Decomposition
R has basic functions for all kinds of matrix decompositions. Here we will examine SVD.

```r
s <- svd(a)
s
```

```
## $d
## [1] 4.1306486 0.9683709
## 
## $u
##            [,1]       [,2]
## [1,] -0.8649101 -0.5019268
## [2,] -0.5019268  0.8649101
## 
## $v
##            [,1]       [,2]
## [1,] -0.6618026  0.7496782
## [2,] -0.7496782 -0.6618026
```

```r
# both U and V are unitary matrices:
s$u %*% t(s$u)
```

```
##               [,1]          [,2]
## [1,]  1.000000e+00 -1.110223e-16
## [2,] -1.110223e-16  1.000000e+00
```

```r
s$v %*% t(s$v)
```

```
##      [,1] [,2]
## [1,]    1    0
## [2,]    0    1
```

```r
# finally, the product gives us back the matrix a
s$u %*% diag(s$d) %*% s$v 
```

```
##           [,1]     [,2]
## [1,] 2.7287642 -2.35666
## [2,] 0.7442084 -2.10859
```

```r
s$u %*% diag(s$d)
```

```
##           [,1]       [,2]
## [1,] -3.572640 -0.4860513
## [2,] -2.073283  0.8375538
```

## PCA
#To illustrate principal component analysis, 
#let's first simulate some data from the multivariate normal distribution.

```r
install.packages('MASS')
library(MASS)
covm <- matrix(c(10,3,3,2),2,2)
covm
```

```
##      [,1] [,2]
## [1,]   10    3
## [2,]    3    2
```

```r
set.seed(0)
x <- mvrnorm(1000,rep(0,2),covm)
plot(x)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

The two variables are well correlated, so the default axes are not very effective. We could rotate the data so that the first axis explains more variation, and the second axis is orthogonal to it. This is what PCA will try to do. 


```r
# let's re-estimate our covariance matrix ... normally we would need to call scale() to center/scale the matrix, but it's already centered in our case
t(x) %*% x / (nrow(x) -1)
```

```
##          [,1]     [,2]
## [1,] 9.995077 2.933384
## [2,] 2.933384 2.033096
```

```r
cov(x)
```

```
##          [,1]     [,2]
## [1,] 9.993314 2.931699
## [2,] 2.931699 2.031485
```

```r
# let's perform an eigenvalue decomposition
e <- eigen(cov(x))
e$vectors
```

```
##            [,1]       [,2]
## [1,] -0.9500555  0.3120809
## [2,] -0.3120809 -0.9500555
```

```r
# are they really eigenvectors? 
cov(x) %*% e$vectors[,1] / e$values[1]
```

```
##            [,1]
## [1,] -0.9500555
## [2,] -0.3120809
```

```r
# check for all the vectors
cov(x) %*% e$vectors - e$vectors %*% diag(e$values)
```

```
##              [,1]         [,2]
## [1,] 1.776357e-15 5.551115e-17
## [2,] 0.000000e+00 0.000000e+00
```

```r
# let's see where they point
plot(x)
segments(0,0,e$vectors[1,1],e$vectors[2,1],col='red',lwd=2)
segments(0,0,e$vectors[1,2],e$vectors[2,2],col='green',lwd=2)
arrows(0, 0, e$vectors[1, ], e$vectors[2, ], length = 0.1, angle = 15, code = 2,col= c("red","green"),lwd=3)
legend(x='topleft',lty=c(1,1),lwd=c(2,2),col=c('red','green'),legend=c('PC1','PC2'),bty='n')
```

![](LinAlgebra_files/figure-html/unnamed-chunk-6-1.png)<!-- -->
  
  ```r
# similarly, we perform PCA using SVD of the matrix x itself
s <- svd(x)
str(s)
```

```
## List of 3
##  $ d: num [1:2] 104.6 32.7
##  $ u: num [1:1000, 1:2] -0.04 0.0104 -0.0422 -0.0404 -0.0132 ...
##  $ v: num [1:2, 1:2] 0.95 0.312 -0.312 0.95
```

```r
s$d ^ 2 / (nrow(x) -1)
```

```
## [1] 10.959087  1.069086
```

```r
e$values
```

```
## [1] 10.95634  1.06846
```

```r
s$v
```

```
##           [,1]       [,2]
## [1,] 0.9500141 -0.3122070
## [2,] 0.3122070  0.9500141
```

```r
# or we can use a built-in R PCA routine
p <- princomp(x)
str(p)
```

```
## List of 7
##  $ sdev    : Named num [1:2] 3.31 1.03
##   ..- attr(*, "names")= chr [1:2] "Comp.1" "Comp.2"
##  $ loadings: loadings [1:2, 1:2] -0.95 -0.312 0.312 -0.95
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:2] "Comp.1" "Comp.2"
##  $ center  : num [1:2] 0.042 0.0401
##  $ scale   : num [1:2] 1 1
##  $ n.obs   : int 1000
##  $ scores  : num [1:1000, 1:2] 4.24 -1.04 4.46 4.28 1.43 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:2] "Comp.1" "Comp.2"
##  $ call    : language princomp(x = x)
##  - attr(*, "class")= chr "princomp"
```

```r
matrix(p$loadings,nrow=2)
```

```
##            [,1]       [,2]
## [1,] -0.9500555  0.3120809
## [2,] -0.3120809 -0.9500555
```

```r
# we can project the points onto a principal component
str(x %*% e$vectors[,1])
```
?str
```
##  num [1:1000, 1] 4.19 -1.09 4.41 4.23 1.38 ...
```

```r
# plot points in the PC1/ PC2 space
plot(x %*% e$vectors[,1] , x %*% e$vectors[,2], xlab='PC1', ylab='PC2',col=1:nrow(x),pch=19,cex=0.8)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-6-2.png)<!-- -->
  
  ```r
# amount of standard deviation explained
p$sdev
```

```
##   Comp.1   Comp.2 
## 3.308381 1.033146
```

```r
sqrt(var(x %*% e$vectors[,1]))
```

```
##          [,1]
## [1,] 3.310036
```

```r
# fraction of variance explained by each component
var(x %*% e$vectors[,1]) / sum(apply(x,2,var))
```

```
##           [,1]
## [1,] 0.9111453
```

```r
var(x %*% e$vectors[,2]) / sum(apply(x,2,var))
```

```
##            [,1]
## [1,] 0.08885471
```

```r
# inline ..
apply(e$vectors,2,function(v) { var(x %*% v)/sum(apply(x,2,var))})
```

```
## [1] 0.91114529 0.08885471
```

# A key advantage fo the PCA is then the ability to reduce the dimensionality by ignoring separation captured by the eigenvectors with small eigenvalues. For intance, taking first few principal components. Let's add more non-informative dimensions to the data and see if we can separate out the original, informative dimensions.

```r
str(x)
```

```
##  num [1:1000, 1:2] -4.06 1.61 -4.23 -4.44 -1.77 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : NULL
##   ..$ : NULL
```

```r
# add 10 more dimensions
set.seed(0)
x2 <- cbind(x, do.call( cbind, lapply(1:10,function(i) rnorm(nrow(x)) ) ))
str(x2)
```

```
##  num [1:1000, 1:12] -4.06 1.61 -4.23 -4.44 -1.77 ...
```

```r
p <- princomp(x2)
barplot(p$sdev,las=2)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
p$loadings[,1:3]
```

```
##             Comp.1       Comp.2       Comp.3
##  [1,]  0.909454042  0.216504834  0.015012565
##  [2,]  0.298351448 -0.667531919 -0.043971624
##  [3,] -0.288585816  0.001717893 -0.000101646
##  [4,]  0.004553584  0.701741222  0.046462535
##  [5,] -0.003459885  0.021064234  0.163837913
##  [6,]  0.014663415 -0.024302666 -0.250025292
##  [7,]  0.001439987 -0.004864176 -0.563782485
##  [8,]  0.007181365  0.073576939 -0.415649695
##  [9,]  0.011824024  0.006775705  0.358204802
## [10,]  0.003977541  0.022510769 -0.128969991
## [11,] -0.010994147  0.080337430 -0.436716635
## [12,]  0.004485492  0.040060010 -0.282924223
```

```r
e$vectors
```

```
##            [,1]       [,2]
## [1,] -0.9500555  0.3120809
## [2,] -0.3120809 -0.9500555
```

```r
plot(x2 %*% p$loadings[,c(1,2)], xlab='PC1', ylab='PC2',col=1:nrow(x),pch=19,cex=0.8)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

## ICA
#Independent component analysis tries to identify a different set of basis (not necessarily orthonormal). The independence is achieved by 1. minimizing the mutual information (i.e. KL dirergence,), and 2. maximizing "non-Gaussianity" of the projections.

```r
install.packages('fastICA')
library(fastICA)
i <- fastICA(x,2)

plot(x)
arrows(0, 0, e$vectors[1, ], e$vectors[2, ], length = 0.1, angle = 15, code = 2,col= c("red","darkred"),lwd=3)
arrows(0, 0, i$A[,1], i$A[,2], length = 0.1, angle = 15, code = 2,col= c("aquamarine1","aquamarine4"),lwd=3)
legend(x='topleft',lty=c(1,1),lwd=c(3,3),col=c('red','aquamarine1'),legend=c('PCA','ICA'),bty='n')
```

![](LinAlgebra_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
# can try on an actual expression data
ed <- read.table("http://pklab.med.harvard.edu/peterk/PRAD_normdata.txt",header=T,nrows=1e3)
str(rownames(ed))
```

```
##  chr [1:1000] "TSPAN6" "TNMD" "DPM1" "SCYL3" "C1orf112" ...
```

```r
ed <- scale(t(ed),center=T,scale=T)
p.ed <- prcomp(ed)

plot(p.ed$sdev)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
str(sort(abs(p.ed$rotation[,1]),decreasing=T))
```

```
##  Named num [1:1000] 0.0647 0.062 0.0611 0.0609 0.0605 ...
##  - attr(*, "names")= chr [1:1000] "CEP68" "METTL24" "ICA1" "IFFO1" ...
```

```r
plot(ed %*% p.ed$rotation[,c(1,2)], xlab='PC1', ylab='PC2',col=1:nrow(x),pch=19,cex=0.8)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

## Other dimensionality reduction methods
R packages implement a variety of other dimensionality techniques. 
For instance, Multidimensional Scaling (MDS) will minimize the discrepancy between high-dimensional and low-dimensional distances between points.


```r
# calcualte Eucledian distance between samples (btw. not a great option for expression magnitudes)
ed.d <- dist(ed)
str(as.matrix(ed.d))
```

```
##  num [1:86, 1:86] 0 24.7 56.1 54.6 43.4 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:86] "TCGA.HC.8260.11A" "TCGA.HC.8259.11A" "TCGA.EJ.7123.11A" "TCGA.G9.6496.01A" ...
##   ..$ : chr [1:86] "TCGA.HC.8260.11A" "TCGA.HC.8259.11A" "TCGA.EJ.7123.11A" "TCGA.G9.6496.01A" ...
```

```r
# run MDS
m.ed <- cmdscale(ed.d,k=2)
plot(m.ed, xlab='MDS dim 1', ylab='MDS dim 2',col=1:nrow(x),pch=19,cex=0.8)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

tSNE also tries to match up high- and low-dimensional distances, but it does so by minimizing the differences between the 'neighbor probabilitiy' distributions in high- and low- dimensional space, modeling the former using normal distribution, and the later using Student t distribution. This gives more stability wrt. to capturing the exact neighbor dependencies.

## tSNE

```r
install.packages('Rtsne')
library(Rtsne)
```

```
## Warning: package 'Rtsne' was built under R version 3.2.5
```

```r
t.ed <- Rtsne(ed.d,is_distance=TRUE,perplexity=15)
plot(t.ed$Y, xlab='dim 1', ylab='dim 2',col=1:nrow(x),pch=19,cex=0.8)
```

![](LinAlgebra_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
# hard to match up, let's use the colors to show PC1 score of each sample
nrow(ed)
```

```
## [1] 86
```

```r
cols <- colorRampPalette(c('blue','red'))(nrow(ed))[rank( ed %*% p.ed$rotation[,c(1)] )]
plot(t.ed$Y, xlab='tSNE dim 1', ylab='tSNE dim 2',pch=19,cex=0.8,col=cols)
¢```

![](LinAlgebra_files/figure-html/unnamed-chunk-10-2.png)<!-- -->
  
  
  