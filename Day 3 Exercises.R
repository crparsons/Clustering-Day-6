## Chase R. Parsons
## 8/17/2016 
### Day 3 Exercises

##Developing a Smith-Waterman Aligner for two sequences
##Generate some sample DNA data

#DNA example
seq_1 <- c("A","C","A","C","A","C","T","A")
seq_2 <- c("A","G","C","A","C","A","C","A")

seq_1

#Score +1 for match, 0 for no matach
dist_nucleotides <- diag(4)
rownames(dist_nucleotides) <- c("A","C","T","G")
colnames(dist_nucleotides) <- c("A","C","T","G")
dist_identity_nucleotides <- dist_nucleotides
dist_identity_nucleotides

## Create a pairwise matrix using the sample sequences

s <- dist_nucleotides[seq_2,seq_1]
colnames(s) <- c(seq_1)
rownames(s) <- c(seq_2)
s

## Constructing the core aligner
# A. Use nested loops to populate a nucleotide similarity matrix where matches are +2 and non-matches are -1. 
#  Name the rows and columns.
##    A  C  T  G
## A  2 -1 -1 -1
## C -1  2 -1 -1
## T -1 -1  2 -1
## G -1 -1 -1  2

for (j in colnames(dist_identity_nucleotides)[j]) {
  for (i in rownames(dist_identity_nucleotides)[i]) {
   if (i == j) {
     dist_identity_nucleotides[i,j] <- 2
   }
   if (i != j) {
     dist_identity_nucleotides[i,j] <- (-1)
   }
  }
}
dist_identity_nucleotides

# B. Using the distance matrix from (A.), create a pairwise distance matrix based on the sequences to align.
s2 <- dist_identity_nucleotides[seq_2,seq_1]
colnames(s2) <- c(seq_1)
rownames(s2) <- c(seq_2)
s2

# C. Populate an empty matrix of zeroes using the dimensions of each sequence plus one additional cell to pad the upper-left corner, as shown.
paddedm <- matrix(0, nrow = 9, ncol = 9)
rownames(paddedm) <- c("","A","C","A","C","A","C","T","A")
colnames(paddedm) <- c("","A","G","C","A","C","A","C","A")
paddedm

# D. Fill the matrix from (C.) with the maximum similarity scores 
#    by way of nested loops to implement the Smith-Waterman algorithm.

i = 2
j = 2
paddedm[i,j] <- max (
  paddedm[i-1,j]-d,
  paddedm[i,j-1]-d,
  paddedm[i-1,j-1]+2
)
paddedm[i,j]

for (j in 2:ncol(paddedm)) {
  for (i in 2:nrow(paddedm)) {
    paddedm[i,j] <- max(
      0,
      paddedm[i-1,j]-1,
      paddedm[i,j-1]-1,
      paddedm[i-1,j-1]+s2[i-1,j-1]
    )
  }
}
paddedm

# E. What is the optimal alignment score for the two sample sequences? 
#    What position in the maximum similarity-score matrix holds the highest value?
max(paddedm)
max(apply(paddedm,1 , max))
which.max(apply(paddedm, 2, max))
which.max(apply(paddedm, 1, max))

## Recursion
#A. Implement a recursive Smith-Waterman algorithm to align the sample DNA sequences. 
#Fill the maximum similarity scores matrix with a function that calls itself, recursively.
#F is the maximum similarity scores matrix, reset to be filled with all zeroes prior to this function call.
#sw_align is the recursive function, starting in the lower-right of F
F <- matrix(0, nrow = 9, ncol = 9)
rownames(F) <- c("","A","C","A","C","A","C","T","A")
colnames(F) <- c("","A","G","C","A","C","A","C","A")
F

sw_align <- function(i,j){
  if(i==0){
    #return(-d * j)
    return(0)
  }
  else if(j==0){
    #return(-d * i)
    return(0)  
  }
  else{
    F[i,j] <<- max(c(0,sw_align(i-1,j-1)+s[i-1,j-1],sw_align(i-1,j)-d,sw_align(i,j-1)-d))
  }
}

#reset F
F <- H
sw_align(nrow(F),ncol(F))
F

sw_align(nrow(F),ncol(F))
F

# B. Time how long the recursive function takes. Find a way to speed it up and show the new completion time. 
#Hint: have the function print which indices it is currently computing for a clue why it may be taking a long time.
#recursive function, sw_align
time_start = proc.time()
sw_align(nrow(F),ncol(F))
proc.time() - time_start

sw_align(2,2)

time_start = proc.time()
sw_align_fast(nrow(F),ncol(F))
proc.time() - time_start

## Extending the algorithm
# A. Modify the algorithm from (D., above; Constructing the core aligner) to store 
#information about the direction of each cell’s most optimal neighbor(left, up, or upper-left).

# In this example the values correspond to:
# 0 or 1 no direction
# 2 upper-left; diagonal; match
# 3 up; deletion 
# 4 left; insertion
##     A C A C A C T A
##   0 0 0 0 0 0 0 0 0
## A 0 2 4 2 4 2 4 1 2
## G 0 3 2 3 2 3 2 1 3
## C 0 1 2 4 2 4 2 4 4
## A 0 2 3 2 4 2 4 4 2
## C 0 3 2 3 2 4 2 4 4
## A 0 2 3 2 3 2 4 4 2
## C 0 3 2 3 2 3 2 4 4
## A 0 2 3 2 3 2 3 2 2

# B. Print the alignment according the optimal path.
