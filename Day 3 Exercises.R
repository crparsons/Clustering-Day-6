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
# A. Use nested loops to populate a nucleotide similarity matrix where matches are +2 and non-matches are -1. Name the rows and columns.


# B. Using the distance matrix from (A.), create a pairwise distance matrix based on the sequences to align.


# C. Populate an empty matrix of zeroes using the dimensions of each sequence plus one additional cell to pad the upper-left corner, as shown.

# D. Fill the matrix from (C.) with the maximum similarity scores by way of nested loops to implement the Smith-Waterman algorithm.
 
# E. What is the optimal alignment score for the two sample sequences? What position in the maximum similarity-score matrix holds the highest value?



## Recursion
#A. Implement a recursive Smith-Waterman algorithm to align the sample DNA sequences. Fill the maximum similarity scores matrix with a function that calls itself, recursively.
#F is the maximum similarity scores matrix, reset to be filled with all zeroes prior to this function call.
#sw_align is the recursive function, starting in the lower-right of F

sw_align(nrow(F),ncol(F))
F

# B. Time how long the recursive function takes. Find a way to speed it up and show the new completion time. Hint: have the function print which indices it is currently computing for a clue why it may be taking a long time.
#recursive function, sw_align
time_start = proc.time()
sw_align(nrow(F),ncol(F))
proc.time() - time_start

sw_align(2,2)

time_start = proc.time()
sw_align_fast(nrow(F),ncol(F))
proc.time() - time_start

## Extending the algorithm
# A. Modify the algorithm from (D., above; Constructing the core aligner) to store information about the direction of each cell’s most optimal neighbor(left, up, or upper-left).

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
