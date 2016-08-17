getwd()

#simple console test
3+5

y <- 3 
x <- 5 

number <- x + y

#creating numeric vector
glengths <- c(4.6, 3000, 50000)

class(glengths)

#creating character vector
species <- c("ecoli", "human", "corn")
class(species)

# factors vectors 
expression <- c("low", "high", "medium", "high", "low", "medium", "high")

expression <- factor(expression)

# data frame
df <- data.frame(species, glengths)

# creating a list

list1 <- list(species, df, number)

# adding 90 to the end of the glengths vector (or 30 to the beginning)

glengths <- c(glengths, 90) # adding at the end 

glengths <- c(30, glengths) # adding at the beginning

sqrt(81)
sqrt(glengths)

round(3.14159, digits = 2)

round(3.145)

# indexes
age <- c(15, 22, 45, 52, 73, 81)

age[-5]

metadata[1,5]

# subset



rpkm_data <- read.csv("counts.rpkm.csv")

head(rpkm_data)

ncol(rpkm_data)

A <- c(2,3,5,8,9,11)   # odd numbers
B <- c(2,4,6,8,10,12)  # even numbers

# test to see if any of A are in B  
any(A %in% B)
all(A %in% B)

A <- c(1,2,3,4,5)
B <- c(5,4,3,2,1)  # same numbers but backwards 

# test to see if any of A are in B
A %in% B

# test to see if any of A is equal to B
A == B

# use all to check if they are a perfect match
all(A == B)