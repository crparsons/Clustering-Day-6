#Data Handling
#Problem B
apply(my_data, 1, mean)

#Problem C
head(apply(my_data, 1, min))
head(apply(my_data, 1, max))

head(apply(my_data, 2, min))
head(apply(my_data, 2, max))

#Problem D
which.max(apply(my_data, 1, sd))
max(apply(my_data, 1, sd))

#Problem E
max(apply(my_data, 1, max))
which.max(apply(my_data, 2, max))
which.max(apply(my_data, 1, max))

#Plotting
#A
hist(as.numeric(my_data["EGFR",]), xlab = "EGFR expression", ylab = "Frequency (# of samples)", breaks = 30, main = paste("Histogram of EGFR expression"))

#B
plot(as.numeric(my_data["EGFR", ]), as.numeric(my_data["IDH1", ]), xlab = "EGFR expression", ylab = "IDH1 expression", main = paste("EGFR versus IDH1"), type = "p", pch = 20)

#C
length(grep("ZNF", rownames(my_data)))
grep("ZNF", rownames(my_data))
names <- rownames(my_data)
idx <- grep("ZNF", rownames(my_data))
names[idx]
znfbox <- (apply(my_data[idx, ], 2, mean))
allbox <- (apply(my_data, 2, mean))
boxplot(znfbox, allbox, main = "ZNF gene expression versus all genes", notch = TRUE, ylab = "Expression", names = c("Average ZNF", "Average of All Genes"))

#Functions
#A and B
my_vector <- c(10,5,2,6,8,4,1,9,3,7)

##naming the loop
selectionsort.loop <- function(v) {
  #N is going to be the full length, this is how far the function should work
  N = length(v)
  Position = 1:N
  #go through the whole thing
  for(i in 1:N) {
    min = i 
    #
    for(j in i:N) {
      if(v[j] < v[min]) {
        min = j
      }
      
    }
    
    temp = v[i]
    v[i] = v[min]
    v[min] = temp
    
    temp2 = Position[i]
    Position[i] = Position[min]
    Position[min] = temp2
  }
  #return(v) - this would give you #A
  return(Position)
}

