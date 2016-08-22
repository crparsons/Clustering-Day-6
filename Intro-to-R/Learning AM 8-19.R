## Learning in AM

runif(10)

n <- str(rbind(runif(10), runif(10), runif(10)))

print(n)

n.genes <- 5
n.samples <- 10

x <- lapply(1:n.genes, function(i) runif(n.samples))

str(do.call(rbind, x))

str(x)

# or call inline
str(do.call(rbind, lapply(1:n.genes, function(i) runif(n.samples)))
    
1:nrow(m) %% 2 ==0
ifelse(1: nrow(m) %%2 == 0, "even", "odd"


rpois(n = 10, lambda = 500)
#remainder operator:
1024 %% 1000

# simulation
n.rounds <- 10000
x <- lapply(1:n.rounds, function(i) {sum(rpois(n = 1e4, lambda = 5000)) %% 1000})
str(unlist(x))}

library(parallel)

system.time( x <- mclapply(1:n.rounds, function(i) {sum(rpois(n=1e5, lambda = 500)) %% 1000}, mc.cores = 5))

ls()

Browse[1] > l


# easy way using built-in functions
k <- 10; a <- 3
sum(a + 0:k)

# sanity check
(k+1)*(a + a+k)/2


#But let's try to design a space-efficient function $f(k,a)$ to calcualte the result

f <- function(k,a=3) {
  if(k==0) { 
  return(a) 
  } else {
    a+k + f(k-1,a=a)
  }
}
f

sum(3,5)

# f(k)
There's also a built-in $Reduce()$ function in R to perform recursive on a binary function:
  
# Reduce(f , c(x1 ,x2 , x3, x4) ) == f( f( f(x1,x2), x3), x4)
  Reduce(sum,a+0:k)