library(BalancedSampling)
library(SamplingBigData)

N <- 1000 
n <- 2 
p <- 2 

for( i in 1:10) {
cat(sprintf("%d\n",i))
set.seed(i)
x <- matrix(rnorm(N*p),ncol=p)

prob <- rep(n/N,N)

set.seed(i)
y <-  lpm2(prob,x) 

set.seed(i)
y_kdtree <- SamplingBigData::lpm2_kdtree(prob,x,resample=1)

print( sum(y != y_kdtree) )

}
