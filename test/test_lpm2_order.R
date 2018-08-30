library(BalancedSampling)
library(SamplingBigData)
N <- 1000
n <- 100
p <- 10

for( i in 1:2) {
i <- 1

set.seed(i)
x <- matrix(rnorm(N*p),ncol=p)

prob <- rep(n/N,N)

set.seed(i)
y <- lpm2(prob,x)

set.seed(i)
y_kdtree_order <- SamplingBigData::lpm2_kdtree(prob,x,inOrder=TRUE)

print( sum(y != sort(y_kdtree_order) ))

}

