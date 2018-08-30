\name{split_sample}
\alias{split_sample}
\title{
Split Sample
}
\description{
  This is a fast implementation of the pivotal method.
}
\usage{

  split_sample(
     prob,
     delta  = exp(-16)
     ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{
    An array of length N such that the sum of prob is equal to the sample size.
  }
  \item{delta}{
    A small real value that is used for tolerance in determining if the value is
    included or excluded from the sample.
  }
}
\value{ 
  An array of indexes from prob.  Indexes with this list are sampled.
}
\references{
Deville, J.-C. and  Tillé, Y. (1998). Unequal probability sampling without replacement through a splitting method. Biometrika 85, 89-101.
}  
\author{
  Jonathan Lisic
}
\examples{


N <- 1000
n <- 100
runif(N)

set.seed(100)
Cprog <- proc.time()
sampled <- split_sample( rep(n/N,N))
print(proc.time() - Cprog) 

}
