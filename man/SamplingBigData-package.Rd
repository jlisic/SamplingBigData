\name{SamplingBigData-package}
\alias{SamplingBigData-package}
\alias{SamplingBigData}
\docType{package}
\title{
Sampling Methods for Big Data
}
\description{
Select sampling methods for probability samples using large data sets.  This 
includes spatially balanced sampling in multi-dimensional spaces with any 
prescribed inclusion probabilities. All implementations are written in C using 
efficient data structures such as k-d trees that easily scale to several million
rows on a modern desktop computer. 
}
\author{
Jonathan Lisic, Anton Grafström 

Maintainer: Jonathan Lisic <jlisic@gmail.com>

Webpage: https://github.com/jlisic/SamplingBigData
}
\references{
Deville, J.-C. and  Tillé, Y. (1998). Unequal probability sampling without 
  replacement through a splitting method. Biometrika 85, 89-101.

Grafström, A. (2012). Spatially correlated Poisson sampling. Journal of 
  Statistical Planning and Inference, 142(1), 139-147.

Grafström, A. and Lundström, N.L.P. (2013). Why well spread probability 
  samples are balanced. Open Journal of Statistics, 3(1).

Grafström, A. and Schelin, L. (2014). How to select representative samples.
  Scandinavian Journal of Statistics.

Grafström, A., Lundström, N.L.P. and Schelin, L. (2012). Spatially balanced
  sampling through the Pivotal method. Biometrics 68(2), 514-520.

Grafström, A. and Tillé, Y. (2013). Doubly balanced spatial sampling with 
  spreading and restitution of auxiliary totals. Environmetrics, 24(2), 120-131.

Lisic, L. and Cruze, N. (2016). Local Pivotal Methods for Large Surveys. In 
  proceedings, ICES V, Geneva Switzerland 2016. 

}
\keyword{ sampling }
\examples{
# *********************************************************
# check inclusion probabilities
# *********************************************************
set.seed(1234567);
p = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
N = length(p);
X = cbind(runif(N),runif(N));
p1 = p2 = p3 = p4 = rep(0,N);
nrs = 1000; # increase for more precision 
for(i in seq(nrs) ){

  # lpm2 kdtree
  s = lpm2_kdtree(p,X);
  p1[s]=p1[s]+1;
  
  # pivotal method 
  s = split_sample(p);
  p2[s]=p2[s]+1;
  
}
print(p);
print(p1/nrs);
print(p2/nrs);
}
