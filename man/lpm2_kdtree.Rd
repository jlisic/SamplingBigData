\name{lpm2_kdtree}
\alias{lpm2_kdtree}
\title{
Local Pivotal Method 
}
\description{
  The local pivotal method provides a way to perform balanced sampling.  This
  implementation replace linear searches in lpm2, with k-d trees.  K-d trees
  are binary trees used to effectively search high dimensional spaces, and 
  reduce the average computational complexity of lpm2 from O(N^2) to 
  O(N log(N)).  Both nearest neighbor and approximate nearest neighbor 
  searching algorithms are provided.
}
\usage{

  lpm2_kdtree(
     prob,
     x, 
     m=40,
     algorithm = "kdtree", 
     maxCheck = 4,
     termDist = 0.1,
     inOrder = FALSE,
     resample = 1,
     probTree = FALSE,
     returnTree = FALSE,
     returnBounds = FALSE
     ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{
    An array of length N such that the sum of prob is equal to the sample size,
    where the N is the number of rows of \code{x}.
  }
  \item{x}{
    A matrix of N rows and p columns, each row is assumed to be a sampling unit. 
  }
  \item{m}{
    Max leaf size used as a terminal condition for building the k-d tree. When
    probTree is \code{FALSE}, \code{m} is the number of rows of \code{x} held within
    each leaf node.  When probTree is \code{TRUE}, \code{m} is the sum of the
    probability held within each node.  
  }
  \item{algorithm}{
    The algorithm used to search the k-d tree.  The algorithms include "kdtree",
    "kdtree-count", and "kdtree-dist".  The "kdtree" algorithm reproduces the
   lpm2 using a k-d tree for nearest neighbor search.  "kdtree-count" and 
   "kdtree-dist" use approximate nearest neighor searches based on number of
   nodes to check and minimal sufficient distance respectfully.  
  }
  \item{maxCheck}{
    A positive integer scalar parameter only used when the algorithm "kdtree-count"
    is specified.  This parameter is the maximum number of non-empty leaf nodes 
    to check for a nearest neighbor.
  }
  \item{termDist}{
    A positive valued scalar parameter only used when the algorithm "kdtree-dist"
    is specified.  This parameter specifies a minimal sufficient distance to be 
    considered a nearest neighbor.  No tie handling is performed; the first 
    nearest neighbor found will be used. 
  }
  \item{inOrder}{
    A boolean value, \code{TRUE} will return results in order of selection.
    \code{FALSE} will return in order of index.
  }
  \item{resample}{
    The number of samples to return. Resampling builds the k-d tree exactly once
    for all samples.  Each sample will be a distinct column in a matrix.
  }
  \item{probTree}{
    A boolean value, \code{TRUE} will split the k-d tree based on the weighted
    median, using the values in \code{prob}. 
  }
  \item{returnTree}{
    A boolean value, \code{TRUE} will return the node assigment.  This assignment
    will be appended to a returned list.. 
  }
  \item{returnBounds}{
    A boolean value, \code{TRUE} will return the bounds for each node assigment.  
    These bounds will be appended to a returned list. 
  }
}
\value{ 
A vector of selected indexes from the matrix x. If using default values for \code{inOrder}, 
\code{resample}, and \code{algorithm}, the results identical to the lpm2 function when no 
ties exist in the distance function exist.  \code{inOrder=TRUE} will return results in order 
of selection, and \code{resample > 1} will return a matrix with each set of samples returned
as a column vector.  

A list is returned includiing this vector if returnTree or returnBounds is set to TRUE.
}
\references{
Lisic, L. and Cruze, N. (2016). Local Pivotal Methods for Large Surveys. In 
  proceedings, ICES V, Geneva Switzerland 2016. 
}
\author{
  Jonathan Lisic
}
\examples{


N <- 1000
n <- 100
x <- cbind( runif(N), runif(N)) 

set.seed(100)
Cprog <- proc.time()
sampled <- lpm2_kdtree( rep(n/N,N), x)
print("lpm2_kdtree running time")
print(proc.time() - Cprog) 

}
