
#include "median.h"

#ifdef CLI
 #define PRINTF printf
 #define RUNIF (double) rand() / RAND_MAX
#else
 #include "R.h"
 #include "Rmath.h"
 #define PRINTF Rprintf
 #define RUNIF runif(0.0,1.0)
#endif




/* partition index */
void partitionIndex(double ** A, size_t * x1, size_t  * x2, const size_t n, const double p) {
  
  size_t i;
  
  size_t m1;

  *x1 = 0;
  *x2 = n;
  
  double ** X =calloc( n, sizeof(double *));

  for( i = 0;i < n; i++) {

    if( *A[i] < p) {
      X[*x1] = A[i]; 
      *x1 = *x1 + 1; 
    } else if( *A[i] > p) {
      *x2 = *x2 - 1;
      X[*x2] = A[i]; 
    } 
  }

  m1 = *x1;

  for( i = 0;i < n; i++) { 
    if( *A[i] == p) {
      X[m1] = A[i];
      m1++;
    } 
  }

  for( i = 0; i < n; i++) A[i] = X[i];

  free(X);

  return;
}




//Quickselect Index
double quantile_quickSelectIndex( double ** A, const size_t k, const size_t n ) {

  double p;
  size_t a1, a2;

  //1. find pivot element p (median)
//  p = (*A[0] + *A[n-1] + *A[n/2])/3.0;
  p = *A[n/2];

  //2. Partition A by p; let A1 < p, A2 = p, A3 > p
  partitionIndex( A, &a1, &a2, n, p); 

  //3. If k in A1 
  if(k < a1) 
    return quantile_quickSelectIndex(A, k, a1);
  //4. Else if k in A3 
  else if(k >= a2 ) 
    return quantile_quickSelectIndex(A + a2, k - a2, n - a2);

  //5. if k in A2 return p 
  return(p);
}


