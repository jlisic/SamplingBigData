/* a quick select median program */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>



void partitionIndex(double ** A, size_t * x1, size_t  * x2, const size_t n, const double p); 


double quantile_quickSelectIndex( double ** A, const size_t k, const size_t n ); 

