/* Copyright (c) 2015-2017 Jonathan Lisic 
 * Last edit: 17/10/19 - 10:57:03 
 * License: GPL (>=2) 
 */  

#ifndef LPM3_HEADER

#define LPM3_HEADER


#include <stdio.h>
#include <stdlib.h>

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"
#include <R_ext/Rdynload.h>
#define PRINTF Rprintf

#include "kdtree_lpm.h"


/********************************/
/* Function Prototypes          */
/********************************/


/* return new prob */
void updateProb( 
    double * xPtr, 
    double * yPtr, 
    double U );


/* function to update the mapping and inverse mapping */
void updateMapping(size_t j,size_t i, size_t * indexMap, size_t * reverseIndexMap);


/* the R interface */
void R_lpm3(
    double * x,
    double * pi,
    int * nPtr,
    int * KPtr, 
    int * mPtr,
    int * algorithmPtr,
    int * maxCountPtr,
    double * termDist,
    int * recordOrder, 
    int * drawsPtr, /* resampling count */
    int * useProbPtr,
    int * nodeAssignment,
    double * bound 
  );



/* split sampling R interface */
void split_sample(
    double * pi,
    size_t  n,
    double * delta,
    size_t * indexMap,
    size_t max_size
); 


// function to find neighbors 
void nn_sample( 
    rootNodePtr r, 
    nodePtr c, 
    double * p,
    double * delta
    ); 


/* split sampling R interface */
void R_split_sample(
    double * pi,
    int * nPtr,
    double * delta
); 


/********** REGISTER FUNCTIONS **********/

static R_NativePrimitiveArgType R_lpm3_t[] = {
      REALSXP, REALSXP, 
      INTSXP, INTSXP, 
      INTSXP, INTSXP, 
      INTSXP, 
      REALSXP,
      INTSXP, INTSXP, 
      INTSXP, INTSXP, 
      REALSXP
};

static R_NativePrimitiveArgType R_split_sample_t[] = {
      REALSXP, INTSXP, REALSXP
};

static const R_CMethodDef cMethods[] = {
     {"R_lpm3", (DL_FUNC) &R_lpm3, 13, R_lpm3_t},
     {"R_split_sample", (DL_FUNC) &R_split_sample, 3, R_split_sample_t},
        {NULL, NULL, 0, NULL}
};

void R_init_myLib(DllInfo *info)
{
     R_registerRoutines(info, cMethods, NULL, NULL, NULL);
     R_useDynamicSymbols(info, TRUE); 
}






#endif
