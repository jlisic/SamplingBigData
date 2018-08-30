/* Copyright (c) 2016  Jonathan Lisic 
 * Last edit: 
 * License: GPL V3 
 */  


#include "kdtree_lpm.h"
#include "median.h"

/* printf fixing */
#ifdef CLI 
 #define PRINTF printf
 #define RUNIF (double) rand() / RAND_MAX  
#else 
 #include "R.h"
 #include "Rmath.h"
 #define PRINTF Rprintf
 #define RUNIF runif(0.0,1.0)
#endif

// Tree
// The tree has a few key components:
// root->pointerIndex: each element points to the i'th element in the tree 


// function to create a new Tree  
rootNodePtr createTree( size_t K, size_t leafSize, size_t n, double * data ) {
  
  rootNodePtr y = malloc( sizeof( rootNode ) );
  
  y->pointerIndex = calloc( n, sizeof( size_t * ) );
  y->nodeIndex = calloc( n, sizeof( size_t * ) );

  // setup root node 
  y->K = K;
  y->leafSize = leafSize;
  y->root = NULL;
  y->data = data;
  y->n = n;

  return(y);
}




// function to create a new Tree 
void deleteTree( rootNodePtr r ) {
 
  if(r->pointerIndex != NULL) free( r->pointerIndex ); 
  if(r->nodeIndex != NULL) free( r->nodeIndex ); 
  r->pointerIndex = NULL;
  r->data = NULL;
    
  deleteNode( r, r->root ); 

  free(r);
  return;
}




// add an element to the tree 
nodePtr buildIndex( 
    rootNodePtr r,      // root pointer 
    size_t dim,         // current dim
    size_t m,           // current length of obs
    size_t * indexPtr,  // pointer to obs indexes 
    int useProb,        // determine if we use probability to build an index
    double * prob,
    size_t * nodeIdentity
  ) {
 
  size_t i,K; 
  size_t * indexLeftPtr = NULL;
  size_t * indexRightPtr = NULL;
  size_t indexLeftSize;
  size_t indexRightSize;
  double probSum = 0;

  nodePtr c = createNode(r);

  // record to the tree structure the new tree 
  c->indexUsed = m;
  c->index = indexPtr;
  c->dim = dim;

  K = r->K;
   
  // do we have too many points? 
  if(!useProb) {
    if( m <= r->leafSize ) {
  
      // save the final pointer locations 
      for( i = 0; i < m; i++) { 
        // go through each element of indexPtr and store a pointer to that indexPtr element in pointerIndex
        r->pointerIndex[ indexPtr[i] ] = &( indexPtr[i] );
        r->nodeIndex[ indexPtr[i] ] = *nodeIdentity;
//        printf(" node assignment %d is %d\n", (int) indexPtr[i], (int) *nodeIdentity );
      }
      *nodeIdentity = *nodeIdentity + 1;
      return c;
    }
  } else {
  // if using probSize we want to figure out how many samples per psu
    for( i = 0; i < m; i++) {
      probSum += prob[ indexPtr[i] ];
      // go through each element of indexPtr and store a pointer to that indexPtr element in pointerIndex
    } 
    if(probSum <= r->leafSize) {
      for( i = 0; i < m; i++) { 
        r->pointerIndex[ indexPtr[i] ] = &( indexPtr[i] );
        r->nodeIndex[ indexPtr[i] ] = *nodeIdentity;
//        printf(" node assignment %d is %d\n", (int) indexPtr[i], (int) *nodeIdentity );
      }
      *nodeIdentity = *nodeIdentity + 1;
      return c;
    }
#ifdef DEBUG_PROB  
    printf("split!\n");
#endif
  } 

  // if we are here we have too many points 
  // create children
  // figure out our new dim
  // split data and give to children 

  if( useProb ) { 
    c->split = splitDataProb( 
      r->data,
      c->index, 
      &indexLeftPtr,
      &indexRightPtr,
      &indexLeftSize,
      &indexRightSize,
      m, 
      K,
      dim,
      prob 
      ); 
#ifdef DEBUG_PROB  
    printf("Left Side Size = %d, Right Side Size = %d split = %f\n", (int) indexLeftSize, (int) indexRightSize, c->split);
    printf("Left\n:");
    for(i=0; i < indexLeftSize; i++) printf("%d ", (int) indexLeftPtr[i]);
    printf("\nRight\n:");
    for(i=0; i < indexRightSize; i++) printf("%d ", (int) indexRightPtr[i]);
    printf("\n");
#endif
  } else {
    c-> split = splitData( 
      r->data,
      c->index, 
      &indexLeftPtr,
      &indexRightPtr,
      &indexLeftSize,
      &indexRightSize,
      m, 
      K,
      dim  
      ); 
  }

  free(c->index);
  c->index = NULL; 

  // move current contents to new children
  c->left  = buildIndex( r, (dim+1) % K, indexLeftSize , indexLeftPtr,  useProb, prob, nodeIdentity);
  c->right = buildIndex( r, (dim+1) % K, indexRightSize, indexRightPtr, useProb, prob, nodeIdentity);

  return c;
}




// function to create a simple node 
nodePtr createNode( rootNodePtr r ) {

  nodePtr c = malloc( sizeof( node ) );

  c->index = NULL;
  c->left  = NULL;
  c->right = NULL;
  c->indexUsed = 0;
  c->split = 0;
  c->dim = 0;

  return(c);
}




// function to delete a node 
void deleteNode( rootNodePtr r, nodePtr c ) {

  if( c == NULL ) return;

  if( c->index != NULL) free(c->index);

  // to the left
  if( c->left != NULL ) {
    deleteNode( r, c->left );
    c->left = NULL;
  }

  // to the right
  if( c->right != NULL ) {
    deleteNode( r, c->right );
    c->right = NULL;
  }

  free(c);
  c = NULL;
  return;
}


int compDblPtr ( const void * aPtr, const void * bPtr ) {
   double a = **(double **) aPtr;
   double b = **(double **) bPtr;
   if (a < b) return -1;
   if (a > b) return  1;
   return 0;
}



// split and create children using prob weighted median
double splitDataProb( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRightSize,
    size_t n, 
    size_t p,
    size_t dim,
    double * prob
    ) {

  double split;
  double prob_sum,pi;
  size_t i;

  // get the median 
  double * x =  calloc( n, sizeof(double) );        //allocate some temporary space for finding the median
  double ** xPtr =  calloc( n, sizeof( double * ) ); //allocate some temporary space for finding the median
  
  // create input for qsort
  for( i = 0; i < n; i++) {
    x[i] = y[ index[i] * p + dim];
    xPtr[i] = &(x[i]);
  }
  
  // do qsort 
  qsort( xPtr, n, sizeof( double * ), compDblPtr);

  // calculate split point for prob_sum
  prob_sum = 0;
  for( i = 0; i < n ; i++) prob_sum += prob[i]; 
  pi = prob_sum / 2.0;

  *indexLeftSize = 0;
  prob_sum = 0;
  for( i = 0; i < n ; i++) {
    if( prob_sum >= pi ) break;  
    *indexLeftSize = *indexLeftSize + 1; 
    prob_sum += prob[xPtr[i]-x];
  }

  split=*xPtr[i-1];
      
  *indexRightSize = n - *indexLeftSize; 

  if( *indexLeftSize > 0 )
    *indexLeft  = calloc(*indexLeftSize , sizeof(size_t) );
  if( *indexRightSize > 0 )
    *indexRight = calloc(*indexRightSize , sizeof(size_t) );
    
  // now let's have some fun with pointer math 
  for( i = 0; i < *indexLeftSize ; i++) (*indexLeft)[i]  = index[xPtr[i] - x];
  for( i = 0; i < *indexRightSize; i++) (*indexRight)[i] = index[xPtr[*indexLeftSize + i] - x];
  
  free(xPtr);
  free(x); 

  return split;
}




// split and create children
double splitData( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRightSize,
    size_t n, 
    size_t p,
    size_t dim 
    ) {

  double split;
  size_t splitIndex,i;

  // get the median 
  double * x =  calloc( n, sizeof(double) );        //allocate some temporary space for finding the median
  double ** xPtr =  calloc( n, sizeof( double * ) ); //allocate some temporary space for finding the median
  
  // create input for qsort
  for( i = 0; i < n; i++) {
    x[i] = y[ index[i] * p + dim];
    xPtr[i] = &(x[i]);
  }
  // use quick sort to find the median 
  // get split
  splitIndex = n/2;

  // return split value
  split = quantile_quickSelectIndex( xPtr, splitIndex, n ); 
    
  *indexLeftSize  = n / 2;
  *indexRightSize = n - *indexLeftSize;
  
  *indexLeft  = calloc(*indexLeftSize , sizeof(size_t) );
  *indexRight = calloc(*indexRightSize , sizeof(size_t) );
    
  // now let's have some fun with pointer math 
  for( i = 0; i < *indexLeftSize ; i++) (*indexLeft)[i]  = index[xPtr[i] - x];
  for( i = 0; i < *indexRightSize; i++) (*indexRight)[i] = index[ xPtr[*indexLeftSize + i] -x ];
  
  free(xPtr);
  free(x); 

  return split;
}




// a function to print the tree 
void printTree( rootNodePtr r, nodePtr c ) {

  size_t i;

  PRINTF("node: %p\n", (void *) c);
  if( c->index != NULL) {
    for( i=0; i < c->indexUsed; i++) PRINTF("%d ", (int) c->index[i]); 
  } 
  PRINTF("\n\tleft: %p right %p (split = %f) \n", (void *) c->left, (void*) c->right, c->split );
  PRINTF("\n");

  if( c->left ) {
    PRINTF("left ");
    printTree( r, c->left);
  }

  if( c->right ) {
    PRINTF("right ");
    printTree( r, c->right);
  }

}



// a function to print the tree 
void printTree2( rootNodePtr r, nodePtr c, double * splitPointLower, double * splitPointUpper ) {

  size_t i;
  double upper_temp;
  double lower_temp;

  if( !c->left  & !c->right) {
    PRINTF("node: %d\n", (int) r->nodeIndex[c->index[0]] );
    if( c->index != NULL) 
      for( i=0; i < r->K; i++) PRINTF("%d: %f, %f\n",(int) i , splitPointLower[i], splitPointUpper[i]); 
    PRINTF("\n");
    return;
  }

  if( c->left ) {
    upper_temp = splitPointUpper[c->dim];

    splitPointUpper[c->dim] = c->split;
    printTree2( r, c->left, splitPointLower, splitPointUpper);

    splitPointUpper[c->dim] = upper_temp;
  }
  
  if( c->right ) {
    lower_temp = splitPointLower[c->dim];

    splitPointLower[c->dim] = c->split;
    printTree2( r, c->right, splitPointLower, splitPointUpper);
    
    splitPointLower[c->dim] = lower_temp;
  }
  
  return;
}



// a function to write the tree 
void recordBounds( rootNodePtr r, nodePtr c, double * splitPointLower, double * splitPointUpper, double * bound ) {

  size_t i,row;
  double upper_temp;
  double lower_temp;

  if( !c->left  & !c->right) {
    if( c->index != NULL) {
      row = r->nodeIndex[c->index[0]]; 
      for( i=0; i < r->K; i++) { 
        bound[row * 2 * r->K + i] = splitPointLower[i];
        bound[row * 2 * r->K + r->K + i] = splitPointUpper[i];
      } 
      return;
    }
  }

  if( c->left ) {
    upper_temp = splitPointUpper[c->dim];
    splitPointUpper[c->dim] = c->split;
    recordBounds( r, c->left, splitPointLower, splitPointUpper, bound);
    splitPointUpper[c->dim] = upper_temp;
  }
  
  if( c->right ) {
    lower_temp = splitPointLower[c->dim];
    splitPointLower[c->dim] = c->split;
    recordBounds( r, c->right, splitPointLower, splitPointUpper, bound);
    splitPointLower[c->dim] = lower_temp;
  }
  
  return;
}



// function to find the minimal Euclidian distance 
size_t getClosestTie( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint,
    double * dist, 
    double * tieBreak 
  ) {

  size_t i,j,d;
  size_t closestIndex = r->n;

  size_t K = r->K;
  double * x = r->data;            
  double currentDist;
  double newTieBreak;
  double tmp;

  // iterate over all obs on the leaf 
  for( i = 0; i < c->indexUsed; i++) {

    // get first index 
    j = c->index[i]; 
    
    // check if it's a valid index 
    if( j  >= r->n) continue;  
    
    // don't match what we are not looking for
    if( j == query ) continue; 

    // calculate distance
    for( d=0, currentDist=0; d < K; d++) { 
      tmp = x[j * K + d] - queryPoint[d];
      tmp *= tmp;
      currentDist += tmp;
    }
   
    // if smaller than prior index update 
    if( currentDist < *dist ) {
      *dist = currentDist; 
      closestIndex = i;
      *tieBreak = -1;  // set no tie

    // if it is a tie
    } else if( currentDist == *dist ) {
    
      // generate a deviate on (0,1) and pick the biggest
      // each obs has the same prob of being largest due
      // to exchangability  
      newTieBreak = RUNIF;

      if( *tieBreak < 0 ) *tieBreak = RUNIF;  // if no tie was set

      if( newTieBreak > *tieBreak) *tieBreak = newTieBreak;
      closestIndex = i;
    }
  }

  // return the new index if a new NN is found
  if( closestIndex < r->n ) 
    return( c->index[closestIndex] );

  // if no NN is found return n
  return( r->n );
}






// find the nearest neighbor that is not a specific index 
// bound should be first set to the value of the node you 
// are trying to find a neighbor for 
size_t find_nn_notMe( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint, 
    double * dist, 
    double * tieBreak
  ) {

  double distMin;
  size_t nn = r->n;
  size_t nnTmp;
  
  // return if c == NULL 
  if( c == NULL ) return (nn);

  // nothing to search for 
  if( query >= (r->n) ) return (nn);
 
  // is there anything here ? 
  // if there is we get the closest item 
  if( c->index != NULL ) 
    return( getClosestTie( r, c, query, queryPoint, dist, tieBreak) ); 
   
  // first check if the query point is less than split 
  if( queryPoint[c->dim] <= c->split ) {
      nn = find_nn_notMe( r, c->left, query, queryPoint, dist, tieBreak );  

      // now check if there is a point in the split that can be close 
      distMin = queryPoint[c->dim] - c->split;
      distMin *= distMin;
      if(distMin < *dist) {
        nnTmp=find_nn_notMe( r, c->right, query, queryPoint, dist, tieBreak ); 
        // if we found a new closer update our neighbor
        if( nnTmp != r->n) nn = nnTmp;
      } 
  // the query point is greater than the split   
  } else {
      nn = find_nn_notMe( r, c->right, query, queryPoint, dist, tieBreak );  
      
      // now check if there is a point in the split that can be close 
      distMin = queryPoint[c->dim] - c->split;
      distMin *= distMin;
      if(distMin < *dist) {
        nnTmp = find_nn_notMe( r, c->left, query, queryPoint, dist, tieBreak );  
        if( nnTmp != r->n) nn = nnTmp;
      } 
    }


  return (nn);
}





// find the nearest neighbor that is not a specific index 
// bound should be first set to the value of the node you 
// are trying to find a neighbor for 
size_t find_nn_notMe_count( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint, 
    double * dist, 
    double * tieBreak,
    size_t * count,
    size_t * maxCount
  ) {

  double distMin;
  size_t nn = r->n;
  size_t nnTmp;
  
  // return if c == NULL 
  if( c == NULL ) return (nn);

  // nothing to search for 
  if( query >= (r->n) ) return (nn);

  // exceeded count
  if( *count >= *maxCount ) return(nn);
 
  // is there anything here ? 
  // if there is we get the closest item 
  if( c->index != NULL ) { 
    nn = getClosestTie( r, c, query, queryPoint, dist, tieBreak);
    if( nn < r->n) *count = *count + 1; 
    return(nn);
  }
   
  // first check if the query point is less than split 
  if( queryPoint[c->dim] <= c->split ) {
      nn = find_nn_notMe_count( r, c->left, query, queryPoint, dist, tieBreak, count, maxCount );  
      // now check if there is a point in the split that can be close 
      distMin = queryPoint[c->dim] - c->split;
      distMin *= distMin;
      if(distMin < *dist) {
        nnTmp=find_nn_notMe_count( r, c->right, query, queryPoint, dist, tieBreak,count, maxCount ); 
        // if we found a new closer update our neighbor
        if( nnTmp != r->n) nn = nnTmp;
      } 
  // the query point is greater than the split   
  } else {
      nn = find_nn_notMe_count( r, c->right, query, queryPoint, dist, tieBreak,count, maxCount );  
      
      // now check if there is a point in the split that can be close 
      distMin = queryPoint[c->dim] - c->split;
      distMin *= distMin;
      if(distMin < *dist) {
        nnTmp = find_nn_notMe_count( r, c->left, query, queryPoint, dist, tieBreak,count, maxCount );  
        if( nnTmp != r->n) nn = nnTmp;
      } 
    }


  return (nn);
}




// find the nearest neighbor that is not a specific index 
// bound should be first set to the value of the node you 
// are trying to find a neighbor for 
size_t find_nn_notMe_dist( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint, 
    double * dist, 
    double * tieBreak,
    double * termDist 
  ) {

  double distMin;
  size_t nn = r->n;
  size_t nnTmp;
  
  // return if c == NULL 
  if( c == NULL ) return (nn);

  // nothing to search for 
  if( query >= (r->n) ) return (nn);

  // exceeded count
  //printf(" dist = %f < tesrmDist = %f (nn=%zu)\n", *dist, *termDist, nn);
  if( *dist < *termDist ) return(nn);
 
  // is there anything here ? 
  // if there is we get the closest item 
  if( c->index != NULL ) 
    return( getClosestTie( r, c, query, queryPoint, dist, tieBreak) );
  
   
  // first check if the query point is less than split 
  if( queryPoint[c->dim] <= c->split ) {
      nn = find_nn_notMe_dist( r, c->left, query, queryPoint, dist, tieBreak, termDist);  
      // now check if there is a point in the split that can be close 
      distMin = queryPoint[c->dim] - c->split;
      distMin *= distMin;
      if(distMin < *dist) {
        nnTmp=find_nn_notMe_dist( r, c->right, query, queryPoint, dist, tieBreak,termDist ); 
        // if we found a new closer update our neighbor
        if( nnTmp != r->n) nn = nnTmp;
      } 
  // the query point is greater than the split   
  } else {
      nn = find_nn_notMe_dist( r, c->right, query, queryPoint, dist, tieBreak,termDist );  
      
      // now check if there is a point in the split that can be close 
      distMin = queryPoint[c->dim] - c->split;
      distMin *= distMin;
      if(distMin < *dist) {
        nnTmp = find_nn_notMe_dist( r, c->left, query, queryPoint, dist, tieBreak,termDist );  
        if( nnTmp != r->n) nn = nnTmp;
      } 
    }


  return (nn);
}




/* test code */
#ifdef CLI 
int main () {

  double x[20] = {
    0.68512126, 0.3251399, // 0 8
    0.05171296, 0.7740967, // 1 7 
    0.74261974, 0.9242472, // 2 4
    0.13036790, 0.3870030, // 3 6 
    0.77980495, 0.7918827, // 4 2
    0.01413735, 0.5849822, // 5 1 
    0.25770368, 0.4773944, // 6 3 
    0.09543018, 0.8095111, // 7 1 
    0.39014922, 0.3908506, // 8 6 
    0.32050716, 0.1994035  // 9 8 
  };


  size_t ncol=2;
  size_t nrow=10;
  size_t i;
  size_t j;
  double * queryPoint;
  double dist;
  double tieBreak;
  size_t nodeIdentity = 0;

  rootNodePtr myTree = NULL;

  myTree = createTree(ncol, 20, nrow, x);

  // this will get freed
  size_t * index = calloc( nrow, sizeof(size_t));
  for( i = 0; i < nrow; i++) index[i] = i;

  myTree->root = buildIndex( 
    myTree,      // root pointer 
    0,           // current dim
    nrow,        // current length of obs
    index,       // pointer to obs indexes 
    0,
    NULL,
    &nodeIdentity;
  ); 

  printTree( myTree, myTree->root );
 
  for( i = 0; i < nrow; i++) { 
    queryPoint = &(myTree->data[i*ncol]);
    dist = INFINITY;
    tieBreak = -1;
    j = find_nn_notMe(myTree, myTree->root, i, queryPoint, &dist, &tieBreak );  
//    printf("%zu: %zu %f\n", i, j, dist );
  }

  deleteTree( myTree );

  return 0;
}

#endif



