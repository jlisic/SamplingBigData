/* Copyright (c) 2015  Jonathan Lisic 
 * Last edit: 15/07/08 - 11:45:28
 * License: GPL (>=2) 
 */  

#ifndef KDTREE_HEADER

#define KDTREE_HEADER


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/* create a node */
struct node {
  size_t dim;         /* dim split across */
  size_t * index;     /* index rows for x */
  size_t indexUsed;   /* since we are using arrays, this indicates the number of elements of each array */
  double split;       /* point split on for dim */
  struct node * left; /* values less than or equal to split */
  struct node * right; /* value greater than or equal to split */ 
};

/* create typedef */
typedef struct node node;
typedef struct node * nodePtr;


/* tree node */
struct rootNode {
  size_t K;                 /* dims of x */
  size_t leafSize;          /* max number of elements on each leaf */
  size_t n;                 /* number of rows in x */ 
  size_t ** pointerIndex;   /* pointer index to allow easy changing of returned values */
  double * data;            /* pointer to x */
  size_t * nodeIndex;            /* pointer to node assignment */
  nodePtr root;             /* root node pointer */
};

/* create typedef */
typedef struct rootNode rootNode;
typedef struct rootNode * rootNodePtr;


/* function to print a tree */
void printTree( rootNodePtr r, nodePtr c ); 

/* function to identify bounds of the tree */
void printTree2( rootNodePtr r, nodePtr c, double * splitPointLower, double * splitPointUpper ); 

/* save bound */
void recordBounds( rootNodePtr r, nodePtr c, double * splitPointLower, double * splitPointUpper, double * bound ); 

/* function to create a new Tree */
rootNodePtr createTree( 
    size_t K,        // number of columns of data
    size_t leafSize, // maximum size of leaf nodes
    size_t n,        // number of rows of data
    double * data    // pointer to the data
  );

/* delete tree */
void deleteTree( rootNodePtr r );

/* build index for the tree */
nodePtr buildIndex( 
    rootNodePtr r,    // root node pointer
    size_t dim,       // dim to base next node on
    size_t m,         // size of index
    size_t * indexPtr,// index array to add to new node
    int useProb,        // determine if we use probability to build an index
    double * prob,
    size_t * nodeIdentity
  ); 

// create Node 
nodePtr createNode( rootNodePtr r );

// delete a node and all of it's children 
void deleteNode( rootNodePtr r, nodePtr c ); 

// split and create children
double splitDataProb( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRighSize,
    size_t n, 
    size_t p,
    size_t dim,
    double * prob 
    ); 

// split and create children
double splitData( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRighSize,
    size_t n, 
    size_t p,
    size_t dim 
    ); 

// split and create children by prob
double splitDataProb( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRighSize,
    size_t n, 
    size_t p,
    size_t dim,
    double * prob 
    ); 

// function to get the closest neighbor, with tie handling 
size_t getClosestTie( 
    rootNodePtr r,    // rootnode ptr
    nodePtr c,        // leaf node to search
    size_t query,      // index of item (used to not pick self)
    double * queryPoint,       // vector of item 
    double * dist,    // current min dist
    double * tieBreak // tie 
    ); 


// function to find neighbors 
size_t find_nn_notMe( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint, 
    double * dist, 
    double * tieBreak  
    ); 

// function to find neighbors fixed count
size_t find_nn_notMe_count( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint, 
    double * dist, 
    double * tieBreak,
    size_t * count,  
    size_t * maxCount
    ); 

// function to find neighbors with min dist 
size_t find_nn_notMe_dist( 
    rootNodePtr r, 
    nodePtr c, 
    size_t query, 
    double * queryPoint, 
    double * dist, 
    double * tieBreak,
    double * termDist  
    ); 

int compDblPtr ( 
    const void * aPtr, 
    const void * bPtr 
  );

#endif
