/* Helpers for bigrf classification forests. */
/* This header file contains short descriptions for each function. For more
   details, see the .cpp file. */
#include <iostream>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"

#ifndef BIGRFC_HELPERS_H
#define BIGRFC_HELPERS_H

/* Prepares the a matrix based on random sample of rows for modelling. */
SEXP moda(SEXP asaveP, SEXP aP, SEXP factorsP, SEXP insampP);

/* Builds a classification tree. */
template <typename xtype>
SEXP buildtree(BigMatrix *x, const int *y, BigMatrix *asave, BigMatrix *a,
    BigMatrix *aOut, SEXP forestP, SEXP insampP, SEXP inweightP, int treenum,
    int trace);

/* Finds the best split at the current node. */
template <typename xtype>
int findbestsplit(BigMatrix *x, const int *y, BigMatrix *asave,
    BigMatrix *a, const int *factors, const int *nlevels, int maxnlevels,
    int nvar, const int *varselect, const int *contvarseq, int nclass,
    int nsplitvar, int maxeslevels, int nrandsplit, const int *ncase,
    const double *inweight, int ndstart, int ndend, const double *ndclasspop,
    int *bestvar, double *decsplit, int *nbest, int *ncatsplit);
    
/* Moves the data in the current node to the left and right children, according
   to the best split on the current node. */
void movedata(BigMatrix *asave, BigMatrix *a, int nsample, const int *factors,
    const int *contvarseq, int ndstart, int *ndendl, int ndend, int *ncase,
    int bestvar, int nbest, const int *bestcatsplit);

/* Similar to movedata, except for out-of-bag samples. */
template <typename xType>
void movedataOut(BigMatrix *x, BigMatrix *asave, BigMatrix *a, int nsample,
    const int *factors, const int *varselect, const int *contvarseq,
    int ndstart, int *ndendl, int ndend, int *ncase, int bestvar,
    double bestnumsplit, const int *bestcatsplit);

/* Worker function for movedata and movedataOut that does the actual moving. */
void movedataWorker(MatrixAccessor<int> aAcc, const int *factors,
    const int *contvarseq, int ndstart, int ndend, const int *idmove,
    int *ncase, int bestvar, int bestvarA, int nCols);
    
/* Utility function for "unpacking"" an integer into an array of 1s and 0s (i.e.
   binary). */
void unpack(unsigned long long npack, int icat[], int l);

/* The reverse of unpack. */
unsigned long long pack(int icat[], int l);

#endif
