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
SEXP moda(SEXP asaveP, SEXP aP, SEXP insampP);

/* Grows a classification tree. */
template <typename xtype>
SEXP growtree(BigMatrix *x, BigMatrix *a, BigMatrix *aOut, SEXP forestP,
    SEXP insampP, SEXP inweightP, int treenum, int trace);

/* Finds the best split at the current node. */
template <typename xtype>
int findbestsplit(BigMatrix *x, const int *y, BigMatrix *a,
    const int *factorvars, const int *varnlevels, int maxnlevels, int nvar,
    const int *varselect, const int *contvarseq, int ynclass, int nsplitvar,
    int maxeslevels, int nrandsplit, const int *ncase, const double *inweight,
    int ndstart, int ndend, const double *ndclasspop, int *bestvar,
    double *decsplit, int *nbest, int *ncatsplit);
    
/* Moves the data in the current node to the left and right children, according
   to the best split on the current node. */
template <typename xtype>
void movedata(BigMatrix *x, BigMatrix *a, int nexamples, const int *factorvars,
    bool haveNumericVar, const int *varselect, const int *contvarseq,
    int ndstart, int *ndendl, int ndend, int *ncase, int bestvar, int nbest,
    const int *bestcatsplit);

/* Similar to movedata, except for out-of-bag examples. */
template <typename xtype>
void movedataOut(BigMatrix *x, BigMatrix *a, int nexamples,
    const int *factorvars, bool haveNumericVar, const int *varselect,
    const int *contvarseq, int ndstart, int *ndendl, int ndend, int *ncase,
    int bestvar, double bestnumsplit, const int *bestcatsplit);

/* Worker function for movedata and movedataOut that does the actual moving. */
void movedataWorker(BigMatrix *a, const int *factorvars, bool haveNumericVar,
    const int *contvarseq, int ndstart, int ndend, const int *idmove,
    int *ncase, int bestvar, int bestvarA);
    
/* Predicts classification for test set. */
template <typename xtype>
SEXP treepredict(BigMatrix* x, int ntest, SEXP forestP, SEXP treeP);
    
/* Predicts classification for out-of-bag examples in training set, with the
   specified variable randomly permuted. */
template <typename xtype>
SEXP treepredictimp(BigMatrix* x, int nout, const int *whichout,
    const int *whichpermout, int permutevar, SEXP forestP, SEXP treeP);

/* Utility function for "unpacking"" an integer into an array of 1s and 0s (i.e.
   binary). */
void unpack(unsigned long npack, int icat[], int l);

/* The reverse of unpack. */
unsigned long pack(int icat[], int l);

#endif
