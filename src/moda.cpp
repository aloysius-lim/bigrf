#include "bigrfc.h"



/*******************************************************************************
 * C interface functions to be called from R. These functions are wrappers to
 * the C++ functions that follow.
 */
extern "C" {
    
    SEXP modaC(SEXP asaveP, SEXP aP, SEXP insampP) {
        return moda(asaveP, aP, insampP);
    }
    
}



/*******************************************************************************
 * C++ functions.
 */

/* Prepares the a matrix based on random sample of examples for modelling. For
   each continuous variable, copies only the in-sample indices from asave to a.
   Data for categorical variables are not copied, as they are stored in x.
   This function should only be called if there are any continuous variables. */
SEXP moda(SEXP asaveP, SEXP aP, SEXP insampP) {
    // Initialize function arguments.
    BigMatrix *asave = (BigMatrix*)R_ExternalPtrAddr(asaveP);
    BigMatrix *a = (BigMatrix*)R_ExternalPtrAddr(aP);
    MatrixAccessor<int> asaveAcc(*asave);
    MatrixAccessor<int> aAcc(*a);
    int *asaveCol, *aCol;
    int *insamp = INTEGER(insampP);
    
    // Set up working variables.
    index_type nCols = asave->ncol();
    index_type nRows = asave->nrow();
    index_type i, ja, jb;
    
    // For each numerical variable, move all the in-sample data to the top rows
    // of a.
    for (i = 0; i < nCols; i++) {
        asaveCol = asaveAcc[i];
        aCol = aAcc[i];
        for (ja = 0, jb = 0; ja < nRows; ja++) {
            if (insamp[asaveCol[ja] - 1] >= 1) {
                aCol[jb++] = asaveCol[ja];
            }
        }
    }
    return R_NilValue;
}
