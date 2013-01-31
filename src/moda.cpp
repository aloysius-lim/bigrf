#include "bigrfc.h"



/*******************************************************************************
 * C interface functions to be called from R. These functions are wrappers to
 * the C++ functions that follow.
 */
extern "C" {
    
    SEXP modaC(SEXP asaveP, SEXP aP, SEXP factorsP, SEXP insampP) {
        return moda(asaveP, aP, factorsP, insampP);
    }
    
}



/*******************************************************************************
 * C++ functions.
 */

/* Prepares the a matrix based on random sample of rows for modelling. For each
   continuous variable, copies only the in-sample indices from asave to a. Data
   for categorical variables are not copied, as they are stored in asave. This
   functions should only be called if there are any continuous variables. */
SEXP moda(SEXP asaveP, SEXP aP, SEXP factorsP, SEXP insampP) {
    // Initialize function arguments.
    BigMatrix *asave = (BigMatrix*)R_ExternalPtrAddr(asaveP);
    BigMatrix *a = (BigMatrix*)R_ExternalPtrAddr(aP);
    MatrixAccessor<int> asaveAcc(*asave);
    MatrixAccessor<int> aAcc(*a);
    int *asaveCol, *aCol;
    int *factors = LOGICAL(factorsP), *insamp = INTEGER(insampP);
    
    // Set up working variables.
    index_type nCols = asave->ncol();
    index_type nRows = asave->nrow();
    index_type ia, ja, ib, jb;
    
    // For each numerical variable, move all the in-sample data to the top rows
    // of a.
    for (ia = 0, ib = 0; ia < nCols; ia++) {
        if (factors[ia] == 0) {
            asaveCol = asaveAcc[ia];
            aCol = aAcc[ib++];
            for (ja = 0, jb = 0; ja < nRows; ja++) {
                if (insamp[asaveCol[ja] - 1] >= 1) {
                    aCol[jb++] = asaveCol[ja];
                }
            }
        }
    }
    return R_NilValue;
}
