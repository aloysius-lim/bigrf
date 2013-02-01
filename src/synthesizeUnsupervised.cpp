#include "bigrfc.h"



/*******************************************************************************
 * C interface functions to be called from R. These functions are wrappers to
 * the C++ functions that follow.
 */
extern "C" {
    
    SEXP synthesizeUnsupervisedC(SEXP xP, SEXP xNewP, SEXP xtypeP) {
        int xtype = *INTEGER(xtypeP);
        switch (xtype) {
            case 1:
                synthesizeUnsupervised<char>(xP, xNewP);
                break;
            case 2:
                synthesizeUnsupervised<short>(xP, xNewP);
                break;
            case 4:
                synthesizeUnsupervised<int>(xP, xNewP);
                break;
            case 8:
                synthesizeUnsupervised<double>(xP, xNewP);
                break;
        }
        return R_NilValue;
    }
    
}



/*******************************************************************************
 * C++ functions.
 */

/* Synthesizes second class for unsupervised learning. */
template <typename xtype>
void synthesizeUnsupervised(SEXP xP, SEXP xNewP) {
    
    // Initialize function arguments.
    BigMatrix *x = (BigMatrix*)R_ExternalPtrAddr(xP);
    BigMatrix *xNew = (BigMatrix*)R_ExternalPtrAddr(xNewP);
    MatrixAccessor<xtype> xAcc(*x);
    MatrixAccessor<xtype> xNewAcc(*xNew);
    
    // Initialize working variables.
    int nrows = x->nrow(), ncols = x->ncol();
    
    // Initialize random number generator.
    GetRNGstate();

    for (int i = 0; i < ncols; i++) {
        xtype *xCol = xAcc[i], *xNewCol = xNewAcc[i];
        
        // Copy original values to first half of xNew.
        for (int j = 0; j < nrows; j++) {
            xNewCol[j] = xCol[j];
        }
        // Select random value for synthesized data.
        for (int j = nrows; j < 2 * nrows; j++) {
            xNewCol[j] = xCol[(int)ftrunc(runif(0, nrows))];
        }
    }
    
    // Clean up after random number generation.
    PutRNGstate();
}
