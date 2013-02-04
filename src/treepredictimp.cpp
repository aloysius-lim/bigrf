#include "bigrfc.h"



/*******************************************************************************
 * C interface functions to be called from R. These functions are wrappers to
 * the C++ functions that follow.
 */
extern "C" {
    
    #define CALL_TREEPREDICTIMP(xtype) {                                       \
        return treepredictimp<xtype>(x, nout, whichout, whichpermout,          \
            permutevar, forestP, treeP);                                       \
    }
    
    SEXP treepredictimpC(SEXP xP, SEXP xtypeP, SEXP noutP, SEXP whichoutP,
        SEXP whichpermoutP, SEXP permutevarP, SEXP forestP, SEXP treeP) {
        
        // Initialize function arguments.
        BigMatrix *x = (BigMatrix*)R_ExternalPtrAddr(xP);
        int xtype = *INTEGER(xtypeP), nout = *INTEGER(noutP),
            *whichout = INTEGER(whichoutP),
            *whichpermout = INTEGER(whichpermoutP),
            permutevar = *INTEGER(permutevarP) - 1;
        
        switch (xtype) {
            case 1:
                CALL_TREEPREDICTIMP(char);
                break;
            case 2:
                CALL_TREEPREDICTIMP(short);
                break;
            case 4:
                CALL_TREEPREDICTIMP(int);
                break;
            case 8:
                CALL_TREEPREDICTIMP(double);
                break;
            default:
                return R_NilValue;
        }
    }
    
}



/*******************************************************************************
 * C++ functions.
 */

/* Predicts classification for out-of-bag examples in training set, with the
   specified variable randomly permuted. */
template <typename xtype>
SEXP treepredictimp(BigMatrix* x, int nout, const int *whichout,
    const int *whichpermout, int permutevar, SEXP forestP, SEXP treeP) {
    
    // Initialize function arguments.
    MatrixAccessor<xtype> xAcc(*x);
    const int *factorvars = INTEGER(GET_SLOT(forestP, install("factorvars")));
    const int *varselect = INTEGER(GET_SLOT(forestP, install("varselect")));
    const int nnodes = *INTEGER(GET_SLOT(treeP, install("nnodes")));
    const int *bestvar = INTEGER(GET_SLOT(treeP, install("bestvar")));
    const int *treemap1 = INTEGER(GET_SLOT(treeP, install("treemap")));
    const int *treemap2 = treemap1 + nnodes;
    const double *bestnumsplit = REAL(GET_SLOT(treeP, install("bestnumsplit")));
    SEXP bestcatsplit = GET_SLOT(treeP, install("bestcatsplit"));
    const int *nodeclass = INTEGER(GET_SLOT(treeP, install("nodeclass")));
    
    // Initialize return variables.
    SEXP retP, retNamesP;
    PROTECT(retP = NEW_LIST(2));
    SET_VECTOR_ELT(retP, 0, NEW_INTEGER(nout));
    SET_VECTOR_ELT(retP, 1, NEW_INTEGER(nout));
    PROTECT(retNamesP = NEW_CHARACTER(2));
    SET_STRING_ELT(retNamesP, 0, mkChar("oobpredclass"));
    SET_STRING_ELT(retNamesP, 1, mkChar("oobprednode"));
    SET_NAMES(retP, retNamesP);
    // Out-of-bag predictions.
    int *oobpredclass = INTEGER(VECTOR_ELT(retP, 0));
    // Final node that each oub-of-bag example ended up in, in a particular tree.
    int *oobprednode = INTEGER(VECTOR_ELT(retP, 1));
    
    // Test each out-of-bag example.
    for (int i = 0; i < nout; i++) {
        // Go down the tree until we reach a terminal node.
        int nd = 0;
        while (treemap1[nd] > 0) {
            int bv = bestvar[nd] - 1;
            double xVal;
            if (bv == permutevar) {
                xVal = xAcc[varselect[bv] - 1][whichpermout[i] - 1];
            } else {
                xVal = xAcc[varselect[bv] - 1][whichout[i] - 1];
            }
            if (!factorvars[bv]) {
                // This node was split on a numerical variable.
                if (xVal <= bestnumsplit[nd]) {
                    nd = treemap1[nd] - 1;
                } else {
                    nd = treemap2[nd] - 1;
                }
            } else {
                // This node was split on a categorical variable.
                if (INTEGER(VECTOR_ELT(bestcatsplit, nd))[(int)xVal - 1] == 1) {
                    nd = treemap1[nd] - 1;
                } else {
                    nd = treemap2[nd] - 1;
                }
            }
        }
		oobpredclass[i] = nodeclass[nd];
		oobprednode[i] = nd + 1;
    }
    
    // Return.
    UNPROTECT(2);
    return retP;
}