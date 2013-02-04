#include "bigrfc.h"



/*******************************************************************************
 * C interface functions to be called from R. These functions are wrappers to
 * the C++ functions that follow.
 */
extern "C" {
    
    #define CALL_TREEPREDICT(xtype) {                                          \
        return treepredict<xtype>(x, ntest, forestP, treeP);                   \
    }
    
    SEXP treepredictC(SEXP xP, SEXP xtypeP, SEXP ntestP, SEXP forestP,
        SEXP treeP) {
        
        // Initialize function arguments.
        BigMatrix *x = (BigMatrix*)R_ExternalPtrAddr(xP);
        int xtype = *INTEGER(xtypeP), ntest = *INTEGER(ntestP);
        
        switch (xtype) {
            case 1:
                CALL_TREEPREDICT(char);
                break;
            case 2:
                CALL_TREEPREDICT(short);
                break;
            case 4:
                CALL_TREEPREDICT(int);
                break;
            case 8:
                CALL_TREEPREDICT(double);
                break;
            default:
                return R_NilValue;
        }
    }
    
}



/*******************************************************************************
 * C++ functions.
 */

/* Predicts classification for test set. */
template <typename xtype>
SEXP treepredict(BigMatrix* x, int ntest, SEXP forestP, SEXP treeP) {
    
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
    SET_VECTOR_ELT(retP, 0, NEW_INTEGER(ntest));
    SET_VECTOR_ELT(retP, 1, NEW_INTEGER(ntest));
    PROTECT(retNamesP = NEW_CHARACTER(2));
    SET_STRING_ELT(retNamesP, 0, mkChar("testpredclass"));
    SET_STRING_ELT(retNamesP, 1, mkChar("testprednode"));
    SET_NAMES(retP, retNamesP);
    // Test set predictions for test cases in a particular tree.
    int *testpredclass = INTEGER(VECTOR_ELT(retP, 0));
    // Final node that each test example ended up in, in a particular tree.
    int *testprednode = INTEGER(VECTOR_ELT(retP, 1));
    
    // Test each test example.
    for (int i = 0; i < ntest; i++) {
        // Go down the tree until we reach a terminal node.
        int nd = 0;
        while (treemap1[nd] > 0) {
            int bv = bestvar[nd] - 1;
            double xVal = xAcc[varselect[bv] - 1][i];
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
		testpredclass[i] = nodeclass[nd];
		testprednode[i] = nd + 1;
    }
    
    // Return.
    UNPROTECT(2);
    return retP;
}