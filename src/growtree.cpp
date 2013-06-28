#include "bigrfc.h"



/*******************************************************************************
 * C interface functions to be called from R. These functions are wrappers to
 * the C++ functions that follow.
 */
extern "C" {
    
    #define CALL_GROWTREE(xtype) {                                            \
        return growtree<xtype>(x, a, aOut, forestP, insampP, inweightP,        \
            treenum, trace);                                                   \
    }
    
    SEXP growtreeC(SEXP xP, SEXP xtypeP, SEXP aP, SEXP aOutP, SEXP forestP,
        SEXP insampP, SEXP inweightP, SEXP treenumP, SEXP traceP) {

        // Initialize function arguments.
        BigMatrix *x = (BigMatrix*)R_ExternalPtrAddr(xP);
        BigMatrix *a = (BigMatrix*)R_ExternalPtrAddr(aP);
        BigMatrix *aOut = (BigMatrix*)R_ExternalPtrAddr(aOutP);
        const int xtype = *INTEGER(xtypeP), treenum = *INTEGER(treenumP),
            trace = *INTEGER(traceP);
        
        switch (xtype) {
            case 1:
                CALL_GROWTREE(char);
                break;
            case 2:
                CALL_GROWTREE(short);
                break;
            case 4:
                CALL_GROWTREE(int);
                break;
            case 8:
                CALL_GROWTREE(double);
                break;
            default:
                return R_NilValue;
        }
    }
    
}



/*******************************************************************************
 * C++ functions.
 */

/* Growtree consists of repeated calls to findbestsplit and movedata.
 * Findbestsplit does just that--it finds the best split of the current node.
 * Movedata moves the data in the split node right and left so that the data
 * corresponding to each child node is contiguous.
 *
 * The growtree bookkeeping is different from that in Friedman's original CART
 * program: 
 *     nnodes is the total number of nodes to date.
 *     treemap[k, ] = child node numbers if the kth node has been split.
 * 	                  -1 if the node exists but has not yet been split.
 *	                  0 if the node is terminal.
 * 
 * A node is terminal if its size is below a threshold value, or if it is all
 * one class, or if all the x-values are equal. If the current node k is split,
 * then its children are numbered nnodes+1 [left], and nnodes+2 [right], nnodes
 * increases to nnodes+2 and the next node to be split is numbered k+1. When no
 * more nodes can be split, growtree returns to the main program.
 */
template <typename xtype>
SEXP growtree(BigMatrix *x, BigMatrix *a, BigMatrix *aOut, SEXP forestP,
    SEXP insampP, SEXP inweightP, int treenum, int trace) {
        
    // Initialize function arguments.
    MatrixAccessor<xtype> xAcc(*x);
    MatrixAccessor<int> aAcc(*a);
    const int nexamples = *INTEGER(GET_SLOT(forestP, install("nexamples")));
    const int *factorvars = INTEGER(GET_SLOT(forestP, install("factorvars")));
    const int *varnlevels = INTEGER(GET_SLOT(forestP, install("varnlevels")));
    const int *varselect = INTEGER(GET_SLOT(forestP, install("varselect")));
    const int nvar = LENGTH(GET_SLOT(forestP, install("varselect")));
    const int *contvarseq = INTEGER(GET_SLOT(forestP, install("contvarseq")));
    const int *y = INTEGER(GET_SLOT(forestP, install("y")));
    const int ynclass = LENGTH(GET_LEVELS(GET_SLOT(forestP, install("y"))));
    const double *yclasswts =
        REAL(GET_SLOT(forestP, install("yclasswts")));
    const int nsplitvar = *INTEGER(GET_SLOT(forestP, install("nsplitvar")));
    const int maxndsize = *INTEGER(GET_SLOT(forestP, install("maxndsize")));
    const int maxeslevels = *INTEGER(GET_SLOT(forestP, install("maxeslevels")));
    const int nrandsplit = *INTEGER(GET_SLOT(forestP, install("nrandsplit")));
    
    // Any numeric variables?
    bool haveNumericVar = false;
    for (int i = 0; i < nvar; i++) {
        if (!factorvars[i]) {
            haveNumericVar = true;
        }
    }
    
    // Theoretically, 2*nexamples - 1 should be enough
    const int maxnodes = 2 * nexamples + 1;
    
    // Initialize tree object.
    SEXP treeP;
    PROTECT(treeP = NEW_OBJECT(MAKE_CLASS("bigctree")));
    SEXP insampS = install("insamp");
    SEXP inweightS = install("inweight");
    SEXP nnodesS = install("nnodes");
    SEXP treemapS = install("treemap");
    SEXP nodeclassS = install("nodeclass");
    SEXP nodewtS = install("nodewt");
    SEXP bestvarS = install("bestvar");
    SEXP bestnumsplitS = install("bestnumsplit");
    SEXP bestcatsplitS = install("bestcatsplit");
    SEXP termincountS = install("termincount");
    SEXP trainprednodeS = install("trainprednode");
    SEXP trainpredclassS = install("trainpredclass");
    SEXP tginiS = install("tgini");
    SET_SLOT(treeP, insampS, duplicate(insampP));
    SET_SLOT(treeP, inweightS, duplicate(inweightP));
    SET_SLOT(treeP, nnodesS, NEW_INTEGER(1));
    SET_SLOT(treeP, treemapS, allocMatrix(INTSXP, maxnodes, 2));
    SET_SLOT(treeP, nodeclassS, NEW_INTEGER(maxnodes));
    SET_SLOT(treeP, nodewtS, NEW_NUMERIC(maxnodes));
    SET_SLOT(treeP, bestvarS, NEW_INTEGER(maxnodes));
    SET_SLOT(treeP, bestnumsplitS, NEW_NUMERIC(maxnodes));
    SET_SLOT(treeP, bestcatsplitS, NEW_LIST(maxnodes));
    SET_SLOT(treeP, termincountS, NEW_INTEGER(maxnodes));
    SET_SLOT(treeP, trainprednodeS, NEW_INTEGER(nexamples));
    SET_SLOT(treeP, trainpredclassS, NEW_INTEGER(nexamples));
    SET_SLOT(treeP, tginiS, NEW_NUMERIC(nvar));
    int *insamp = INTEGER(GET_SLOT(treeP, insampS));
    double *inweight = REAL(GET_SLOT(treeP, inweightS));
    int *nnodes = INTEGER(GET_SLOT(treeP, nnodesS));
    int *treemap1 = INTEGER(GET_SLOT(treeP, treemapS));
    int *treemap2 = treemap1 + maxnodes;
    int *nodeclass = INTEGER(GET_SLOT(treeP, nodeclassS));
    double *nodewt = REAL(GET_SLOT(treeP, nodewtS));
    int *bestvar = INTEGER(GET_SLOT(treeP, bestvarS));
    double *bestnumsplit = REAL(GET_SLOT(treeP, bestnumsplitS));
    SEXP bestcatsplitP = GET_SLOT(treeP, bestcatsplitS);
    int *termincount = INTEGER(GET_SLOT(treeP, termincountS));
    int *trainprednode = INTEGER(GET_SLOT(treeP, trainprednodeS));
    int *trainpredclass = INTEGER(GET_SLOT(treeP, trainpredclassS));
    double *tgini = REAL(GET_SLOT(treeP, tginiS));
    *nnodes = 0;
    for (int i = 0; i < maxnodes; i++) {
        treemap1[i] = -1;
        treemap2[i] = -1;
        nodeclass[i] = 0;
        nodewt[i] = 0;
        bestvar[i] = 0;
        bestnumsplit[i] = 0;
        termincount[i] = 0;
    }
    for (int i = 0; i < nexamples; i++) {
        trainprednode[i] = 0;
        trainpredclass[i] = 0;
    }
    for (int v = 0; v < nvar; v++) {
        tgini[v] = 0;
    }
    
    // Initialize working variables.
    int maxnlevels = 0;
    for (int j = 0; j < nvar; j++) {
        if (varnlevels[j] > maxnlevels) {
            maxnlevels = varnlevels[j];
        }
    }
    int *ndstart = (int*) R_alloc(maxnodes, sizeof(int)),
        ndendl = -1,
        *ndend = (int*) R_alloc(maxnodes, sizeof(int)),
        *ncase = (int*) R_alloc(nexamples, sizeof(int)),
        *ndstartOut = (int*) R_alloc(maxnodes, sizeof(int)),
        ndendlOut = -1,
        *ndendOut = (int*) R_alloc(maxnodes, sizeof(int)),
        *ncaseOut = (int*) R_alloc(nexamples, sizeof(int));
    double **classpop = (double**) R_alloc(ynclass, sizeof(double*)),
        *dgini = (double*) R_alloc(maxnodes, sizeof(double));
    for (int c = 0; c < ynclass; c++) {
        classpop[c] = (double*) R_alloc(maxnodes, sizeof(double));
        for (int n = 0; n < maxnodes; n++) {
            classpop[c][n] = 0;
        }
    }
    ndstart[0] = 0;
    ndend[0] = -1;
    ndstartOut[0] = 0;
    ndendOut[0] = -1;
    for (int i = 0; i < nexamples; i++) {
        if (insamp[i]) {
            int cl = y[i] - 1;
            classpop[cl][0] += insamp[i] * yclasswts[cl];
            ndend[0]++;
            ncase[ndend[0]] = i;
        } else {
            ndendOut[0]++;
            ncaseOut[ndendOut[0]] = i;
        }
    }
    // Variables for storing split results.
    int splitBestvar, splitNumsplit,
        *splitCatsplit = (int*) R_alloc(maxnlevels, sizeof(int));
    double splitDecgini,
        *ndclasspop = (double*) R_alloc(ynclass, sizeof(double));
;
    
    // Main loop.
    int nd;
    for (nd = 0; nd < maxnodes; nd++) {
        if (trace >= 2)
            Rprintf("Tree %d: Processing node %d.\n", treenum, nd + 1);
        
        if (nd > *nnodes) {
            break;
        }
        
        // Initialize.
        for (int c = 0; c < ynclass; c++) {
            ndclasspop[c] = classpop[c][nd];
        }
        
        // If this node has not been split, find the best split.
        if (treemap1[nd] == -1) {
            int stat = findbestsplit<xtype>(x, y, a, factorvars, varnlevels,
                maxnlevels, nvar, varselect, contvarseq, ynclass, nsplitvar,
                maxeslevels, nrandsplit, ncase, inweight, ndstart[nd],
                ndend[nd], ndclasspop, &splitBestvar, &splitDecgini,
                &splitNumsplit, splitCatsplit);
            if (stat) {
                treemap1[nd] = 0;
                treemap2[nd] = 0;
            }
        }
        
        // If this is a terminal node, calculate some statistics and continue.
        if (treemap1[nd] == 0) {
            termincount[nd] = 0;
            for (int i = ndstart[nd]; i <= ndend[nd]; i++) {
                termincount[nd] += insamp[ncase[i]];
            }
            
            double ndclasspopTot = 0, max = 0;
            for (int c = 0; c < ynclass; c++) {
                ndclasspopTot += ndclasspop[c];
                if (ndclasspop[c] > max) {
                    max = ndclasspop[c];
                    nodeclass[nd] = c + 1;
                }
            }
            nodewt[nd] = ndclasspopTot / termincount[nd];
            
            for (int i = ndstart[nd]; i <= ndend[nd]; i++) {
                trainprednode[ncase[i]] = nd + 1;
                trainpredclass[ncase[i]] = nodeclass[nd];
            }
            
            if (ndendOut[nd] - ndstartOut[nd] + 1 > 0) {
                for (int i = ndstartOut[nd]; i <= ndendOut[nd]; i++) {
                    trainprednode[ncaseOut[i]] = nd + 1;
                    trainpredclass[ncaseOut[i]] = nodeclass[nd];
                }
            }
            continue;
        }
        
        // Extract split results.
        bestvar[nd] = splitBestvar + 1;
        dgini[nd] = splitDecgini;
        if (!factorvars[splitBestvar]) {
            // Continuous. Take the split translate it back into x-values
            xtype *xCol = xAcc[varselect[splitBestvar] - 1];
            int *aCol = aAcc[contvarseq[splitBestvar] - 1];
            bestnumsplit[nd] = (xCol[aCol[splitNumsplit] - 1] +
                xCol[aCol[splitNumsplit + 1] - 1]) / 2;
        } else {
            // Categorical
            int lcat = varnlevels[splitBestvar];
            SET_VECTOR_ELT(bestcatsplitP, nd, NEW_INTEGER(lcat));
            int *bestcatsplit = INTEGER(VECTOR_ELT(bestcatsplitP, nd));
            for (int l = 0; l < lcat; l++) {
                bestcatsplit[l] = splitCatsplit[l];
            }
        }
        
        // Move data in a and ncase to reflect split.
        movedata<xtype>(x, a, nexamples, factorvars, haveNumericVar, varselect,
            contvarseq, ndstart[nd], &ndendl, ndend[nd], ncase, splitBestvar,
            splitNumsplit, splitCatsplit);
        if (ndendOut[nd] - ndstartOut[nd] + 1 > 0) {
            movedataOut<xtype>(x, aOut, nexamples, factorvars, haveNumericVar,
                varselect, contvarseq, ndstartOut[nd], &ndendlOut, ndendOut[nd],
                ncaseOut, splitBestvar, bestnumsplit[nd], splitCatsplit);
        }

        // Update tree data. The left node will be nnodes + 1, and the right
        // node will be nnodes + 2.
        ndstart[*nnodes + 1] = ndstart[nd];
        ndend[*nnodes + 1] = ndendl;
        ndstart[*nnodes + 2] = ndendl + 1;
        ndend[*nnodes + 2] = ndend[nd];
        if (ndendOut[nd] - ndstartOut[nd] + 1 > 0) {
            ndstartOut[*nnodes + 1] = ndstartOut[nd];
            ndendOut[*nnodes + 1] = ndendlOut;
            ndstartOut[*nnodes + 2] = ndendlOut + 1;
            ndendOut[*nnodes + 2] = ndendOut[nd];
        } else {
            ndstartOut[*nnodes + 1] = ndstartOut[nd];
            ndendOut[*nnodes + 1] = ndendOut[nd];
            ndstartOut[*nnodes + 2] = ndstartOut[nd];
            ndendOut[*nnodes + 2] = ndendOut[nd];
        }
        
        // Find class populations in both nodes.
        for (int i = ndstart[nd]; i <= ndendl; i++) {
            int nc = ncase[i];
            int c = y[nc] - 1;
            classpop[c][*nnodes + 1] += inweight[nc];
        }
        for (int i = ndendl + 1; i <= ndend[nd]; i++) {
            int nc = ncase[i];
            int c = y[nc] - 1;
            classpop[c][*nnodes + 2] += inweight[nc];
        }
        
        // Check for terminal child nodes.
        int nd1classes = 0, nd2classes = 0;
        for (int c = 0; c < ynclass; c++) {
            if (classpop[c][*nnodes + 1] > 0) {
                nd1classes++;
            }
            if (classpop[c][*nnodes + 2] > 0) {
                nd2classes++;
            }
        }
        if ((ndend[*nnodes + 1] - ndstart[*nnodes + 1] + 1) <= maxndsize ||
            nd1classes == 1) {
            treemap1[*nnodes + 1] = 0;
            treemap2[*nnodes + 1] = 0;
        }
        if ((ndend[*nnodes + 2] - ndstart[*nnodes + 2] + 1) <= maxndsize ||
            nd2classes == 1) {
            treemap1[*nnodes + 2] = 0;
            treemap2[*nnodes + 2] = 0;
        }
        
        // Update treemap.
        treemap1[nd] = *nnodes + 1 + 1;
        treemap2[nd] = *nnodes + 2 + 1;
        *nnodes += 2;
        // if (*nnodes >= maxnodes - 1) {
        //   break;
        // }
        
        // Check for R interrupt.
        R_CheckUserInterrupt();
    }
    
    // Increment nnodes by 1 because R starts counting at 1.
    (*nnodes)++;
    
    // Calculate out-of-bag estimates.
    for (int n = 0; n < *nnodes; n++) {
        if (treemap1[n]> 0) {
            int v = bestvar[n] - 1;
            tgini[v] += dgini[n] * ((ndend[n] - ndstart[n] + 1) +
                (ndendOut[n] - ndstartOut[n] + 1));
        }
    }
    for (int v = 0; v < nvar; v++) {
        tgini[v] /= nexamples;
    }
    
    // Truncate vectors and matrices to nnodes rows, to save memory.
    for (int i = 0; i < *nnodes; i++) {
        treemap1[*nnodes + i] = treemap2[i];
    }
    SETLENGTH(GET_SLOT(treeP, treemapS), 2 * *nnodes);
    SEXP treemapDimP = GET_DIM(GET_SLOT(treeP, treemapS));
    int* treemapDim = INTEGER(treemapDimP);
    treemapDim[0] = *nnodes;
    treemapDim[1] = 2;
    SEXP treemapDimnamesP, treemapDimnamesNamesP;
    PROTECT(treemapDimnamesP = NEW_LIST(2));
    PROTECT(treemapDimnamesNamesP = NEW_CHARACTER(2));
    SET_STRING_ELT(treemapDimnamesNamesP, 0, mkChar("Tree"));
    SET_STRING_ELT(treemapDimnamesNamesP, 1, mkChar("Child"));
    SET_NAMES(treemapDimnamesP, treemapDimnamesNamesP);
    SET_VECTOR_ELT(treemapDimnamesP, 1, NEW_CHARACTER(2));
    SEXP treemapDim2NamesP = VECTOR_ELT(treemapDimnamesP, 1);
    SET_STRING_ELT(treemapDim2NamesP, 0, mkChar("Left"));
    SET_STRING_ELT(treemapDim2NamesP, 1, mkChar("Right"));
    SET_DIMNAMES(GET_SLOT(treeP, treemapS), treemapDimnamesP);
    SETLENGTH(GET_SLOT(treeP, nodeclassS), *nnodes);
    SETLENGTH(GET_SLOT(treeP, nodewtS), *nnodes);
    SETLENGTH(GET_SLOT(treeP, bestvarS), *nnodes);
    SETLENGTH(GET_SLOT(treeP, bestnumsplitS), *nnodes);
    SETLENGTH(GET_SLOT(treeP, bestcatsplitS), *nnodes);
    SETLENGTH(GET_SLOT(treeP, termincountS), *nnodes);
    
    // Return.
    UNPROTECT(3);
    return treeP;
}



/* Finds the best split at the current node. For the best split, bestvar is the
   variable split on. decsplit is the decrease in impurity. If bestvar is
   numerical, nsplit is the example number of value of bestvar split on, and
   nsplitnext is the example number of the next larger value of bestvar. If
   bestvar is categorical, then nsplit is the coding into an integer of the
   categories going left. */
template <typename xtype>
int findbestsplit(BigMatrix *x, const int *y, BigMatrix *a,
    const int *factorvars, const int *varnlevels, int maxnlevels, int nvar,
    const int *varselect, const int *contvarseq, int ynclass, int nsplitvar,
    int maxeslevels, int nrandsplit, const int *ncase, const double *inweight,
    int ndstart, int ndend, const double *ndclasspop, int *bestvar,
    double *decsplit, int *nbest, int *ncatsplit) {
    
    // Mark top of memory stack.
    void *vmax = vmaxget();
    
    // Initialize function arguments.
    MatrixAccessor<xtype> xAcc(*x);
    MatrixAccessor<int> aAcc(*a);
    
    // Set up working variables.
    double tol = sqrt(DOUBLE_EPS);
    double **tclasscat = (double**) R_alloc(ynclass, sizeof(double*));
    for (int c = 0; c < ynclass; c++) {
        tclasscat[c] = (double*) R_alloc(maxnlevels, sizeof(double));
    }
    int *icat = (int*) R_alloc(maxnlevels, sizeof(int));
    double *wl = (double*) R_alloc(ynclass, sizeof(double)),
        *wr = (double*) R_alloc(ynclass, sizeof(double));
    
    // Set up return variables.
    int stat = 0;
    *bestvar = -1;
    *nbest = -1;
    double critmax = 0;
    *decsplit = 0;
    int *ncatsplitbest = (int*) R_alloc(maxnlevels, sizeof(int));
        
    // Initialize random number generator.
    GetRNGstate();

    // Compute initial values of numerator and denominator of Gini, and Gini.
    double pno = 0;
    double pdo = 0;
    for (int c = 0; c < ynclass; c++) {
        pno += R_pow_di(ndclasspop[c], 2);
        pdo += ndclasspop[c];
    }
    double crit0 = pno / pdo;
    
    // Start main loop through variables to find best split.
    for (int v = 0; v < nsplitvar; v++) {
        int var = ftrunc(runif(0, nvar));
        
        if (factorvars[var] == 0) {
            // Continuous variable. Loop through each possible split and find
            // the one that will give the greatest decrease in gini.
            int varX = varselect[var] - 1;
            xtype *xCol = xAcc[varX];
            int varA = contvarseq[var] - 1;
            int *aCol = aAcc[varA];
            
            // Get baseline measures before any splits.
        	double rrn = pno;
			double rrd = pdo;
			double rln = 0;
			double rld = 0;
            for (int c = 0; c < ynclass; c++) {
                wl[c] = 0;
                wr[c] = ndclasspop[c];
            }
            
            // Split one example at a time, starting with smallest values of
            // this variable, and find the split with the largest crit.
            for (int i = ndstart; i <= (ndend - 1); i++) {
                int n = aCol[i] - 1;
                double nw = inweight[n];
                int ny = y[n] - 1;
				rln += nw * (2 * wl[ny] + nw);
				rrn += nw * (-2 * wr[ny] + nw);
				rld += nw;
				rrd -= nw;
				wl[ny] += nw;
				wr[ny] -= nw;
                // Check if this split was better than others before.
  				if (xCol[n] < xCol[aCol[i + 1] - 1]) {
					if(rld > tol && rrd > tol) {
						double crit = (rln / rld) + (rrn / rrd);
						if (crit > critmax) {
    						*bestvar = var;
							*nbest = i;
							critmax = crit;
						}
					}
  				}
            }
        } else {
            // Categorical variable.
            // Compute the decrease in impurity given by categorical splits.
            int varX = varselect[var] - 1;
            xtype *xCol = xAcc[varX];

            // Initialize.
    		int lcat = varnlevels[var];
            for (int c = 0; c < ynclass; c++) {
                for (int l = 0; l < lcat; l++) {
                    tclasscat[c][l] = 0;
                }
            }
            
            // Determine total weights of samples in each class and factor
            // level.
            for (int i = ndstart; i <= ndend; i++) {
    			int n = ncase[i];
                int nl = xCol[n] - 1;
                tclasscat[y[n] - 1][nl] += inweight[n];
            }
            
            // If this categorical variable has more levels than maxeslevels,
            // take the best of nrandsplit random splits. Otherwise, perform an
            // exhaustive search for the best split over all partitions of the
            // variable levels.
            unsigned long maxlsearch;
            if (lcat <= maxeslevels) {
                // Exhaustive search.
                maxlsearch = R_pow_di(2, lcat - 1) - 1;
            } else {
                // Random splits.
                maxlsearch = nrandsplit;
            }
            // Test each split (random or exhustive).
            for (unsigned long n = 1; n <= maxlsearch; n++) {
                // Set up icat to indicate which variables to split on.
                if (lcat <= maxeslevels) {
                    // Exhaustive search.
                    unpack(n, icat, lcat);
                } else {
                    // Generate random split. icat[k] is bernouilli.
                    for (int l = 0; l < lcat; l++) {
                        icat[l] = rbinom(1, 0.5);
                    }
                }
                
                // Try split and calculate decrease in error.
                double pln = 0;
                double pld = 0;
                double prn = 0;
                for (int c = 0; c < ynclass; c++) {
            		double tweight = 0;
                    for (int l = 0; l < lcat; l++) {
                        if (icat[l] == 1) {
                            tweight += tclasscat[c][l];
                        }
                    }
                	pln += R_pow_di(tweight, 2);
            		pld += tweight;
                	tweight = ndclasspop[c] - tweight;
        			prn += R_pow_di(tweight, 2);
                }
                
                // Check if this split was better than others before.
                if (pld > tol && pdo > pld + tol) {
                    double tdec = (pln / pld) + (prn / (pdo - pld));
                    if (tdec > critmax) {
            			*bestvar = var;
                        critmax = tdec;
                        for (int l = 0; l < lcat; l++) {
                            ncatsplitbest[l] = icat[l];
                        }
                    }
                }
            }
        }
    }
    
    // Calculate decrease in gini with the best split.
    *decsplit = critmax - crit0;
    if (critmax == 0) {
        stat = 1;
    }
    
    // If the best split was categorical, retrieve the split.
    if (factorvars[*bestvar] == 1) {
        int lcat = varnlevels[*bestvar];
        for (int l = 0; l < lcat; l++) {
            ncatsplit[l] = ncatsplitbest[l];
        }
    }
    
    // Clean up after random number generation.
    PutRNGstate();
    
    // Release memory.
    vmaxset(vmax);
    
    // Return.
    return stat;    
}



/* Moves the data in the current node to the left and right children, according
   to the best split on the current node. Modifies a, ndendl and ncase. */
template <typename xtype>
void movedata(BigMatrix *x, BigMatrix *a, int nexamples, const int *factorvars,
    bool haveNumericVar, const int *varselect, const int *contvarseq,
    int ndstart, int *ndendl, int ndend, int *ncase, int bestvar, int nbest,
    const int *bestcatsplit) {

    // Mark top of memory stack.
    void *vmax = vmaxget();

    // Initialize function arguments.
    MatrixAccessor<xtype> xAcc(*x);
    MatrixAccessor<int> aAcc(*a);
    
    // Set up working variables.
    int bestvarX = varselect[bestvar] - 1, bestvarA = contvarseq[bestvar] - 1;
    
    // compute idmove = indicator of example nos. going left
    int *idmove = (int*) R_alloc(nexamples, sizeof(int));
    if (factorvars[bestvar] == 0) {
        int *aCol = aAcc[bestvarA];
        for (int j = ndstart; j <= nbest; j++) {
            idmove[aCol[j] - 1] = 1;
        }
        for (int j = nbest + 1; j <= ndend; j++) {
            idmove[aCol[j] - 1] = 0;
        }
        *ndendl = nbest;
    } else {
        xtype *xCol = xAcc[bestvarX];
        *ndendl = ndstart - 1;
        for (int j = ndstart; j <= ndend; j++) {
            int nc = ncase[j];
            if (bestcatsplit[(int)xCol[nc] - 1] == 1) {
                idmove[nc] = 1;
                (*ndendl)++;
            } else {
                idmove[nc] = 0;
            }
        }
    }
    
    movedataWorker(a, factorvars, haveNumericVar, contvarseq, ndstart, ndend,
        idmove, ncase, bestvar, bestvarA);
        
    // Release memory.
    vmaxset(vmax);
}



/* Similar to movedata, except for out-of-bag samples. Modifies a, ndendl and
   ncase. */
template <typename xtype>
void movedataOut(BigMatrix *x, BigMatrix *a, int nexamples,
    const int *factorvars, bool haveNumericVar, const int *varselect,
    const int *contvarseq, int ndstart, int *ndendl, int ndend, int *ncase,
    int bestvar, double bestnumsplit, const int *bestcatsplit) {
    
    // Mark top of memory stack.
    void *vmax = vmaxget();

    // Initialize function arguments.
    MatrixAccessor<xtype> xAcc(*x);
    MatrixAccessor<int> aAcc(*a);
            
    // Set up working variables.
    int bestvarX = varselect[bestvar] - 1, bestvarA = contvarseq[bestvar] - 1;
    
    // compute idmove = indicator of example nos. going left
    int *idmove = (int*) R_alloc(nexamples, sizeof(int));
    if (!factorvars[bestvar]) {
        xtype *xCol = xAcc[bestvarX];
        int *aCol = aAcc[bestvarA];
        *ndendl = ndstart - 1;
        for (int j = ndstart; j <= ndend; j++) {
            int nc = aCol[j] - 1;
            if (xCol[nc] <= bestnumsplit) {
                idmove[nc] = 1;
                (*ndendl)++;
            } else {
                idmove[nc] = 0;
            }
        }
    } else {
        xtype *xCol = xAcc[bestvarX];
        *ndendl = ndstart - 1;
        for (int j = ndstart; j <= ndend; j++) {
            int nc = ncase[j];
            if (bestcatsplit[(int)xCol[nc] - 1] == 1) {
                idmove[nc] = 1;
                (*ndendl)++;
            } else {
                idmove[nc] = 0;
            }
        }
    }
    
    if (!(*ndendl < ndstart || *ndendl == ndend)) {
        // At least one example is split from the rest of the group.
        movedataWorker(a, factorvars, haveNumericVar, contvarseq, ndstart,
            ndend, idmove, ncase, bestvar, bestvarA);
    }
    
    // Release memory.
    vmaxset(vmax);
}



/* Worker function for movedata and movedataOut that does the actual moving. */
void movedataWorker(BigMatrix *a, const int *factorvars, bool haveNumericVar,
    const int *contvarseq, int ndstart, int ndend, const int *idmove,
    int *ncase, int bestvar, int bestvarA) {
    
    // Initialize function arguments.
    MatrixAccessor<int> aAcc(*a);
    int *tmp = (int*) R_alloc(ndend - ndstart + 1, sizeof(int));
    
    // shift example nos. right and left for numerical variables.
    if (haveNumericVar) {
        for (int i = 0; i < a->ncol(); i++) {
            int *aCol = aAcc[i];
            int k = 0;
            for (int j = ndstart; j <= ndend; j++) {
                if (idmove[aCol[j] - 1] == 1) {
                    tmp[k++] = aCol[j];
                }
            }
            for (int j = ndstart; j <= ndend; j++) {
                if (idmove[aCol[j] - 1] == 0) {
                    tmp[k++] = aCol[j];
                }
            }
            k = 0;
            for (int j = ndstart; j <= ndend; j++) {
                aCol[j] = tmp[k++];
            }
        }
    }
    
    // Compute example nos. for right and left nodes.
    if (!factorvars[bestvar]) {
        int *aCol = aAcc[bestvarA];
        for (int j = ndstart; j <= ndend; j++) {
            ncase[j] = aCol[j] - 1;
        }
    } else {
        int k = 0;
        for (int j = ndstart; j <= ndend; j++) {
            if (idmove[ncase[j]] == 1) {
                tmp[k++] = ncase[j];
            }
        }
        for (int j = ndstart; j <= ndend; j++) {
            if (idmove[ncase[j]] == 0) {
                tmp[k++] = ncase[j];
            }
        }
        k = 0;
        for (int j = ndstart; j <= ndend; j++) {
            ncase[j] = tmp[k++];
        }
    }
}



/* Utility function for "unpacking"" an integer into an array of 1s and 0s (i.e.
   binary). Used to represent best categorical splits (ncatsplit). Requires
   an empty array icat to be passed in. */
void unpack(unsigned long npack, int *icat, int l) {
    icat[0] = npack % 2;
    for (int k = 1; k < l; k++) {
        npack = (npack - icat[k - 1]) / 2;
        icat[k] = npack % 2;
    }
}



/* The reverse of unpack. */
unsigned long pack(int *icat, int l) {
    unsigned long n = 0;
    for (int i = l - 1; i >= 0; i--) {
        n *= 2;
        n += icat[i];
    }
    return n;
}
