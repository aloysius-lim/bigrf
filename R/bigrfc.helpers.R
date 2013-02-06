# ------------------------------------------------------------------------------
# Converts input variables x into a big.matrix. Argument x can be a matrix or
# data.frame.
makex <- function(x, backingfile="", cachepath=NULL) {
    # Get appropriate type for new big.matrix.
    if (class(x) == "matrix") {
        xtype <- switch(typeof(x), double="double", integer="integer",
                        logical="char", NULL)
        if (is.null(xtype)) {
            stop("Matrix x can only be a numeric, integer or logical ",
                 "matrix.")
        }
    } else if (class(x) == "data.frame") {
        xclasses <- sapply(x, class)
        if (any(xclasses == "numeric")) {
            xtype <- "double"
        } else if (any(xclasses %in% c("integer", "factor"))) {
            xtype <- "integer"
        } else if (any(xclasses == "logical")) {
            xtype <- "char"
        } else {
            stop("Data.frame x can only contain numeric, integer, factor ",
                 "or logical data.")
        }
    }
    
    # Create big.matrix.
    if (is.null(cachepath)) {
        xnew <- big.matrix(nrow(x), ncol(x), type=xtype,
                           dimnames=list(NULL, colnames(x)))
    } else {
        xnew <- big.matrix(nrow(x), ncol(x), type=xtype,
                           dimnames=list(NULL, colnames(x)),
                           backingfile=backingfile,
                           descriptorfile=paste0(backingfile, ".desc"),
                           backingpath=cachepath)
    }
    
    # Copy data.
    old.opt <- options(bigmemory.typecast.warning=FALSE)
    for (j in seq_len(ncol(x))) {
        if (xtype %in% c("integer", "char")) {
            xnew[, j] <- as.integer(x[, j])
        } else {
            xnew[, j] <- x[, j]
        }
    }
    options(old.opt)
    
    return(xnew)
}



# ------------------------------------------------------------------------------
# makea constructs the nexamples x nvarx integer matrix a. For each numerical 
# variable with values x[n,m],n=1,...,nexamples, the x-values are sorted from 
# lowest to highest. Denote these by xs[n,m]. Then asave[n,m] is the example
# number in which xs[n,m] occurs. If the mth variable is categorical, then
# asave[n,m] is the category of the nth example number. asave is a big.matrix
# passed by reference.
makea <- function(x, asave, factorvars, varselect) {
    v5 <- numeric(length(varselect))
    v95 <- numeric(length(varselect))
    
    for (var in which(!factorvars)) {
        asave[, var] <- order(x[, varselect[var]])
        # v5[var] <- quantile(x[, varselect[var]], probs=0.05)
        # v95[var] <- quantile(x[, varselect[var]], probs=0.95)
    }
    for (var in which(factorvars)) {
        asave[, var] <- as.integer(x[, varselect[var]])
    }
    
    # return(list(v5=v5, v95=v95))
    return()
}



# ------------------------------------------------------------------------------
moda <- function(asave, a, factorvars, insamp) {
    return(.Call("modaC", asave@address, a@address, factorvars,
                 as.integer(insamp), PACKAGE="bigrf"))
}



# ------------------------------------------------------------------------------
# Buildtree consists of repeated calls to findbestsplit and movedata. 
# Findbestsplit does just that--it finds the best split of the current node.
# Movedata moves the data in the split node right and left so that the data
# corresponding to each child node is contiguous.
# 
# The buildtree bookkeeping is different from that in Friedman's original CART
# program: 
#     nnodes is the total number of nodes to date.
# 	treemap[k, ] = child node numbers if the kth node has been split.
# 	               -1 if the node exists but has not yet been split.
# 		              0 if the node is terminal.
# 
# A node is terminal if its size is below a threshold value, or if it is all one
# class, or if all the x-values are equal. If the current node k is split, then
# its children are numbered nnodes+1 [left], and nnodes+2 [right], nnodes
# increases to nnodes+2 and the next node to be split is numbered k+1. When no
# more nodes can be split, buildtree returns to the main program.
buildtree <- function(x, y, asave, a, a.out, forest, insamp, inweight, treenum,
                      trace) {
    xtype <- as.integer(.Call("CGetType", x@address, PACKAGE="bigmemory"))
    return(.Call("buildtreeC", x@address, xtype, y, asave@address, a@address,
                 a.out@address, forest, insamp, inweight, treenum, trace))
}



# ------------------------------------------------------------------------------
# Combine results of tree builds. To be used only as a .combine function in
# foreach().
combine.treeresults <- function(forest, newtree) {
    treenum <- forest@ntrees + 1L
    oldntrees <- newtree$oldntrees
    ntrees <- newtree$ntrees
    y <- newtree$y
    tree <- newtree$tree
    printerrfreq <- newtree$printerrfreq
    printclserr <- newtree$printclserr
    rm(newtree)
    
    forest[[treenum]] <- tree
    forest@ntrees <- treenum
    
    forest@oobtimes[tree@insamp == 0L] <-
        forest@oobtimes[tree@insamp == 0L] + 1L
    
    # Get out-of-bag estimates -------------------------------------------------
    
    for (c in seq_len(forest@ynclass)) {
        # Out-of-bag examples with votes for this class.
        w <- which(tree@trainpredclass == c & tree@insamp == 0L)
        forest@oobvotes[w, c] <- forest@oobvotes[w, c] +
            tree@nodewt[tree@trainprednode[w]]
    }
    rm(c, w)
    
    # Get training set error estimates -----------------------------------------
    
    forest@oobpred[forest@oobtimes > 0L] <-
        max.col(forest@oobvotes[forest@oobtimes > 0L, ])
    
    for (c in seq_len(forest@ynclass)) {
        forest@trainclserr[treenum, c] <- sum(y == c & forest@oobpred != c)
    }
    
    forest@trainerr[treenum] <- sum(forest@trainclserr[treenum, ]) /
        forest@nexamples
    forest@trainclserr[treenum, ] <- forest@trainclserr[treenum, ] /
        as.numeric(forest@ytable)
    
    # Accumulate Gini decreases for each variable ------------------------------
    
    forest@varginidec <- forest@varginidec + tree@tgini
    
    # Give running output ------------------------------------------------------
    
    if (treenum == oldntrees + 1L) {
        cat("OOB errors:\n")
        if (printclserr) {
            cat(" Tree  Overall error  Error by class\n")
            cat("                      ")
            cat(format(names(forest@ytable), justify="right", width=5),
                sep="  ")
            cat("\n")
        } else {
            cat(" Tree  Overall error\n")
        }
    }
    
    if ((treenum - oldntrees) %% printerrfreq == 0L ||
            treenum == oldntrees + ntrees) {
        cat(format(treenum, justify="right", width=5),
            format(100 * forest@trainerr[treenum], justify="right", width=13,
                   digits=3, nsmall=2), sep="  ")
        if (printclserr) {
            cat("",
                format(100 * forest@trainclserr[treenum, ], justify="right",
                       width=max(nchar(forest@ylevels), 5), digits=3, nsmall=2),
                sep="  ")
        }
        cat("\n")
    }
    
    return(forest)
}



# ------------------------------------------------------------------------------
# Combine results of tree builds. To be used only as a .combine function in
# foreach().
combine.treepredictresults <- function(prediction, treepredict.result) {
    y <- treepredict.result$y
    forest <- treepredict.result$forest
    tree <- treepredict.result$tree
    t <- treepredict.result$t
    printerrfreq <- treepredict.result$printerrfreq
    printclserr <- treepredict.result$printclserr
    
    # Compute votes.
    for (c in seq_len(prediction@ynclass)) {
        w <- which(treepredict.result$testpredclass == c)
        prediction@testvotes[w, c] <- prediction@testvotes[w, c] +
            tree@nodewt[treepredict.result$testprednode[w]]
    }
    rm(c, w)
    prediction[seq_len(prediction@ntest)] <- max.col(prediction@testvotes)
    
    # If test set labels were given, compute test error.
    if (!is.null(y)) {
        prediction@testclserr <- integer(prediction@ynclass)
        for (c in seq_len(prediction@ynclass)) {
            prediction@testclserr[c] <-
                sum(y == c & prediction[] != c)
        }
        prediction@testerr <- sum(prediction@testclserr) / prediction@ntest
        prediction@testclserr <-
            prediction@testclserr / as.numeric(prediction@testytable)
    }
    
    # Give running output --------------------------------------------------
    if (t == 1L) {
        if (printclserr && !is.null(y)) {
            cat("Test errors:\n")
            cat(" Tree  Overall error  Error by class\n")
            cat("                      ")
            cat(format(names(prediction@testytable), justify="right",
                       width=5),
                sep="  ")
            cat("\n")
        } else {
            cat("Processing tree number:\n")
        }
    }
    
    if (t %% printerrfreq == 0L || t == forest@ntrees) {
        cat(format(t, justify="right", width=5))
        if (!is.null(y) && printclserr) {
            cat("",
                format(100 * prediction@testerr, justify="right", width=13,
                       digits=3, nsmall=2),
                format(100 * prediction@testclserr, justify="right",
                       width=max(nchar(names(prediction@testytable)), 5),
                       digits=3,
                       nsmall=2),
                sep="  ")
        }
        cat("\n")
    }
    
    return(prediction)
}



# -------------------------------------------------------
# preprox <- function(near, maxnodes, ntrees, treemap, ncount, treenum, trainprednode, nodexb,
#                     insamp, jinb, ndbegin, rinpop, npcase, termincount) {
#   nod <- integer(maxnodes)
#   ncn <- integer(near)
#   ncount[, treenum] <- 0L
#   ndbegin[, treenum] <- 0L
#   
# 	ntt <- 0L
# 	for (k in seq_len(maxnodes)) {
# 		if (nodestatus[k, 1] == 0L) {
# 			ntt <- ntt + 1L
# 			nod[k] <- ntt
# 		}
# 	}
# 	nterm <- ntt
# 	for (n in seq_len(near)) {
# 		rinpop[n, treenum] <- termincount[trainprednode[n]]
# 		nodexb[n, treenum] <- nod[trainprednode[n]]
# 		jinb[n, treenum] <- insamp[n]
# 		k <- nodexb[n, treenum]
# 		ncount[k, treenum] <- ncount[k, treenum] + 1L
# 		ncn[n] <- ncount[k, treenum]
# 	}
# 	ndbegin[1, treenum] <- 1L
# 	for (k in 2:(nterm + 1)) {
# 		ndbegin[k, treenum] <- ndbegin[k - 1, treenum] + ncount[k - 1, treenum]
# 	}
# 	for (n in seq_len(near)) {
# 		kn <- ndbegin[nodexb[n, treenum], treenum] + ncn[n] - 1L
# 		npcase[kn, treenum] <- n
# 	}
# }



# -------------------------------------------------------
# comprox <- function(prox, nodexb, jinb, ndbegin, npcase, rinpop, near, ntrees,
#                     noutlier, y, loz, nrnn, wtx) {
#   outtr <- numeric(near)
#   
# 	for (n in seq_len(near)) {
#     ppr <- numeric(near)
#     for (treenum in seq_len(ntrees)) {
#       k <- nodexb[n, treenum]
#       if (jinb[n, treenum] > 0L) {
#         for (j in ndbegin[k, treenum]:(ndbegin[k + 1, treenum] - 1L)) {
#           kk <- npcase[j, treenum]
#           if (jinb[kk, treenum] == 0L) {
#             ppr[kk] <- ppr[kk] + (wtx[n] / rinpop[n, treenum])
#           }
#         }
#       }
#       if (jinb[n, treenum] == 0L) {
#         for (j in ndbegin[k, treenum]:(ndbegin[k + 1, treenum] - 1)) {
#           kk <- npcase[j, treenum]
#           if (jinb[kk, treenum] > 0L) {
#             ppr[kk] <- ppr[kk] + (wtx[kk] / rinpop[kk, treenum])
#           }
#         }
#       }
#     }  # ntrees
#     
#     if (noutlier == 1L) {
# 			rsq <- 0L
# 			for (k in seq_len(near)) {
# 			  if (ppr[k] > 0 && y[k] == y[n]) {
# 			    rsq <- rsq + ppr[k] ^ 2
# 			  }
# 			}
# 			if (rsq == 0L) {
#         rsq <- 1L
# 			}
# 			outtr[n] <- near / rsq
# 		}
# 	
# 		if (nrnn == near) {
# 			for (k in seq_len(near)) {
# 				prox[n, k] <- ppr[k]
# 				loz[n, k] <- k
# 			}
# 		} else {
# 			ibest <- biggest(ppr, near, nrnn)
# 			for (k in seq_len(nrnn)) {
# 				prox[n, k] <- ppr[ibest[k]]
# 				loz[n, k] <- ibest[k]
# 			}
# 		}
# 	} # n
#   
#   return(outtr)
# }



# ------------------------------------------------------
# finds the nrnn largest values in the vector x and 
# returns their positions in the vector ibest[1],...,ibest[nrnn]:
#   x[ibest[1]] is the largest
#   ...
# 	x[ibest[nrnn]] is the nrnn-th-largest
# the vector x is not disturbed 
# the vector iwork is used as workspace
# biggest <- function(x, n, nrnn) {
#   iwork <- integer(n)
#   
# 	ihalfn <- as.integer(n / 2)
# 	for (i in seq_len(n)) {
# 		iwork[i] <- i
# 	}
# 	for (j in seq_len(ihalfn)) {
# 		i <- ihalfn + 1L - j
# 		iwork <- sift(x, iwork, n, n, i)
# 	}
# 	for (j in seq_len(nrnn - 1L)) {
# 		i <- n - j + 1L
# 		ibest[j] <- iwork[1]
# 		jsave <- iwork[i]
# 		iwork[i] <- iwork[1]
# 		iwork[1] <- jsave
# 		iwork <- sift(x, iwork, n, n - j, 1L)
# 	}
# 	ibest[nrnn] <- iwork[1]
#   
#   return(ibest)
# }



# ------------------------------------------------------
# used by subroutine biggest,to bring the largest element to the 
# top of the heap
# sift <- function(x, iwork, n, m, i) {
# 	xsave <- x[iwork[i]]
# 	ksave <- iwork[i]
# 	jsave <- i
# 	j <- i + i
# 	for (k in seq_len(m)) {
# 		if  (j > m) {
#       break
# 		}
# 		if (j < m) {
# 			if (x[iwork[j]] < x[iwork[j + 1]]) {
#         j <- j + 1L
# 	    }
# 		}
# 		if (xsave >= x[iwork[j]]) {
#       break
# 		}
# 		iwork[jsave] <- iwork[j]
# 		jsave <- j
# 		j <- j + j
# 	}
# 	iwork[jsave] <- ksave
#   
#   return(iwork)
# }



# # ------------------------------------------------------
# 	subroutine locateout[y,tout,outtr,ncp,isort,devout, 	near,nexamples,ynclass,rmedout]
# # 
# 	real outtr[near],tout[near],devout[ynclass],rmedout[ynclass]
#      	integer y[nexamples],isort[nexamples],ncp[near],near,nexamples, 	ynclass
# 	real rmed, dev
# 	integer jp,nt,n,i
#      
# 	for (jp in seq_len(ynclass) {
# 		nt <- 0
# 		for (n in seq_len(near) {
# 			if (y[n] == jp) {
# 				nt <- nt + 1
# 				tout[nt]=outtr[n]
# 				ncp[nt]=n
# 			}
# 		}
# 		quicksort.result <- quicksort(tout,isort,1,nt,nexamples)
# 		rmed <- tout[[1 + nt]/2]
# 		dev <- 0
# 		for (i in seq_len(nt) {
# 			dev <- dev + amin1[abs[tout[i]-rmed],5 * rmed]
# 		}
# 		dev <- dev / nt
# 		devout[jp]=dev
# 		rmedout[jp]=rmed
# 		for (i in seq_len(nt) {
# 			outtr[ncp[i]]=amin1[[outtr[ncp[i]]-rmed]/dev,20.0]
# 		}
# 	} # jp
# 	end



# # -------------------------------------------------------
# 	subroutine myscale[loz,prox,xsc,y,u,near,nscale,red,nrnn, 	ee,ev,dl]
# # 
# 	double precision prox[near,nrnn],y[near],u[near],dl[nscale], 	xsc[near,nscale],red[near],ee[near],ev[near,nscale],bl[10]
# 	integer loz[near,nrnn]
# # 
# 	integer near,nscale,nrnn,j,i,it,n,jit,k
# 	double precision dotd,y2,sred,eu,ru,ra,ynorm,sa
# 	for (j in seq_len(near) {
# 		ee[j]=dble[1.]
# 	}
# # 
# 	for (j in seq_len(near) {
# 		red[j]=0
# 		for (i in seq_len(nrnn) {
# 			red[j]=red[j]+prox[j,i]
# 		}
# 		red[j]=red[j]/near
# 	}
# 	sred <- dotd[ee,red,near]
# 	sred <- sred / near
# 	for (it in seq_len(nscale) {
# 		for (n in seq_len(near) {
# 			if (mod[n,2] == 0) {
# 				y[n]=1
# 			} else {
# 				y[n]=-1
# 			}
# 		}
# 		for (jit in seq_len(1000) {
# 			y2 <- dotd[y,y,near]
# 			y2 <- dsqrt[y2]
# 			for (n in seq_len(near) {
# 				u[n]=y[n]/y2
# 			}
# 			for (n in seq_len(near) {
# 				y[n]=0
# 				for (k in seq_len(nrnn) {
# 					y[n]=y[n]+prox[n,k]*u[loz[n,k]]
# 				}
# 			}
# 			eu <- dotd[ee,u,near]
# 			ru <- dotd[red,u,near]
# 			for (n in seq_len(near) {
# 				y[n]=y[n]-[red[n]-sred]*eu - ru
# 				y[n]=.5 * y[n]
# 			}
# 			if (it > 1) {
# 				for (j in seq_len(it) {-1
# 					bl[j]=0
# 					for (n in seq_len(near) {
# 						bl[j]=bl[j]+ev[n,j]*u[n]
# 					}
# 					for (n in seq_len(near) {
# 						y[n]=y[n]-bl[j]*dl[j]*ev[n,j]
# 					}
# 				}
# 			}
# 			ra <- dotd[y,u,near]
# 			ynorm <- 0
# 			for (n in seq_len(near) {
# 				ynorm <- ynorm+[y[n]-ra * u[n]]**2
# 			}
# 			sa <- dabs[ra]
# 			if (ynorm < sa * 1.0e - 7) {
# 				for (n in seq_len(near) {
# 					xsc[n,it]=[dsqrt[sa]]*u[n]
# 					ev[n,it]=u[n]
# 				}
# 				dl[it]=ra
# 				goto 101
# 			}
# 		}
# 101		continue
# 	} # nn
# 	end



# # -------------------------------------------------------
# 
# 	subroutine xfill[x,nexamples,nvarx,fill,code]
# # 
# # input:
# 	real code,fill[nvarx]
# 	integer nvarx,nexamples
# # output:
# 	real x[nvarx,nexamples]
# # local:
# 	integer n,m
# 	for (n in seq_len(nexamples) {
# 		for (m in seq_len(nvarx) {
# 			if[abs[x[n, m]-code] < 8.232D - 11]  				x[n, m]=fill[m]
# 		} # m
# 	}
# 	end



# # -------------------------------------------------------
# 	subroutine roughfix[x,v,ncase,nvarx,nexamples,cat,code, 	nrcat,maxcat,fill]
# # 
# 	real x[nvarx,nexamples],v[nexamples],fill[nvarx],code
# 	  integer ncase[nexamples],cat[nvarx],nrcat[maxcat]
# 	integer nvarx,nexamples,maxcat
# 	integer m,n,nt,j,jmax,lcat,nmax
# 	real rmed
# # 
# 	for (m in seq_len(nvarx) {
# 		if (cat[m] == 1) {
# # 		continuous variable
# 			nt <- 0
# 			for (n in seq_len(nexamples) {
# 				if (abs[x[n, m]-code] >= 8.232D - 11) {
# 					nt <- nt + 1
# 					v[nt]=x[n, m]
# 				}
# 			}
# 	   		call quicksort [v,ncase,1,nt,nexamples]
# 			if (nt > 0) {
# 				rmed <- v[[nt + 1]/2]
# 			} else {
# 				rmed <- 0
# 			}
# 			fill[m]=rmed
# 		} else {
# # 		categorical variable
# 			lcat <- cat[m]
# 			nrcat <- integer(maxcat)
# 			for (n in seq_len(nexamples) {
# 				if (abs[x[n, m]-code] >= 8.232D - 11) {
# 					j <- nint[x[n, m]]
# 					nrcat[j]=nrcat[j]+1
# 				}
# 			}
# 			nmax <- 0
# 			jmax <- 1
# 			for (j in seq_len(lcat) {
# 				if (nrcat[j] > nmax) {
# 					nmax <- nrcat[j]
# 					jmax <- j
# 				}
# 			}
# 			fill[m]=real[jmax]
# 		}
# 	} # m
# 	for (n in seq_len(nexamples) {
# 		for (m in seq_len(nvarx) {
# 			if[abs[x[n, m]-code] < 8.232D - 11] x[n, m]=fill[m]
# 		}
# 	}
# 	end



# # -------------------------------------------------------
# 	subroutine impute[x,prox,near,nvarx, 	maxcat,votecat,cat,nrnn,loz,missing]
# # 
# 	real x[nvarx,near],votecat[maxcat]
# 	double precision  prox[near,nrnn]
# 	integer near,nvarx,maxcat,nrnn
# 	integer cat[nvarx],loz[near,nrnn], 	missing[nvarx,near]
# 	integer i,j,jmax,m,n,k
# 	real sx,dt,rmax
# # 
# 	for (m in seq_len(nvarx) {
# 		if (cat[m] == 1) {
# 			for (n in seq_len(near) {
# 				if (missing[n, m] == 1) {
# 					sx <- 0
# 					dt <- 0
# 					for (k in seq_len(nrnn) {
# 						if (missing[loz[n,k], m] != 1) {
# 							sx <- sx + real[prox[n,k]]*x[loz[n,k], m]
# 							dt <- dt + real[prox[n,k]]
# 						}
# 					}
# 					if[dt > 0] x[n, m]=sx / dt
# 				}
# 			} # n
# 		}
# 	} # m
# 	for (m in seq_len(nvarx) {
# 		if (cat[m] > 1) {
# 			for (n in seq_len(near) {
# 				if (missing[n, m] == 1) {
# 					votecat <- numeric(maxcat)
# 					for (k in seq_len(nrnn) {
# 						if (missing[loz[n,k], m] != 1) {
# 							j <- nint[x[loz[n,k], m]]
# 							votecat[j]=votecat[j]+real[prox[n,k]]
# 						}
# 					} # k
# 					rmax=-1
# 					for (i in seq_len(cat] {[m)
# 						if (votecat[i] > rmax) {
# 							rmax <- votecat[i]
# 							jmax <- i
# 						}
# 					}
# 					x[n, m]=real[jmax]
# 				}
# 			} # n
# 		}
# 	} # m
# 	end
# # 
