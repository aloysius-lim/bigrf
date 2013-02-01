# ------------------------------------------------------------------------------
# makea constructs the nexamples x nvarx integer matrix a. For each numerical 
# variable with values x[n,m],n=1,...,nexamples, the x-values are sorted from 
# lowest to highest. Denote these by xs[n,m]. Then asave[n,m] is the example
# number in which xs[n,m] occurs. If the mth variable is categorical, then
# asave[n,m] is the category of the nth example number. asave is a big.matrix
# passed by reference.
makea <- function(x, asave, factors, varselect) {
    v5 <- numeric(length(varselect))
    v95 <- numeric(length(varselect))
    
    for (var in which(!factors)) {
        asave[, var] <- order(x[, varselect[var]])
        # v5[var] <- quantile(x[, varselect[var]], probs=0.05)
        # v95[var] <- quantile(x[, varselect[var]], probs=0.95)
    }
    for (var in which(factors)) {
        asave[, var] <- as.integer(x[, varselect[var]])
    }
    
    # return(list(v5=v5, v95=v95))
    return()
}



# ------------------------------------------------------------------------------
moda <- function(asave, a, factors, insamp) {
    return(.Call("modaC", asave@address, a@address, factors, as.integer(insamp),
                PACKAGE="bigrf"))
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
    trace <- newtree$trace
    treenum <- forest@ntrees + 1L
    oldntrees <- newtree$oldntrees
    ntrees <- newtree$ntrees
    y <- newtree$y
    insamp <- newtree$insamp
    tree <- newtree$tree
    printerrfreq <- newtree$printerrfreq
    printclserr <- newtree$printclserr
    rm(newtree)
    
    forest[[treenum]] <- tree
    forest@ntrees <- treenum
    
    forest@oobtimes[insamp == 0L] <- forest@oobtimes[insamp == 0L] + 1L
    
    # Get out-of-bag estimates -------------------------------------------------
    
    for (c in seq_len(forest@nclass)) {
        # Process for out-of-bag cases for this class.
        w <- which(tree@trainpredclass == c & insamp == 0L)
        forest@oobvotes[w, c] <- forest@oobvotes[w, c] +
            tree@nodewt[tree@trainprednode[w]]
    }
    rm(c, w)
    forest@avgini <- forest@avgini + tree@tgini
    
    # Get training set error estimates -----------------------------------------
    
    forest@oobpred[forest@oobtimes > 0L] <-
        max.col(forest@oobvotes[forest@oobtimes > 0L, ])
    
    for (c in seq_len(forest@nclass)) {
        forest@trainclserr[c] <- sum(y == c & forest@oobpred != c)
    }
    
    forest@trainerr <- sum(forest@trainclserr) / forest@nexamples
    forest@trainclserr <- forest@trainclserr / as.numeric(forest@ytable)
    
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
            format(100 * forest@trainerr, justify="right", width=13,
                   digits=3, nsmall=2), sep="  ")
        if (printclserr) {
            cat("",
                format(100 * forest@trainclserr, justify="right",
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
    for (c in seq_len(prediction@nclass)) {
        w <- which(treepredict.result$testpredclass == c)
        prediction@testvotes[w, c] <- prediction@testvotes[w, c] +
            tree@nodewt[treepredict.result$testprednode[w]]
    }
    rm(c, w)
    prediction[seq_len(prediction@ntest)] <- max.col(prediction@testvotes)
    
    # If test set labels were given, compute test error.
    if (!is.null(y)) {
        prediction@testclserr <- integer(prediction@nclass)
        for (c in seq_len(prediction@nclass)) {
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
# Predicts the class of out-of-bag-cases for variable importance. Also computes
# nodexvr
# testreeimp <- function(x, joob, nout, mr, treemap, bestnumsplit,
#                        bestvar, nodeclass, nnodes, factors, bestcatsplit) {
#   jvr <- integer(nout)
#   nodexvr <- integer(nout)
# 	pjoob <- sample(joob[seq_len(nout)], nout)
# 	n <- 1L
# 	while (n <= nout) {
# 	  kt <- 1L
# 	  for (k in seq_len(nnodes)) {
# 	    if (treemap[kt, 1] == 0L) {
# 	      jvr[n] <- nodeclass[kt]
# 	      nodexvr[n] <- kt
#         next
# 	    }
# 	    m <- bestvar[kt]
# 	    if (varselect[m] == mr) {
# 	      xmn <- x[pjoob[n], varselect[m]] # permuted value
# 	    } else {
# 	      xmn <- x[joob[n], varselect[m]]
# 	    }
# 	    if (!factors[m]) {
# 	      if (xmn <= bestnumsplit[kt]) { 
# 	        kt <- treemap[kt, 1]
# 	      } else {
# 	        kt <- treemap[kt, 2]
# 	      }
# 	    }
# 	    if (factors[m]) {
# 	      jcat <- xmn
# 	      if (bestcatsplit[[kt]][jcat] == 1L) {
# 	        kt <- treemap[kt, 1]
# 	      } else {
# 	        kt <- treemap[kt, 2]
# 	      }
# 	    }
# 	  } # k
# 	  n <- n + 1L
# 	}
#   
#   return(list(jvr=jvr, nodexvr=nodexvr))
# }



# # -------------------------------------------------------
# 	subroutine permobmr[joob,pjoob,nout]
# # 
# # randomly permute the elements of joob and put them in pjoob
# # 
# # input:
# 	integer joob[nout],nout
# # output:
# 	integer pjoob[nout]
# # local:
# 	integer j,k,jt
# 	real rnd,randomu
# # 
# 	for (j in seq_len(nout) {
# 		pjoob[j]=joob[j]
# 	} 
# 	j <- nout
# 11	rnd <- randomu[]
# 	k <- int[j * rnd]
# 	if[k < j] k <- k + 1
# # switch j and k
# 	jt <- pjoob[j]
# 	pjoob[j]=pjoob[k]
# 	pjoob[k]=jt
# 	j <- j - 1
# 	if[j > 1] go to 11
# 	end
# # 



# # -------------------------------------------------------
# 	subroutine comperrts[qts,clts,ntest,nclass,errts,   tmissts,ncts,jests,labelts]
# # 
# # 
# 	integer clts[ntest],ncts[nclass],jests[ntest]
# 	real qts[nclass,ntest],tmissts[nclass]
# 	integer ntest,nclass,labelts
# 	real errts,cmax
# 	integer n,j,jmax
# 	tmissts <- numeric(nclass)
# 	errts <- 0
# 	for (n in seq_len(ntest) {
# 		cmax <- 0
# 		for (j in seq_len(nclass) {
# 			if (qts[j,n] > cmax) {
# 				jmax <- j
# 				cmax <- qts[j,n]
# 			}
# 		}
# 		jests[n]=jmax
# 		if (labelts == 1) {
# 			if (jmax != clts[n]) {
# 				tmissts[clts[n]]=tmissts[clts[n]]+1
# 				errts <- errts + 1
# 			}
# 		}
# 	}
# 	if (labelts == 1) {
# 		errts <- errts / ntest
# 		for (j in seq_len(nclass) {
# 			tmissts[j]=tmissts[j]/ncts[j]
# 		}
# 	}
# 	end



# -------------------------------------------------------
# varimp <- function(x, nexamples, nvarx, y, insamp, trainpredclass, impn, varselect, qimpm,
#                    treemap, bestnumsplit, bestvar, nodeclass,
#                    nnodes, factors, bestcatsplit, nodewt, trainprednode, mimp) {
#   sqsd <- numeric(mimp)
#   avimp <- numeric(mimp)
#   
#   # nout = number of obs out-of-bag for this tree.
#   joob <- which(insamp == 0L)
#   nout <- length(joob)
#   # Update count of correct oob classifications.
#   # trainpredclass[n]=y[n] if example n is correctly classified.
#   w <- which(trainpredclass[joob] == y[joob])
#   right <- sum(nodewt[trainprednode[joob[w]]])
# 
#   if (impn == 1L) {
#     qimp <- numeric(nexamples)
#     qimp[joob[w]] <- nodewt[trainprednode[joob[w]]] / nout
#   }
#   # iv[j]=1 if variable j was used to split on
#   iv <- integer(length(varselect))
#   iv[unique(bestvar[treemap[, 1] > 0L])] <- 1L
#   for (k in seq_len(length(varselect))) {
#     mr <- varselect[k]
#     # choose only those that used a split on variable mr 
#     if (iv[k] == 1L) {
#       testreeimp.result <- testreeimp(x, joob, nout, mr, treemap, 
#                                       bestnumsplit, bestvar, nodeclass, nnodes,
#                                       factors, bestcatsplit)
#       rightimp <- 0
#       for (n in seq_len(nout)) {
#         # the nth out-of-bag example is the nnth original example
#         nn <- joob[n]
#         if (impn == 1L) {
#           if (testreeimp.result$jvr[n] == y[nn]) {
#             qimpm[nn, k] <- qimpm[nn, k] +
#               nodewt[testreeimp.result$nodexvr[n]] / nout
#           }
#         }
#         if(testreeimp.result$jvr[n] == y[nn]) {
#           rightimp <- rightimp + nodewt[testreeimp.result$nodexvr[n]]
#         }
#       }
#       avimp[k] <- avimp[k] + (right - rightimp) / nout
#       sqsd[k] <- sqsd[k] + ((right - rightimp) ^ 2) / (nout ^ 2)
#     } else {
#       for (n in seq_len(nout)) {
#         # the nth out-of-bag example is the nnth original example
#         nn <- joob[n]
#         if (impn == 1L) {
#           if (trainpredclass[nn] == y[nn]) {
#             qimpm[nn, k] <- qimpm[nn, k] + nodewt[trainprednode[nn]] / nout
#           }
#         }
#       }
#     }
#   } # k
#   
#   if (impn == 1L) {
#     return(list(avimp=avimp, sqsd=sqsd, qimp=qimp))
#   } else {
#     return(list(avimp=avimp, sqsd=sqsd))
#   }
# }



# -------------------------------------------------------
# finishimp <- function(nvarx, sqsd, avimp, ntrees, varselect) {
#   zscore <- numeric(length(varselect))
#   signif <- numeric(length(varselect))
#   
# 	for (k in seq_len(length(varselect))) {
# 		m1 <- varselect[k]
# 		avimp[k] <- avimp[k] / ntrees
# 		av <- avimp[k]
# 		se <- (sqsd[k] / ntrees) - av ^ 2
# 		se <- sqrt(se / ntrees)
# 		if (se > 0) {
# 			zscore[k] <- avimp[k] / se
# 			v <- zscore[k]
# 			signif[k] <- erfcc(v)
# 		} else {
# 			zscore[k] <- -5
# 			signif[k] <- 1
# 		}
# 	}
#   
#   return(list(avimp=avimp, zscore=zscore, signif=signif))
# }



# # -------------------------------------------------------
# 	subroutine compinteract[votes,effect,varselect,nvarx,	ntrees,g,iv,irnk,hist,teffect]
# # 
# 	real votes[nvarx,ntrees],effect[nvarx,nvarx],g[nvarx], 	hist[0:nvarx,nvarx],teffect[nvarx,nvarx]
# 	integer varselect[length(varselect)],iv[nvarx],irnk[nvarx,ntrees]
# 	integer nvarx,length(varselect),ntrees,treenum,i,nt,irk,j,ii,   jj,ij,mmin,m,k
# 	real gmin,rcor
# 	for (treenum in seq_len(ntrees) {
# 		nt <- 0
# 		iv <- integer(nvarx)
# 		for (i in seq_len(length(varselect)) {
# 			m <- varselect[i]
# 			g[m]=votes[m,treenum]
# 			if (abs[g[m]] < 8.232D - 11) {
# 				irnk[m,treenum]=0
# 				iv[m]=1
# 				nt <- nt + 1
# 			}
# 		}
# 		irk <- 0
# 		for (j in seq_len(8000) {
# 			gmin <- 10000
# 			for (i in seq_len(length(varselect)) {
# 				m <- varselect[i]
# 				if (iv[m] == 0 && g[m] < gmin) {
# 					gmin <- g[m]
# 					mmin <- m
# 				}
# 			}
# 			iv[mmin]=1
# 			irk <- irk + 1
# 			irnk[mmin,treenum]=irk
# 			nt <- nt + 1
# 			if[nt >= length(varselect)] goto 79
# 		} # j
# 79		continue
# 	} # treenum
# 	for (j in 0:length(varselect)) {
# 		for (i in seq_len(length(varselect)) {
# 			m <- varselect[i]
# 			hist[j,m]=0
# 		}
# 	}
# 	for (i in seq_len(length(varselect)) {
# 		m <- varselect[i]
# 		for (treenum in seq_len(ntrees) {
# 			hist[irnk[m,treenum],m]=hist[irnk[m,treenum],m]+1
# 		}
# 		for (j in 0:length(varselect)) {
# 			hist[j,m]=hist[j,m]/ntrees
# 		}
# 	} # m
# 	
# 	effect <- matrix(0, nvarx, nvarx)
# 	for (i in seq_len(length(varselect)) {
# 		for (j in seq_len(length(varselect)) {
# 			m <- varselect[i]
# 			k <- varselect[j]
# 			for (treenum in seq_len(ntrees) {
# 				effect[m,k]=effect[m,k]+iabs[irnk[m,treenum]-irnk[k,treenum]]
# 			}
# 			effect[m,k]=effect[m,k]/ntrees
# 		}
# 	}
# 	teffect <- matrix(0, nvarx, nvarx)
# 	for (i in seq_len(length(varselect)) {
# 		for (j in seq_len(length(varselect)) {
# 			m <- varselect[i]
# 			k <- varselect[j]
# 			for (ii in 0:length(varselect)) {
# 				for (jj in 0:length(varselect)) {
# 					teffect[m,k]=teffect[m,k]+abs[ii - jj]*hist[jj,m]*hist[ii,k]
# 				}
# 			}
# 			rcor <- 0
# 			for (ij in seq_len(length(varselect)) {
# 				rcor <- rcor + hist[ij,m]*hist[ij,k]
# 			}
# 			teffect[m,k]=teffect[m,k]/[1 - rcor]
# 		}
# 	}
# 
# 	for (i in seq_len(length(varselect)) {
# 		for (j in seq_len(length(varselect)) {
# 			m <- varselect[i]
# 			k <- varselect[j]
# 			effect[m,k]=100*[effect[m,k]-teffect[m,k]]
# 		}
# 	}
# 	end



# # -------------------------------------------------------
# 	subroutine compprot[loz,nrnn,ns,nvarx,its,  	oobpred,wc,nclass,x,varselect,temp,cat,maxcat,  	jpur,inear,nprot,protlow,prothigh,prot,protfreq,  	protvlow,protvhigh,protv,popclass,npend,freq,v5,v95]
# # 
# 	integer nrnn,ns,nvarx,nclass,length(varselect),nprot,maxcat
# 	integer loz[ns,nrnn],oobpred[ns],varselect[nvarx],its[ns],  	jpur[nrnn],inear[nrnn],npend[nclass],cat[nvarx]
# 	real wc[ns],prot[nvarx,nprot,nclass],  	protlow[nvarx,nprot,nclass],prothigh[nvarx,nprot,nclass], 	protfreq[nvarx,nprot,nclass,maxcat],  	protvlow[nvarx,nprot,nclass],protvhigh[nvarx,nprot,nclass], 	x[nvarx,ns],temp[nrnn],protv[nvarx,nprot,nclass],  	popclass[nprot,nclass],freq[maxcat],v5[nvarx],v95[nvarx]
# 	integer ii,i,k,n,jp,npu,mm,m,nclose,jj,ll,jmax,nn
# 	real fmax,dt
# 
# 
# 	for (jp in seq_len(nclass) {
# 		its <- integer(ns)
# # 	we try to find nprot prototypes for this class:
# 		npend[jp]=nprot
# 		for (i in seq_len(nprot) {
# 			wc <- numeric(ns)
# 			for (n in seq_len(ns) {
# 				if (its[n] == 0) {
# # 				wc[n] is the number of unseen neighbors 
# # 				of example n that are predicted to be 
# # 				in class jp
# # 				loz[n,1],...,loz[n,nrnn] point to 
# # 				the nrnn nearest neighbors of 
# # 				example n 
# 					for (k in seq_len(nrnn) {
# 						nn <- loz[n,k]
# 						if (its[nn] == 0) {
# 							ii <- oobpred[nn]
# 							if[ii == jp] wc[n]=wc[n]+1
# 						}
# 					}
# 				}
# 			}
# # 		find the unseen example with the largest number 
# # 		of unseen predicted - class - jp neighbors
# 			nclose <- 0
# 			npu <- 0
# 			for (n in seq_len(ns) {
# 				if (wc[n] >= nclose && its[n] == 0) {
# 					npu <- n
# 					nclose <- wc[n]
# 				}
# 			}
# # 		if nclose <- 0,no example has any unseen predicted - class - jp neighbors
# # 		can't find another prototype for this class - reduce npend by 1 and
# # 		start finding prototypes for the next class
# 			if (nclose == 0) {
# 				npend[jp]=i - 1
# 				goto 93
# 			}
# # 		example npu has the largest number 
# # 		of unseen predicted - class - jp neighbors
# # 		put these neighbors in a list of length nclose
# 			ii <- 0	
# 			for (k in seq_len(nrnn) {
# 				nn <- loz[npu,k]
# 				if (its[nn] == 0 && oobpred[nn] == jp) {
# 					ii <- ii + 1
# 					inear[ii]=nn	
# 				}
# 			}
# # 		popclass is a measure of the size of the cluster around 
# # 		this prototype
# 			popclass[i,jp]=nclose
# 			for (mm in seq_len(length(varselect)) {
# # 			m is the index of the mmth variable
# 				m <- varselect[mm]	
# 				if (cat[m] == 1) {
# 					dt <- v95[m]-v5[m]
# 					for (ii in seq_len(nclose) {
# # 					put the value of the mmth variable into the list
# 						temp[ii]=x[inear[ii], m]
# 					} # ii
# # 				sort the list
# 					quicksort.result <- quicksort(temp,jpur,1,nclose,nclose)
# 					ii <- nclose / 4
# 					if[ii == 0] ii <- 1
# # 				find the 25th percentile
# 					protvlow[m,i,jp]=temp[ii]
# 					protlow[m,i,jp]=[temp[ii]-v5[m]]/dt
# 					ii <- nclose / 4
# 					ii=[3 * nclose]/4
# 					if[ii == 0] ii <- 1
# # 				find the 75th percentile
# 					protvhigh[m,i,jp]=temp[ii]
# 					prothigh[m,i,jp]=[temp[ii]-v5[m]]/dt
# 					ii <- nclose / 2
# 					if[ii == 0] ii <- 1
# # 				find the median
# 					protv[m,i,jp]=temp[ii]
# 					prot[m,i,jp]=[temp[ii]-v5[m]]/dt
# 				}	
# 				if (cat[m] >= 2) {
# # 				for categorical variables,choose the most frequent class
# 					freq <- numeric(maxcat)	
# 					for (k in seq_len(nclose) {
# 						jj <- nint[x[loz[npu,k], m]]	
# 						freq[jj]=freq[jj]+1
# 					}	
# 					jmax <- 1
# 					fmax <- freq[1]
# 					for (ll in 2:cat] {[m)
# 						if (freq[ll] > fmax) {
# 							jmax <- ll
# 							fmax <- freq[ll]
# 						}
# 					}	
# 					protv[m,i,jp]=jmax
# 					protvlow[m,i,jp]=jmax
# 					protvhigh[m,i,jp]=jmax
# 					for (ll in seq_len(cat] {[m)
# 						protfreq[m,i,jp,ll]=freq[ll]
# 					}
# 				}
# 			} # m
# # 		record that npu and it's neighbors have been 'seen'
# 			its[npu]=1
# 			for (k in seq_len(nclose) {
# 				nn <- loz[npu,k]
# 				its[nn]=1
# 			}
# 		} # nprot
# 93		continue
# 	} # jp
# 
# 	end



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
# 	subroutine locateout[y,tout,outtr,ncp,isort,devout, 	near,nexamples,nclass,rmedout]
# # 
# 	real outtr[near],tout[near],devout[nclass],rmedout[nclass]
#      	integer y[nexamples],isort[nexamples],ncp[near],near,nexamples, 	nclass
# 	real rmed, dev
# 	integer jp,nt,n,i
#      
# 	for (jp in seq_len(nclass) {
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



# -------------------------------------------------------
# R version. Used only when the R version of findbestsplit is used.
# unpack <- function(l,npack) {
#     icat <- integer(l)
#     
#     n <- npack
#     icat[1L] <- n %% 2L
#     for (k in 2L:l) {
#         n <- (n - icat[k - 1L]) / 2L
#         icat[k] <- n %% 2L
#     }
#     
#     as.logical(icat)
# }



# -------------------------------------------------------
# erfcc <- function(x) {
# 	z <- abs(x) / 1.41421356
# 	t <- 1 / (1 + 0.5 * z)
# 	erfcc <- t * exp(-z * z - 1.26551223 + t*(1.00002368 + t*(.37409196 + t*
#     (.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + t*
#     (1.48851587 + t*(-.82215223 + t*.17087277)))))))))
# 	erfcc <- erfcc / 2
# 	if (x < 0) {
#     erfcc <- 2 - erfcc
# 	}
# 	return(erfcc)
# }
