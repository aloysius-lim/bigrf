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
    }
    for (var in which(factorvars)) {
        asave[, var] <- as.integer(x[, varselect[var]])
    }
    
    return()
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
buildtree <- function(x, asave, a, a.out, forest, insamp, inweight, treenum,
                      trace) {
    xtype <- as.integer(.Call("CGetType", x@address, PACKAGE="bigmemory"))
    return(.Call("buildtreeC", x@address, xtype, asave@address, a@address,
                 a.out@address, forest, insamp, inweight, treenum, trace))
}



# ------------------------------------------------------------------------------
# Combine results of tree builds. To be used only as a .combine function in
# foreach().
combine.treeresults <- function(forest, newtree) {
    y <- forest@y
    treenum <- forest@ntrees + 1L
    oldntrees <- newtree$oldntrees
    ntrees <- newtree$ntrees
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
