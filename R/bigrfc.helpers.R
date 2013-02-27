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
        if (class(x) == "data.frame" && class(x[[j]]) == "logical") {
            # Logical columns in data.frames need to be converted to integers,
            # with 2 for TRUE and 1 for FALSE.
            xcol <- as.integer(x[[j]]) + 1L
        } else if (xtype %in% c("integer", "char")) {
            xcol <- as.integer(x[, j])
        } else {
            xcol <- x[, j]
        }
        xnew[, j] <- xcol
    }
    options(old.opt)
    
    return(xnew)
}



# ------------------------------------------------------------------------------
# makea constructs the big.matrix asave with each column corresponding to a
# numerical variable in x. Each column stores the index number of the training
# examples, sorted in increasing order of the corresponding numerical variable.
# makea <- function(x, asave, factorvars, varselect) {
makea <- function(forest, x) {
    
    if (is.null(forest@cachepath)) {
        asave <- big.matrix(forest@nexamples, sum(!forest@factorvars),
                            type="integer")
    } else {
        asave <- big.matrix(forest@nexamples, sum(!forest@factorvars),
                            type="integer",
                            backingfile="asave",
                            descriptorfile="asave.desc",
                            backingpath=forest@cachepath)
    }
    
    w <- which(!forest@factorvars)
    for (i in seq_along(w)) {
        asave[, i] <- as.integer(morder(x, forest@varselect[w[i]]))
    }
    
    return(asave)
}



# ------------------------------------------------------------------------------
# Combine results of tree-growing. To be used only as a .combine function in
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
    
    for (c in seq_along(levels(y))) {
        # Out-of-bag examples with votes for this class.
        w <- which(tree@trainpredclass == c & tree@insamp == 0L)
        forest@oobvotes[w, c] <- forest@oobvotes[w, c] +
            tree@nodewt[tree@trainprednode[w]]
    }
    rm(c, w)
    
    # Get training set error estimates -----------------------------------------
    
    forest@oobpred[forest@oobtimes > 0L] <-
        max.col(forest@oobvotes[forest@oobtimes > 0L, ])
    
    for (c in seq_along(levels(y))) {
        forest@trainclserr[treenum, c] <- sum(as.integer(y) == c &
                                                  forest@oobpred != c)
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
                       width=max(nchar(levels(y)), 5), digits=3, nsmall=2),
                sep="  ")
        }
        cat("\n")
    }
    
    return(forest)
}



# ------------------------------------------------------------------------------
# Combine results of tree predictions. To be used only as a .combine function in
# foreach().
combine.treepredictresults <- function(prediction, treepredict.result) {
    y <- treepredict.result$y
    forest <- treepredict.result$forest
    tree <- treepredict.result$tree
    t <- treepredict.result$t
    printerrfreq <- treepredict.result$printerrfreq
    printclserr <- treepredict.result$printclserr
    
    # Compute votes.
    for (c in seq_along(levels(forest@y))) {
        w <- which(treepredict.result$testpredclass == c)
        prediction@testvotes[w, c] <- prediction@testvotes[w, c] +
            tree@nodewt[treepredict.result$testprednode[w]]
    }
    rm(c, w)
    prediction[seq_len(prediction@ntest)] <- max.col(prediction@testvotes)
    
    # If test set labels were given, compute test error.
    if (!is.null(y)) {
        prediction@testclserr <- integer(length(levels(forest@y)))
        names(prediction@testclserr) <- levels(forest@y)
        for (c in seq_along(levels(forest@y))) {
            prediction@testclserr[c] <-
                sum(as.integer(y) == c & prediction[] != c)
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
