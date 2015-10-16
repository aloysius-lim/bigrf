setGeneric("varimp", function(forest, ...) standardGeneric("varimp"))



setMethod("varimp", signature(forest="bigcforest"),  function(
    forest,
    x=NULL,
    impbyexample=FALSE,
    reuse.cache=FALSE,
    trace=0L) {
    
    # Check arguments ----------------------------------------------------------
    
    # Check trace.
    if (!is.numeric(trace) ||
            abs(trace - round(trace)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument trace must be an integer.")
    }
    trace <- as.integer(round(trace))
    if (trace < 0L || trace > 1L) {
        stop("Argument trace must be 0 or 1.")
    }
    
    if (trace >= 1L) message("Checking arguments.")
    
    # Check forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }
    
    # Check reuse.cache.
    if (!is.logical(reuse.cache)) {
        stop ("Argument reuse.cache must be a logical.")
    }
    if (reuse.cache) {
        if (is.null(forest@cachepath)) {
            stop("Cache was not used to build this random forest. Cannot reuse ",
                 "cache.")
        }
        if (!file.exists(forest@cachepath)) {
            stop('Cache path "', forest@cachepath,
                 '" does not exist. Cannot reuse cache.')
        }
        
        if (is.null(x)) {
            if (!file.exists(paste0(forest@cachepath, "/", "x"))) {
                stop('File "', paste0(forest@cachepath, "/", "x"),
                     '" does not exist. Cannot reuse cache.')
            }
            if (!file.exists(paste0(forest@cachepath, "/",
                                    "x.desc"))) {
                stop('File "',
                     paste0(forest@cachepath, "/", "x.desc"),
                     '" does not exist. Cannot reuse cache.')
            }
            x <- attach.resource("x.desc", path=forest@cachepath)
            if (nrow(x) != forest@nexamples) {
                stop("Big.matrix x in the cache path does not have the same ",
                     "number of rows as the training data used to build this ",
                     "random forest. Cannot reuse cache.")
            }
        }
    }
    
    # Check x.
    if (!reuse.cache) {
        if (is.null(x)) {
            stop("Argument x must be specified if reuse.cache is FALSE.")
        }
        if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
            stop("Argument x must be a big.matrix, matrix or data.frame.")
        }
        if (nrow(x) != forest@nexamples) {
            stop("Number of rows in argument x does not match number of rows ",
                 "in the training data used to build this forest.")
        }
    }
    
    # Check impbyexample.
    if (!is.logical(impbyexample)) {
        stop ("Argument impbyexample must be a logical.")
    }
    
    
    
    # Initialize ---------------------------------------------------------------
    
    y.int <- as.integer(forest@y)
    
    # Convert x to big.matrix, as C functions in bigrf only support this.
    if (class(x) != "big.matrix") {
        if (trace >= 1L) message("Converting x into a big.matrix.")
        x <- makex(x, "x", forest@cachepath)
    }
    
    nvar <- length(forest@varselect)
    ntest <- as.integer(nrow(x));
    xtype <- as.integer(bigmemory:::getCType(x))
    
    importance <- numeric(nvar)
    importance.sq <- numeric(nvar)
    names(importance) <- names(forest@varselect)
    if (impbyexample) {
        importance.ex.total <- numeric(forest@nexamples)
        importance.ex <- matrix(0, forest@nexamples, nvar,
                                dimnames=list(NULL, names(forest@varselect)))
    } else {
        importance.ex <- NULL
    }
    
    
    
    # Calculate variable importance --------------------------------------------
    
    for (treenum in seq_along(forest)) {
        if (trace >= 1L) message("Processing tree ", treenum, ".")
        
        tree <- forest[[treenum]]
        
        # Count correct oob classifications.
        w <- which(tree@insamp == 0L & tree@trainpredclass == y.int)
        correct <- sum(tree@nodewt[tree@trainprednode[w]])
        nout <- sum(tree@insamp == 0L)
        whichout <- which(tree@insamp == 0L)
        rm(w)
        
        # Variable importance for each example.
        if (impbyexample) {
            w <- which(tree@trainpredclass[whichout] == y.int[whichout])
            w <- whichout[w]
            importance.ex.total[w] <- importance.ex.total[w] +
                tree@nodewt[tree@trainprednode[w]] / nout
            rm(w)
        }
        
        # Which variables were split on in this tree.
        tbestvars <- sort(unique(tree@bestvar[tree@treemap[, 1] > 0L]))
        
        for (var in tbestvars) {
            # Predict class of out-of-bag samples for importance computation.
            whichpermout <- sample(whichout, length(whichout))
            oobpredict.result <- .Call("treepredictimpC", x@address, xtype,
                                       nout, whichout, whichpermout, var,
                                       forest, tree);
            
            w <- which(oobpredict.result$oobpredclass == y.int[whichout])
            impcorrect <- sum(tree@nodewt[oobpredict.result$oobprednode[w]])
            rm(w)
            
            if (impbyexample) {
                w <- which(oobpredict.result$oobpredclass == y.int[whichout])
                importance.ex[whichout[w], var] <-
                    importance.ex[whichout[w], var] +
                    tree@nodewt[oobpredict.result$oobprednode[w]] / nout
                rm(w)
            }
            
            importance[var] <- importance[var] + (correct - impcorrect) / nout
            importance.sq[var] <- importance.sq[var] +
                ((correct - impcorrect) / nout) ^ 2
        }
        if (impbyexample) {
            w <- which(tree@trainpredclass[whichout] == y.int[whichout])
            w <- whichout[w]
            weights <- tree@nodewt[tree@trainprednode[w]] / nout
            for (var in seq_len(nvar)[-tbestvars]) {
                importance.ex[w, var] <- importance.ex[w, var] + weights
            }
            rm(w)
        }
    }
    
    importance <- importance / forest@ntrees
    importance.sq <- importance.sq / forest@ntrees
    se <- sqrt((importance.sq - importance ^ 2) / forest@ntrees)
    zscore <- importance / se
    # Complementary error function
    significance <- 2 * pnorm(zscore * sqrt(2), lower.tail = FALSE)
    
    if (impbyexample) {
        for (var in seq_len(nvar)) {
            importance.ex[, var] <- 100 *
                (importance.ex.total - importance.ex[, var]) / forest@ntrees
        }
    }
    
    return(list(importance=importance, importance.ex=importance.ex,
                zscore=zscore, significance=significance))
})
          
