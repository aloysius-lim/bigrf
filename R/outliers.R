setGeneric("outliers", function(forest, ...) standardGeneric("outliers"))



setMethod("outliers", signature(forest="bigcforest"), function(
    forest,
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

    
    
    # Initialize ---------------------------------------------------------------
    
    y <- forest@y
    outscore <- numeric(forest@nexamples)
    
    
    
    # Compute proximities ------------------------------------------------------
    if (trace >= 1L) message("Computing raw outlier scores.")
    
    for (i in seq_len(forest@nexamples)) {
        
        p <- numeric(forest@nexamples)
        for (t in seq_along(forest)) {
            tree <- forest[[t]]
            
            # Example numbers of all examples in the same node as example i.
            w <- which(tree@trainprednode == tree@trainprednode[i])
            p[w] <- p[w] + 1L
        }
        
        w <- which(p > 0 & y == y[i])
        rsq <- sum(p[w] ^ 2)
        if(rsq == 0) {
            rsq <- 1
        }
        outscore[i] <- forest@nexamples / rsq
    }
    
    
    
    # Compute outlier scores ---------------------------------------------------
    if (trace) message("Computing outlier scores.")
    
    for (c in seq_along(levels(forest@y))) {
        w <- which(as.integer(y) == c)
        classout <- outscore[w]
        med <- median(classout)
        dev <- sum(abs(classout - med)) / length(w)
        outscore[w] <- (outscore[w] - med) / dev
    }
    
    # Return.
    return(outscore)
})
