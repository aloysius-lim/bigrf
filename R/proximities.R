setGeneric("proximities", function(forest, ...) standardGeneric("proximities"))



setMethod("proximities", signature(forest="bigcforest"), function(
    forest,
    nnearest=forest@nexamples,
    cachepath=tempdir(),
    trace=0L) {
    
    # Check arguments ----------------------------------------------------------
    
    # Check trace.
    if (!is.numeric(trace) ||
            abs(trace - round(trace)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument trace must be an integer.")
    }
    trace <- as.integer(round(trace))
    if (trace < 0L || trace > 2L) {
        stop("Argument trace must be 0, 1 or 2.")
    }
    
    if (trace >= 1L) message("Checking arguments.")
    
    # Check forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }

    # Check nnearest.
    if (!is.numeric(nnearest) ||
            abs(nnearest - round(nnearest)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument nnearest must be an integer.")
    }
    nnearest <- as.integer(round(nnearest))
    if (nnearest < 1L || nnearest > forest@nexamples) {
        stop("Argument trace must be at least 1 and no greater than ",
             "forest@nexamples.")
    }

    # Check cachepath.
    if (!(is.null(cachepath) || is.character(cachepath))) {
        stop("Argument cachepath must be a character string, or NULL.")
    }
    if (!is.null(cachepath)) {
        if (!file.exists(cachepath)) {
            if (!dir.create(cachepath)) {
                stop("Cannot create directory ", cachepath, ".")
            }
        }
    }

    
    
    # Initialize ---------------------------------------------------------------
    
    if (trace >= 1L) message("Initializing.")
    
    # Create big.matrix for proximities.
    if (is.null(cachepath)) {
        prox <- big.matrix(forest@nexamples, nnearest, type="double")
    } else {
        prox <- big.matrix(forest@nexamples, nnearest, type="double",
                           backingfile="prox",
                           descriptorfile="prox.desc",
                           backingpath=cachepath)
    }
    
    # Create big.matrix for indices of examples in proximity matrix (only if
    # nnearest < forest@nexamples).
    if (nnearest < forest@nexamples) {
        if (is.null(cachepath)) {
            examples <- big.matrix(forest@nexamples, nnearest,
                                   type="integer")
        } else {
            examples <- big.matrix(forest@nexamples, nnearest,
                                   type="integer",
                                   backingfile="prox.examples",
                                   descriptorfile="prox.examples.desc",
                                   backingpath=cachepath)
        }
    } else {
        examples <- NULL
    }
    
    # Create promixities object.
    proximities <- new("bigrfprox",
                       address=prox@address,
                       examples=examples,
                       cachepath=cachepath)
        
    
    
    # Pre-compute which examples fall in each node -----------------------------
    
    if (trace >= 1L) message("Pre-computing examples in each node of each ",
                             "tree.")
    
    node.examples <- list()
    length(node.examples) <- forest@ntrees
    
    node.examples <- foreach(t=seq_along(forest)) %dopar% {
        if (trace >= 2L) message("Processing tree ", t, ".")
        
        tree <- forest[[t]]
        nex <- list()
        length(nex) <- tree@nnodes
        
        for (i in which(tree@treemap[, 1] == 0L)) {
            nex[[i]] <- which(tree@trainprednode == i)
        }
        
        nex
    }
    
    
    
    # Compute proximities ------------------------------------------------------
    
    if (trace >= 1L) message("Computing proximities.")
    # Process in batchsize batches of batchsize examples.
    batchsize <- as.integer(round(sqrt(forest@nexamples)))
    batches <- data.frame(start=(1:batchsize) * batchsize - (batchsize - 1L),
                          end=(1:batchsize) * batchsize)
    batches <- batches[batches$start <= forest@nexamples, ]
    batches$end[nrow(batches)] <- forest@nexamples
    
    foreach(b=seq_len(nrow(batches)), .inorder=FALSE) %dopar% {
        if (trace >= 2L) message("Computing proximities for examples ",
                                 batches$start[b], " to ", batches$end[b], ".")
        
        for (i in batches$start[b]:batches$end[b]) {
            p <- numeric(forest@nexamples)
            for (t in seq_along(forest)) {
                tree <- forest[[t]]
                # Example numbers of all examples in the same node as example i.
                w <- node.examples[[t]][[tree@trainprednode[i]]]
                p[w] <- p[w] + 1L
            }
            
            if (nnearest == forest@nexamples) {
                prox[i, ] <- p / forest@ntrees
            } else {
                examples[i, ] <- order(p, decreasing=TRUE)[seq_len(nnearest)]
                prox[i, ] <- p[examples[i, ]] / forest@ntrees
            }
        }
        
        NULL
    }
    
    # Return.
    return(proximities)
})
