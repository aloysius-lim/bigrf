setGeneric("proximities", function(forest, ...) standardGeneric("proximities"))



setMethod("proximities", signature(forest="bigcforest"), function(
    forest,
    y,
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
    if (trace < 0L || trace > 1L) {
        stop("Argument trace must be 0 or 1.")
    }
    
    if (trace >= 1L) message("Checking arguments.")
    
    # Check forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }

    # Check y.
    if (is.integer(y)) {
        if (min(y) < 1L) {
            stop("Elements in argument y must not be less than 1. The class ",
                 "labels coded in y should start with 1.")
        }
        y <- factor(y, seq_len(max(y)))
    } else if (!is.factor(y)) {
        stop("Argument y must be a factor or integer vector.")
    }
    if (length(y) != nrow(x)) {
        stop("Argument y must have as many elements as there are rows in x.")
    }
    if (!identical(forest@ylevels, levels(y)) ||
            !identical(forest@ynclass, length(levels(y))) ||
            !identical(forest@ytable, table(y, deparse.level=0))) {
        stop("Argument y is different than that used for building the random ",
             "forest.")
    }
    y <- as.integer(y)
    
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
    if (trace >= 1L) message("Creating big.matrix(s) for proximities.")
    if (nnearest < forest@nexamples) {
        if (is.null(cachepath)) {
            examples <- big.matrix(forest@nexamples, nnearest,
                                   type="double")
        } else {
            examples <- big.matrix(forest@nexamples, nnearest,
                                   type="double",
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
        
    
    
    # Compute proximities ------------------------------------------------------
    
    for (i in seq_len(forest@nexamples)) {
        if (trace >= 1L) message("Computing proximities for example ", i, ".")
        
        p <- numeric(forest@nexamples)
        for (t in seq_along(forest)) {
            tree <- forest[[t]]
            
            # Example numbers of all examples in the same node as example i.
            insame <- which(tree@trainprednode == tree@trainprednode[i])
            p[insame] <- p[insame] + 1L
            
            # Computation in original Fortran code, which does not make sense.
            # if (tree@insamp[i] > 0L) {
            #     w <- tree@insamp[insame] == 0L
            #     p[w] <- p[w] + wtx[i] /
            #         tree@termincount[tree@trainprednode[i]]
            # } else {
            #     w <- tree@insamp[insame] > 0L
            #     p[w] <- p[w] + wtx[w] / tree@termincount[tree@trainprednode[w]]
            # }
        }
        
        # if(noutlier.eq.1) then
        # 	rsq=0
        # 	do k=1,near
        # 	if(p(k).gt.0.and.cl(k).eq.cl(n)) rsq=rsq+p(k)*p(k)
        # 	enddo
        # 	if(rsq.eq.0) rsq=1
        # 	outtr(n)=near/rsq
        # endif
	
        if (nnearest == forest@nexamples) {
            prox[i, ] <- p / forest@ntrees
        } else {
            examples[i, ] <- order(p, decreasing=TRUE)[seq_len(nnearest)]
            prox[i, ] <- p[examples[i, ]] / forest@ntrees
        }
    }
    
    return(proximities)
})
