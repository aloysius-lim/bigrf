setGeneric("prototypes",
           function(forest, prox, ...) standardGeneric("prototypes"))



setMethod("prototypes", signature(forest="bigcforest", prox="bigrfprox"),
          function(
    forest,
    prox,
    nprot=1L,
    x=NULL,
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
    
    # Check prox.
    if (!class(prox) == "bigrfprox") {
        stop("Argument prox must be an object of class bigrfprox")
    }
    
    # Check nprot.
    if (!is.numeric(nprot) ||
            abs(nprot - round(nprot)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument nprot must be an integer.")
    }
    nprot <- as.integer(round(nprot))
    if (nprot < 1L) {
        stop("Argument trace must be at least 1.")
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
    
    
    
    # Initialize ---------------------------------------------------------------
    
    # Convert x to big.matrix, as C functions in bigrf only support this.
    if (class(x) != "big.matrix") {
        if (trace >= 1L) message("Converting x into a big.matrix.")
        x <- makex(x, "x", forest@cachepath)
    }
    
    # Parameters
    nvar <- length(forest@varselect)
    nnearest <- ncol(prox)
    ynclass <- length(levels(forest@y))
    
    # Quantiles of numeric variables. First row represents 5% quantile, second
    # row represents 95% quantile.
    numvarq <- matrix(numeric(), 2, nvar)
    numvarq[, !forest@factorvars] <-
        sapply(which(!forest@factorvars),
               function(v)
                   quantile(x[, forest@varselect[v]], probs=c(0.05, 0.95),
                            type=8))
    
    # Output variables.
    nprotfound <- integer(ynclass)
    names(nprotfound) <- levels(forest@y)
    
    clustersize <- matrix(integer(), ynclass, nprot,
                          dimnames=list(Class=levels(forest@y), Prototype=NULL))
    
    dimnames=list(Class=levels(forest@y),
                  Prototype=NULL,
                  Variable=names(forest@varselect),
                  Value=c("1st quartile", "median", "2nd quartile"))
    prot <- array(numeric(), dim=c(ynclass, nprot, nvar, 3L), dimnames)
    
    prot.std <- array(numeric(), dim=c(ynclass, nprot, nvar, 3L), dimnames)
    
    levelsfreq <- list()
    length(levelsfreq) <- nvar
    names(levelsfreq) <- names(forest@varselect)
    dimnames=list(Class=levels(forest@y),
                  Prototype=NULL,
                  Levels=NULL)
    for (i in which(forest@factorvars)) {
        levelsfreq[[i]] <- array(0L,
                                 c(ynclass, nprot, forest@varnlevels[i]),
                                 dimnames)
    }
    
    
    
    # Compute prototypes -------------------------------------------------------
    
    if (trace >= 1L) message("Computing prototypes.")
    
    for (c in seq_len(ynclass)) {
        
        # Compute nprot prototypes for this class.
        seen <- logical(forest@nexamples)
        
        for (p in seq_len(nprot)) {
            unseen <- which(!seen)
            
            # wc[n] is the number of unseen neighbors of case n that are 
            # predicted to be in class c.
            wc <- integer(forest@nexamples)
            wc[unseen] <- sapply(unseen, function(n) {
                if (nnearest < forest@nexamples) {
                    neighbors <- prox@examples[n, ] 
                    return(sum(!seen[neighbors] &
                                   forest@oobpred[neighbors] == c))
                } else {
                    return(sum(!seen & forest@oobpred == c))
                }
            })
            
            # If wc contains all zeros, no case has any unseen predicted-class-c
            # neighbors. Can't find another prototype for this class. Start
            # finding prototypes for the next class.
            if (all(wc == 0L)) {
                break
            }
            
            # Find the unseen case with the largest number of unseen
            # predicted-class-c neighbors, and put them in inear.
            npu <- which.max(wc)
            if (nnearest < forest@nexamples) {
                neighbors <- prox@examples[npu, ]
                inear <- neighbors[!seen[neighbors] &
                                       forest@oobpred[neighbors] == c]
            } else {
                inear <- which(!seen & forest@oobpred == c)
            }
            
            # clustersize is a measure of the size of the cluster around this
            # prototype
            clustersize[c, p] <- length(inear)
            for (var in seq_along(forest@varselect)) {
				if (!forest@factorvars[var]) {
                    # Find 25th, 50th and 75th percentiles.
                    prot[c, p, var, ] <-
                        quantile(x[inear, forest@varselect[var]],
                                 probs=c(0.25, 0.50, 0.75), type=8)
                    prot.std[c, p, var, ] <-
                        (prot[c, p, var, ] - numvarq[1, var]) /
                        (numvarq[2, var] - numvarq[1, var])
				} else {
                    # Choose the most frequent class.
                    t <- table(x[inear, forest@varselect[var]])
                    jmax <- as.integer(names(t)[which.max(t)])
                    prot[c, p, var, ] <- rep.int(jmax, 3L)
                    prot.std[c, p, var, ] <- prot[c, p, var, ] /
                        forest@varnlevels[var]
                    levelsfreq[[var]][c, p, as.integer(names(t))] <-
                        as.integer(t)
				}
            }
            
            # Record that npu and it's neighbors have been 'seen'
			seen[npu] <- TRUE
            seen[inear] <- TRUE
            nprotfound[c] <- p
		}
    }
    
    # Return.
    return(list(nprotfound=nprotfound,
                clustersize=clustersize,
                prot=prot,
                prot.std=prot.std,
                levelsfreq=levelsfreq))
})
