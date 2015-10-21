setGeneric("grow", function(forest, ...) standardGeneric("grow"))



setMethod("grow", signature(forest="bigcforest"), function(
    forest,
    x=NULL,
    ntrees=50L,
    printerrfreq=10L,
    printclserr=TRUE,
    reuse.cache=FALSE,
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
    
    if (trace >= 1L) message("Checking arguments in grow.")
    
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
            stop("Cache was not used to build this random forest. Cannot ",
                 "reuse cache.")
        }
        if (!file.exists(forest@cachepath)) {
            stop('Cache path "', forest@cachepath,
                 '" does not exist. Cannot reuse cache.')
        }
        
        if (any(!forest@factorvars)) {
            # There are numeric variables. Check big.matrix asave.
            if (!file.exists(paste0(forest@cachepath, "/", "asave"))) {
                stop('File "', paste0(forest@cachepath, "/", "asave"),
                     '" does not exist. Cannot reuse cache.')
            }
            if (!file.exists(paste0(forest@cachepath, "/", "asave.desc"))) {
                stop('File "', paste0(forest@cachepath, "/", "asave.desc"),
                     '" does not exist. Cannot reuse cache.')
            }
            asave <- attach.resource("asave.desc", path=forest@cachepath)
            if (nrow(asave) != forest@nexamples ||
                    ncol(asave) != sum(!forest@factorvars)) {
                stop("Dimensions of big.matrix asave do not match the number ",
                     "of rows of training data and the number of numeric ",
                     "variables. Cannot reuse cache.")
            }
        } else {
            # No numeric variables. Set up dummy big.matrix.
            asave <- big.matrix(1L, 1L, type="integer")
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
    
    # Check ntrees.
    if (!is.numeric(ntrees) ||
            abs(ntrees - round(ntrees)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument ntrees must be an integer.")
    }
    ntrees <- as.integer(round(ntrees))
    if (ntrees < 1L) {
        stop("Argument ntrees must be at least 1.")
    }
    
    # Check printerrfreq.
    if (!is.numeric(printerrfreq) ||
            abs(printerrfreq - round(printerrfreq)) >=
            .Machine$double.eps ^ 0.5) {
        stop ("Argument printerrfreq must be an integer.")
    }
    printerrfreq <- as.integer(round(printerrfreq))
    if (printerrfreq < 1L) {
        stop("Argument printerrfreq cannot be less than 1.")
    }
    
    # Check printclserr.
    if (!is.logical(printclserr)) {
        stop ("Argument printclserr must be a logical.")
    }
    
    

    # Initialize for run -------------------------------------------------------
    
    y <- forest@y
    
    # Convert x to big.matrix, as C functions in bigrf only support this.
    if (class(x) != "big.matrix") {
        if (trace >= 1L) message("Converting x into a big.matrix.")
        x <- makex(x, "x", forest@cachepath)
    }
    
    oldntrees <- forest@ntrees
    length(forest@trainerr) <- oldntrees + ntrees
    forest@trainclserr <- rbind(forest@trainclserr,
                                matrix(0, ntrees, length(levels(forest@y))))
    
    # Set up asave big.matrix.
    if (!reuse.cache) {
        if (any(!forest@factorvars)) {
            if (trace >= 1L) message("Setting up asave big.matrix.")
            asave <- makea(forest, x)
        } else {
            # No numeric variables. Set up dummy big.matrix.
            asave <- big.matrix(1L, 1L, type="integer")
        }
    }
    

    
    # Begin main loop ----------------------------------------------------------

    forest <- foreach(treenum=(forest@ntrees + 1):(forest@ntrees + ntrees),
                      .combine=combine.treeresults, .init=forest,
                      .inorder=FALSE, .verbose=FALSE) %dopar% {
        if (trace >= 1L) message("Growing tree ", treenum, " of ",
                                 oldntrees + ntrees, ".")
        
        # Take a bootstrap sample ----------------------------------------------
        
        if (trace >= 2L) message("Tree ", treenum, ": Taking bootstrap sample.")
        
        insamp <- integer(forest@nexamples)
        inweight <- numeric(forest@nexamples)
        
        k <- as.integer(runif(forest@nexamples, 0, forest@nexamples) + 1)
        t <- table(k)
        insamp[as.integer(names(t))] <- as.integer(t)
        inweight[as.integer(names(t))] <-
          as.integer(t) * forest@yclasswts[as.integer(y[as.integer(names(t))])]
        rm(k, t)
        
        # Set up a and a.out big.matrix's for caching --------------------------
        
        if (trace >= 2L) message("Tree ", treenum,
                                 ": Setting up a and a.out big.matrix's.")
        if (any(!forest@factorvars)) {
            if (is.null(forest@cachepath)) {
                a <- big.matrix(sum(insamp > 0L), sum(!forest@factorvars),
                                type="integer")
                a.out <- big.matrix(sum(insamp == 0L), sum(!forest@factorvars),
                                    type="integer")
            } else {
                a <- big.matrix(sum(insamp > 0L), sum(!forest@factorvars),
                                type="integer",
                                backingfile=paste0("a-", treenum),
                                descriptorfile=paste0("a-", treenum, ".desc"),
                                backingpath=forest@cachepath)
                a.out <- big.matrix(sum(insamp == 0L), sum(!forest@factorvars),
                                    type="integer",
                                    backingfile=paste0("a.out-", treenum),
                                    descriptorfile=paste0("a.out-", treenum,
                                                          ".desc"),
                                    backingpath=forest@cachepath)
                cachefiles <- paste0(forest@cachepath, "/",
                                     c("a-", "a-", "a.out-", "a.out-"), treenum,
                                     c("", ".desc", "", ".desc"))
                on.exit(file.remove(cachefiles))
            }
            .Call("modaC", asave@address, a@address, as.integer(insamp),
                  PACKAGE="bigrf")
            .Call("modaC", asave@address, a.out@address,
                  as.integer(insamp == 0L), PACKAGE="bigrf")
        } else {
            # No numeric variables. Set up dummy big.matrix's.
            a <- big.matrix(1L, 1L, type="integer")
            a.out <- big.matrix(1L, 1L, type="integer")
        }
        
        # Grow tree ------------------------------------------------------------
        
        if (trace >= 2L) message("Tree ", treenum, ": Growing tree.")
        
        xtype <- as.integer(bigmemory:::getCType(x))
        tree <- .Call("growtreeC", x@address, xtype, a@address, a.out@address,
                      forest, insamp, inweight, treenum, trace)
        
        # Clean up to free memory.
        rm(insamp, inweight, a, a.out, xtype)
        gc()
        
        list(treenum=treenum,
             oldntrees=oldntrees,
             ntrees=ntrees,
             tree=tree,
             printerrfreq=printerrfreq,
             printclserr=printclserr)
    }
  
    
    
    # Normalize votes ----------------------------------------------------------
    
    # if (trace >= 1L) message("Normalizing votes.")
    # w <- which(forest@oobtimes > 0)
    # for (c in seq_along(levels(forest@y))) {
    #   forest@oobvotes[w, c] <- forest@oobvotes[w, c] / forest@oobtimes[w]
    #   # if (ntest > 0L) {
    #   #   for (n in seq_len(ntest0)) {
    #   #     qts[c, n] <- qts[c, n] / ntrees
    #   #   }
    #   # }
    # }
    # rm(w)
    
    
    
    # Calculate confusion matrix -----------------------------------------------
    pred <- forest@oobpred
    pred[pred == 0L] <- length(levels(forest@y)) + 1L
    pred <- factor(pred, levels=seq_len(length(levels(forest@y)) + 1),
                   labels=c(levels(forest@y), "Never out-of-bag"))
    forest@trainconfusion <- table(y, pred, dnn=c("Actual", "Predicted"))
    rm(pred)
    
    
    
    # Print results ------------------------------------------------------------
    if (trace >= 1L) {
        cat("\n")
        summary(forest)
    }
    
    return(forest)
})
