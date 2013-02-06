setGeneric("grow", function(forest, ...) standardGeneric("grow"))



grow.bigcforest <- function(forest,
                            x=NULL,
                            y,
                            ntrees=500L,
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
            stop("Cache was not used to grow this random forest. Cannot reuse ",
                 "cache.")
        }
        if (is.null(forest@cachepath)) {
            stop("Cache was not used to grow this random forest. Cannot reuse ",
                 "cache.")
        }
        if (!file.exists(forest@cachepath)) {
            stop('Cache path "', forest@cachepath,
                 '" does not exist. Cannot reuse cache.')
        }
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
                ncol(asave) != length(forest@varselect)) {
            stop("Big.matrix asave in the cache path does not have the same ",
                 "dimensions as the training data used to grow this random ",
                 "forest. Cannot reuse cache.")
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
                     "number of rows as the training data used to grow this ",
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
                 "in the training data used to grow this forest.")
        }
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
    ytable <- table(y, deparse.level=0)
    y <- as.integer(y)
    
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
    
    # Convert x to big.matrix, as C functions in bigrf only support this.
    if (class(x) != "big.matrix") {
        if (trace >= 1L) message("Converting x into a big.matrix.")
        x <- makex(x, "x", forest@cachepath)
    }
    
    oldntrees <- forest@ntrees
    length(forest@trainerr) <- oldntrees + ntrees
    forest@trainclserr <- rbind(forest@trainclserr,
                                matrix(0, ntrees, forest@ynclass))
    
    # Set up asave big.matrix.
    if (!reuse.cache) {
        if (trace >= 1L) message("Setting up asave big.matrix.")
        if (is.null(forest@cachepath)) {
            asave <- big.matrix(forest@nexamples, length(forest@varselect),
                                type="integer")
        } else {
            asave <- big.matrix(forest@nexamples, length(forest@varselect),
                                type="integer",
                                backingfile="asave",
                                descriptorfile="asave.desc",
                                backingpath=forest@cachepath)
        }
        makea(x, asave, forest@factorvars, forest@varselect)
    }
    
    
    
    # Begin main loop ----------------------------------------------------------

    forest <- foreach(treenum=(forest@ntrees + 1):(forest@ntrees + ntrees),
                      .combine=combine.treeresults, .init=forest,
                      .inorder=FALSE, .verbose=FALSE) %dopar% {
        if (trace >= 1L) message("Building tree ", treenum, " of ",
                                 oldntrees + ntrees, ".")
        
        # Take a bootstrap sample ----------------------------------------------
        
        if (trace >= 2L) message("Tree ", treenum, ": Taking bootstrap sample.")
        
        insamp <- integer(forest@nexamples)
        inweight <- numeric(forest@nexamples)
        
        k <- as.integer(runif(forest@nexamples, 0, forest@nexamples) + 1)
        t <- table(k)
        insamp[as.integer(names(t))] <- as.integer(t)
        inweight[as.integer(names(t))] <-
          as.integer(t) * forest@yclasswts[y[as.integer(names(t))]]
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
            }
            moda(asave, a, forest@factorvars, insamp)
            moda(asave, a.out, forest@factorvars, insamp == 0L)
        } else {
            if (is.null(forest@cachepath)) {
                a <- big.matrix(1L, 1L, type="integer")
                a.out <- big.matrix(1L, 1L, type="integer")
            } else {
                a <- big.matrix(1L, 1L, type="integer",
                                backingfile=paste0("a-", treenum),
                                descriptorfile=paste0("a-", treenum, ".desc"),
                                backingpath=forest@cachepath)
                a.out <- big.matrix(1L, 1L, type="integer",
                                    backingfile=paste0("a.out-", treenum),
                                    descriptorfile=paste0("a.out-", treenum,
                                                          ".desc"),
                                    backingpath=forest@cachepath)
            }
        }
        if (!is.null(forest@cachepath)) {
            cachefiles <- paste0(forest@cachepath, "/",
                                 c("a-", "a-", "a.out-", "a.out-"), treenum,
                                 c("", ".desc", "", ".desc"))
            on.exit(file.remove(cachefiles))
        }
        
        # Build tree -----------------------------------------------------------
        
        if (trace >= 2L) message("Tree ", treenum, ": Building tree.")
        tree <- buildtree(x, y, asave, a, a.out, forest, insamp, inweight,
                          treenum, trace)
        list(treenum=treenum,
             oldntrees=oldntrees,
             ntrees=ntrees,
             y=y,
             tree=tree,
             printerrfreq=printerrfreq,
             printclserr=printclserr)
    }
  
    
    
    # Normalize votes ----------------------------------------------------------
    
    # if (trace >= 1L) message("Normalizing votes.")
    # w <- which(forest@oobtimes > 0)
    # for (c in seq_len(forest@ynclass)) {
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
    pred[pred == 0L] <- forest@ynclass + 1L
    pred <- factor(pred, levels=seq_len(forest@ynclass + 1),
                   labels=c(forest@ylevels, "Never out-of-bag"))
    y <- factor(y, levels=seq_len(forest@ynclass), labels=forest@ylevels)
    forest@trainconfusion <- table(y, pred, dnn=c("Actual", "Predicted"))
    rm(pred)
    
    
    
    # Print results ------------------------------------------------------------
    cat("\n")
    summary(forest)
    

    
    # # -------------------------------------------------------
    # # SEND FILL TO FILE [ROUGH FILL ONLY]
    # # 
    # 	if (isavefill == 1 && missfill == 1) {
    # 		write[3,*] [fill[m],m=1,ncol(x)]
    # 	}
    
    return(forest)
}

setMethod("grow", signature(forest="bigcforest"), grow.bigcforest)
