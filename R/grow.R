setGeneric("grow", function(forest, ...) standardGeneric("grow"))



grow.bigcforest <- function(forest,
                            x,
                            y=NULL,
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
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix, matrix or data.frame.")
    }
    
    # Check y.
    if (forest@supervised) {
        if (is.null(y)) {
            stop("Argument y must be specified for supervised learning.")
        }
        if (is.factor(y)) {
            y <- as.integer(y)
        } else if (!is.integer(y)) {
            stop("Argument y must be a factor or integer vector of class ",
                 "labels.")
        }
        if (length(y) != nrow(x)) {
            stop("Argument y must have as many elements as there are rows in ",
                 "x.")
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
    
    # Check reuse.cache.
    if (!is.logical(reuse.cache)) {
        stop ("Argument reuse.cache must be a logical.")
    }
    if (reuse.cache) {
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
        if (!forest@supervised) {
            if (!file.exists(paste0(forest@cachepath, "/", "x.unsupervised"))) {
                stop('File "', paste0(forest@cachepath, "/", "x.unsupervised"),
                     '" does not exist. Cannot reuse cache.')
            }
            if (!file.exists(paste0(forest@cachepath, "/",
                                    "x.unsupervised.desc"))) {
                stop('File "',
                     paste0(forest@cachepath, "/", "x.unsupervised.desc"),
                     '" does not exist. Cannot reuse cache.')
            }
        }
    }
    
    

    # Convert x to big.matrix, as C functions only support this at the moment.
    if (class(x) != "big.matrix") {
        if (is.null(forest@cachepath)) {
            x <- as.big.matrix(x)
        } else {
            x <- as.big.matrix(x, backingfile="x", descriptorfile="x.desc",
                               backingpath=forest@cachepath)
        }
    }
    
    # Add synthetic class for unsupervised learning ----------------------------
    
    if (!forest@supervised) {
        if (reuse.cache) {
            x <- attach.resource("x.unsupervised.desc", path=forest@cachepath)
        } else {
            if (trace >= 1L) message("Creating a synthetic class for ",
                                     "unsupervised learning.")
            x.old <- x
            if (is.null(forest@cachepath)) {
                x <- big.matrix(forest@nexamples, ncol(x), type=typeof(x.old))
            } else {
                x <- big.matrix(forest@nexamples, ncol(x), type=typeof(x.old),
                                backingfile="x.unsupervised",
                                descriptorfile="x.unsupervised.desc",
                                backingpath=forest@cachepath)
            }
            .Call("synthesizeUnsupervisedC", x.old@address, x@address,
                  as.integer(.Call("CGetType", x.old@address,
                                   PACKAGE="bigmemory")))
            rm(x.old)
        }
        
        y <- c(rep.int(1L, nrow(x) / 2),
               rep.int(2L, forest@nexamples - nrow(x) / 2))
    }
    
    
    
    # Tabulate y ---------------------------------------------------------------
    class(y) <- "factor"
    levels(y) <- forest@ylevels
    forest@ytable <- table(y, deparse.level=0);
    y <- as.integer(y)
    
    

    # Initialize for run -------------------------------------------------------
    
    oldntrees <- forest@ntrees
    length(forest@trainerr) <- oldntrees + ntrees
    forest@trainclserr <- rbind(forest@trainclserr,
                                matrix(0, ntrees, forest@ynclass))
    
    # Set up asave big.matrix.
    if (reuse.cache) {
        if (trace >= 1L) message("Loading asave big.matrix.")
        asave <- attach.resource("asave.desc", path=forest@cachepath)
    } else {
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
    class(pred) <- "factor"
    levels(pred) <- c(forest@ylevels, "Never out-of-bag")
    class(y) <- "factor"
    levels(y) <- forest@ylevels
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
