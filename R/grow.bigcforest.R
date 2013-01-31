grow.bigcforest <- function(forest,
                            x,
                            y=NULL,
                            ntrees=500L,
                            reuse.cache=FALSE,
                            trace=FALSE) {
    
    # Check arguments ----------------------------------------------------------
    
    # Check trace.
    if (!is.logical(trace)) {
        stop ("Argument trace must be a logical.")
    }
    
    if (trace) message("Checking arguments.")
    
    # Check forest.
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix.")
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
            if (!file.exists(paste0(forest@cachepath, "/", "x"))) {
                stop('File "', paste0(forest@cachepath, "/", "x"),
                     '" does not exist. Cannot reuse cache.')
            }
            if (!file.exists(paste0(forest@cachepath, "/", "x.desc"))) {
                stop('File "', paste0(forest@cachepath, "/", "x.desc"),
                     '" does not exist. Cannot reuse cache.')
            }
        }
    }
    
    

    # Convert x to big.matrix, as C functions only support this at the moment.
    if (class(x) != "big.matrix") {
        x <- as.big.matrix(x)
    }
    
    # Add synthetic class for unsupervised learning ----------------------------
    
    if (!forest@supervised) {
        if (reuse.cache) {
            x <- attach.resource("x.desc", path=forest@cachepath)
        } else {
            if (trace) message("Creating a synthetic class for unsupervised ",
                               "learning.")
            x.old <- x
            x <- big.matrix(forest@nsample, ncol(x), type="double",
                            backingfile="x.unsupervised",
                            descriptorfile="x.unsupervised.desc",
                            backingpath=forest@cachepath)
            for (m in seq_len(ncol(x))) {
                k <- as.integer(runif(forest@nsample - nrow(x), 0, nrow(x)) + 1)
                x[seq_len(nrow(x)), m] <- x.old[seq_len(nrow(x)), m]
                x[(nrow(x) + 1):forest@nsample, m] <- x.old[k, m]
            }
            rm(x.old, m, k)
        }
        
        y <- c(rep.int(1L, nrow(x)), rep.int(2L, forest@nsample - nrow(x)))
        forest@nclass <- 2L
    }
    
    
    
    # Initialize for run -------------------------------------------------------
    
    oldntrees <- forest@ntrees
    
    # Set up asave big.matrix.
    if (reuse.cache) {
        if (trace) message("Loading asave big.matrix.")
        asave <- attach.resource("asave.desc", path=forest@cachepath)
    } else {
        if (trace) message("Setting up asave big.matrix.")
        asave <- big.matrix(forest@nsample, length(forest@varselect),
                            type="integer",
                            backingfile="asave", descriptorfile="asave.desc",
                            backingpath=forest@cachepath)
        makea(x, asave, forest@factors, forest@varselect)
    }
    
    
    
    # Begin main loop ----------------------------------------------------------

    forest <- foreach(treenum=(forest@ntrees + 1):(forest@ntrees + ntrees),
                      .combine=combine.treeresults, .init=forest,
                      .inorder=FALSE, .verbose=FALSE) %dopar% {
        if (trace) message("Building tree ", treenum, " of ",
                           oldntrees + ntrees, ".")
        
        # Take a bootstrap sample ----------------------------------------------
        
        if (trace) message("Tree ", treenum, ": Taking bootstrap sample.")
        
        insamp <- integer(forest@nsample)
        inweight <- numeric(forest@nsample)
        
        k <- as.integer(runif(forest@nsample, 0, forest@nsample) + 1)
        t <- table(k)
        insamp[as.integer(names(t))] <- as.integer(t)
        inweight[as.integer(names(t))] <-
          as.integer(t) * forest@classweights[y[as.integer(names(t))]]
        rm(k, t)
        
        # Set up a and a.out big.matrix's for caching --------------------------
        
        if (trace) message("Tree ", treenum,
                           ": Setting up a and a.out big.matrix's.")
        cachefiles <- paste0(forest@cachepath, "/",
                             c("a-", "a-", "a.out-", "a.out-"), treenum,
                             c("", ".desc", "", ".desc"))
        if (any(!forest@factors)) {
            a <- big.matrix(sum(insamp > 0L), sum(!forest@factors),
                            type="integer",
                            backingfile=paste0("a-", treenum),
                            descriptorfile=paste0("a-", treenum, ".desc"),
                            backingpath=forest@cachepath)
            a.out <- big.matrix(sum(insamp == 0L), sum(!forest@factors),
                                type="integer",
                                backingfile=paste0("a.out-", treenum),
                                descriptorfile=paste0("a.out-", treenum, ".desc"),
                                backingpath=forest@cachepath)
            moda(asave, a, forest@factors, insamp)
            moda(asave, a.out, forest@factors, insamp == 0L)
        } else {
            a <- big.matrix(1L, 1L,
                            type="integer",
                            backingfile=paste0("a-", treenum),
                            descriptorfile=paste0("a-", treenum, ".desc"),
                            backingpath=forest@cachepath)
            a.out <- big.matrix(1L, 1L,
                                type="integer",
                                backingfile=paste0("a.out-", treenum),
                                descriptorfile=paste0("a.out-", treenum, ".desc"),
                                backingpath=forest@cachepath)
        }
        on.exit(file.remove(cachefiles))
        
        # Build tree -----------------------------------------------------------
        
        if (trace) message("Tree ", treenum, ": Building tree.")
        tree <- buildtree(x, y, asave, a, a.out, forest, insamp, inweight,
                          treenum, trace)
        list(treenum=treenum,
           oldntrees=oldntrees,
           ntrees=ntrees,
           y=y,
           insamp=insamp,
           tree=tree,
           trace=trace)
        
    }
  
    
    
    # Normalize votes ----------------------------------------------------------
    
    # if (trace) message("Normalizing votes.")
    # w <- which(forest@oobtimes > 0)
    # for (c in seq_len(forest@nclass)) {
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
    pred[pred == 0L] <- forest@nclass + 1L
    if (length(forest@ylevels)) {
        class(pred) <- "factor"
        levels(pred) <- c(forest@ylevels, "Never out-of-bag")
        class(y) <- "factor"
        levels(y) <- forest@ylevels
    }
    forest@trainconfusion <- table(y, pred, dnn=c("Actual", "Predicted"))
    rm(pred)
    
    
    
    # Print results ------------------------------------------------------------
    cat("\n")
    show(forest)
    

    
    # # -------------------------------------------------------
    # # SEND INFO ON TRAINING AND/OR TEST SET DATA TO FILE 
    #
    # if (idataout == 2 && ntest > 0) {
    # 	if (labelts == 1) {
    # 		for (n in seq_len(ntest0) {
    # 			write[7,'[3i5,1000f10.3]'] n,clts[n],jests[n],
    #       [qts[j,n],j=1,nclass],[xts[m,n],m=1,ncol(x)]
    # 		}
    # 	} else {
    # 		for (n in seq_len(ntest0) {
    # 			write[7,'[2i5,1000f10.3]'] n,jests[n],
    #       [qts[j,n],j=1,nclass],[xts[m,n],m=1,ncol(x)]
    # 		}
    # 	}
    # }
    
    # # -------------------------------------------------------
    # # SEND FILL TO FILE [ROUGH FILL ONLY]
    # # 
    # 	if (isavefill == 1 && missfill == 1) {
    # 		write[3,*] [fill[m],m=1,ncol(x)]
    # 	}
    
    return(forest)
}

setMethod("grow", signature(forest="bigcforest"), grow.bigcforest)
