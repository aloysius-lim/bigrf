setGeneric("varimp", function(forest, ...) standardGeneric("varimp"))



varimp.bigcforest <- function(forest, x=NULL, y=NULL,
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
    
    # Check impbyexample.
    if (!is.logical(impbyexample)) {
        stop ("Argument impbyexample must be a logical.")
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
    
    
    
    # Initialize ---------------------------------------------------------------
    
    nvar <- length(forest@varselect)
    ntest <- as.integer(nrow(x));
    xtype <- as.integer(.Call("CGetType", x@address, PACKAGE="bigmemory"))
    
    if (!is.null(y)) {
        class(y) <- "factor"
        levels(y) <- forest@ylevels
        ytable <- table(y, deparse.level=0)
        y <- as.integer(y)
    } else {
        ytable <- NULL
    }
    
    importance <- numeric(nvar)
    importance.sq <- numeric(nvar)
    if (impbyexample) {
        importance.ex.total <- numeric(forest@nexamples)
        importance.ex <- matrix(0, forest@nexamples, nvar)
    } else {
        importance.ex <- NULL
    }
    
    
    
    # Calculate variable importance --------------------------------------------
    
    for (tree in forest) {
        # Count correct oob classifications.
        w <- which(tree@insamp == 0L & tree@trainpredclass == y)
        correct <- sum(tree@nodewt[tree@trainprednode[w]])
        nout <- sum(tree@insamp == 0L)
        whichout <- which(tree@insamp == 0L)
        rm(w)
        
        # Variable importance for each example.
        if (impbyexample) {
            w <- which(tree@trainpredclass[whichout] == y[whichout])
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
            
            w <- which(oobpredict.result$oobpredclass == y[whichout])
            impcorrect <- sum(tree@nodewt[oobpredict.result$oobprednode[w]])
            rm(w)
            
            if (impbyexample) {
                w <- which(oobpredict.result$oobpredclass == y[whichout])
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
            w <- which(tree@trainpredclass[whichout] == y[whichout])
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
    significance <- 2 * pnorm(zscore * sqrt(2), lower = FALSE)
    
    if (impbyexample) {
        for (var in seq_len(nvar)) {
            importance.ex[, var] <- 100 *
                (importance.ex.total - importance.ex[, var]) / forest@ntrees
        }
    }
    
    return(list(importance=importance, importance.ex=importance.ex,
                zscore=zscore, significance=significance))
}

setMethod("varimp", signature(forest="bigcforest"),  varimp.bigcforest)
