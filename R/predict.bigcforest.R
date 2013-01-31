predict.bigcforest <- function(object, x, y=NULL,
                               printerrfreq=10L,
                               printclserr=TRUE,
                               cachepath=tempdir(),
                               trace=FALSE) {
    # Check arguments ----------------------------------------------------------
    
    # Check trace.
    if (!is.logical(trace)) {
        stop ("Argument trace must be a logical.")
    }
    
    if (trace) message("Checking arguments.")
    
    # Check forest.
    forest <- object
    if (!class(forest) == "bigcforest") {
        stop("Argument forest must be a bigcforest created with bigrfc.")
    }
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix, matrix or data.frame.")
    }
    
    # Check y.
    if (!is.null(y)) {
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
    
    # Convert x to big.matrix, as C functions only support this at the moment.
    if (class(x) != "big.matrix") {
        if (is.null(forest@cachepath)) {
            x <- as.big.matrix(x)
        } else {
            x <- as.big.matrix(x, backingfile="xtest",
                               descriptorfile="xtest.desc",
                               backingpath=cachepath)
        }
    }
    
    ntest <- as.integer(nrow(x));
    xtype <- as.integer(.Call("CGetType", x@address, PACKAGE="bigmemory"))
    if (!is.null(y)) {
        class(y) <- "factor"
        levels(y) <- forest@ylevels
        ytable <- table(y, deparse.level=0)
        y <- as.integer(y)
    } else {
        ytable <- table(0);
    }
    
    # fast fix on the test data
    # if(missfill.eq.1) then
    # 	read(1,*) (fill(m),m=1,mdim)
    # 	call xfill(xts,ntest,mdim,fill,code)
    # endif

    prediction <- new("bigcprediction",
                      ntest=ntest,
                      testlabelled=!is.null(y),
                      nclass=forest@nclass,
                      ntrees=forest@ntrees,
                      testytable=ytable,
                      testvotes=matrix(0, ntest, forest@nclass),
                      testclserr=numeric(forest@nclass),
                      testerr=0,
                      testconfusion=table(0),
                      printerrfreq=printerrfreq,
                      printclserr=printclserr,
                      cachepath=cachepath
    )
    rm(ntest, ytable, printerrfreq, printclserr, cachepath)
    
    
    
    # Compute test results -----------------------------------------------------
    
    # Loop through all trees.
    for (t in seq_len(forest@ntrees)) {
        if (trace) message("Running tree ", t, " on test cases.")
        tree <- forest[[t]]

        treepredict.result <- .Call("treepredictC", x@address, xtype,
                                    prediction@ntest, forest, tree);
        
        # Compute votes.
        for (c in seq_len(prediction@nclass)) {
            w <- which(treepredict.result$testpredclass == c)
            prediction@testvotes[w, c] <- prediction@testvotes[w, c] +
                tree@nodewt[treepredict.result$testprednode[w]]
        }
        rm(c, w)
        prediction[seq_len(prediction@ntest)] <- max.col(prediction@testvotes)
        
        # If test set labels were given, compute test error.
        if (!is.null(y)) {
            prediction@testclserr <- integer(prediction@nclass)
            for (c in seq_len(prediction@nclass)) {
                prediction@testclserr[c] <-
                    sum(y == c & prediction[] != c)
            }
            prediction@testerr <- sum(prediction@testclserr) / prediction@ntest
            prediction@testclserr <-
                prediction@testclserr / as.numeric(prediction@testytable)
        }
        
        # Give running output --------------------------------------------------
        if (t == 1L) {
            cat("Test errors:\n")
            if (prediction@printclserr && !is.null(y)) {
                cat(" Tree  Overall error  Error by class\n")
                cat("                      ")
                cat(format(names(prediction@testytable), justify="right",
                           width=5),
                    sep="  ")
                cat("\n")
            } else {
                cat(" Tree  Overall error\n")
            }
        }
        
        if (t %% prediction@printerrfreq == 0L || t == forest@ntrees) {
            cat(format(t, justify="right", width=5),
                format(100 * prediction@testerr, justify="right", width=13,
                       digits=3, nsmall=2), sep="  ")
            if (!is.null(y) && prediction@printclserr) {
                cat("",
                    format(100 * prediction@testclserr, justify="right",
                           width=max(nchar(names(prediction@testytable)), 5),
                           digits=3,
                           nsmall=2),
                    sep="  ")
            }
            cat("\n")
        }
    }
    cat("\n")
    
    # Calculate confusion matrix -----------------------------------------------
    
    if (!is.null(y)) {
        pred <- prediction[]
        pred[pred == 0L] <- prediction@nclass + 1L
        if (length(forest@ylevels)) {
            class(pred) <- "factor"
            levels(pred) <- c(forest@ylevels, "Never out-of-bag")
            class(y) <- "factor"
            levels(y) <- forest@ylevels
        }
        prediction@testconfusion <- table(y, pred, dnn=c("Actual", "Predicted"))
        rm(pred)
    }
    
    summary(prediction)
    
    return(prediction)
}

setMethod("predict", signature(object="bigcforest"),  predict.bigcforest)
