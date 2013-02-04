predict.bigcforest <- function(object, x, y=NULL,
                               printerrfreq=10L,
                               printclserr=TRUE,
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
        if (!forest@supervised && any(y != 0L | y != 1L)) {
            stop("Forest was built with unsupervised learning. y must contain ",
                 "only 0s (for original data) or 1s (for synthesized data.")
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
        if (is.null(cachepath)) {
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
        ytable <- NULL
    }
    
    # fast fix on the test data
    # if(missfill.eq.1) then
    # 	read(1,*) (fill(m),m=1,mdim)
    # 	call xfill(xts,ntest,mdim,fill,code)
    # endif

    prediction <- new("bigcprediction",
                      ntest=ntest,
                      testlabelled=!is.null(y),
                      ynclass=forest@ynclass,
                      ntrees=forest@ntrees,
                      testytable=ytable,
                      testvotes=matrix(0, ntest, forest@ynclass),
                      testclserr=if(is.null(y)) NULL else
                          numeric(forest@ynclass),
                      testerr=if(is.null(y)) NULL else 0,
                      testconfusion=NULL
    )
    rm(ntest, ytable)
    
    
    
    # Compute test results -----------------------------------------------------
    
    # Loop through all trees.
    prediction <- foreach(t=seq_len(forest@ntrees),
                          .combine=combine.treepredictresults, .init=prediction,
                          .inorder=FALSE, .verbose=FALSE) %dopar% {
        if (trace >= 1L) message("Running tree ", t, " on test examples.")
        tree <- forest[[t]]

        treepredict.result <- .Call("treepredictC", x@address, xtype,
                                    prediction@ntest, forest, tree);
        treepredict.result$t <- t
        treepredict.result$y <- y
        treepredict.result$forest <- forest
        treepredict.result$tree <- tree
        treepredict.result$printerrfreq <- printerrfreq
        treepredict.result$printclserr <- printclserr
        treepredict.result
    }
    cat("\n")
    
    # Calculate confusion matrix -----------------------------------------------
    
    if (!is.null(y)) {
        pred <- prediction[]
        if (length(forest@ylevels)) {
            class(pred) <- "factor"
            levels(pred) <- forest@ylevels
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
