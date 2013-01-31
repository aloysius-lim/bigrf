bigrfc <- function(x,
                   y=NULL,
                   ntrees=500L,
                   supervised=!is.null(y),
                   nlevels=NULL,
                   varselect=NULL,
                   nsplitvar=round(sqrt(ifelse(is.null(varselect), ncol(x),
                                               length(varselect)))),
                   maxeslevels=11L,
                   nrandsplit=1023L,
                   maxndsize=1L,
                   classweights=NULL,
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
    
    # Check supervised.
    if (!is.logical(supervised)) {
        stop("Argument supervised must be a logical.")
    }
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix, matrix or data.frame.")
    }
    
    # Check y, and set nclass, the number of classes in the response variable.
    # Also, if y is a factor vector, set ylevels, the original factor labels.
    if (supervised) {
        if (is.null(y)) {
            stop("Argument y must be specified for supervised learning.")
        }
        if (is.factor(y)) {
            ylevels <- levels(y)
            nclass <- length(ylevels)
            y <- as.integer(y)
        } else if (is.integer(y)) {
            ylevels <- character()
            nclass <- length(unique(y))
        } else {
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
    
    # Check varselect.
    if (is.null(varselect)) {
        varselect <- seq_len(ncol(x))
    } else {
        if (!is.integer(varselect)) {
            stop("Argument varselect must be an integer vector.")
        }
        if (any(duplicated(varselect))) {
            stop("Argument varselect cannot contain duplicate entries.")
        }
        if (max(varselect) > ncol(x)) {
            stop("Argument varselect cannot contain values greater than the ",
                 "number of columns in x.")
        }
        if (min(varselect) < 1L) {
            stop("Argument varselect cannot contain values less than 1.")
        }
    }

    # Check nlevels, and set factors, a logical vector indicating which
    # variables of x are factors.
    if (class(x) == "data.frame") {
        if (!is.null(nlevels)) {
            warning("Argument nlevels argument is ignored for data.frames. ",
                    "nlevels will be inferred from columns of x instead.")
        }
        factors <- sapply(x[varselect], is.factor)
        nlevels <- integer(length(varselect))
        nlevels[factors] <- as.integer(sapply(x[varselect][factors],
                                              function(n) length(levels(n))))
    } else if (is.null(nlevels)) {
        warning("x is not a data.frame and argument nlevels is not specified. ",
                "Number of levels cannot be inferred so all variables will be ",
                "treated as numeric.")
        factors <- logical(length(varselect))
        nlevels <- numeric(length(varselect))
    } else {
        if (!is.integer(nlevels)) {
            stop("Argument nlevels must be an integer vector.")
        }
        if (length(nlevels) != length(varselect)) {
            stop("Argument nlevels must have the same number of elements as ",
                 "the number of variables being used.")
        }
        if (any(nlevels == 1L)) {
            stop("Variable(s) ", paste(which(nlevels == 1L), collapse=", "),
                 " have only 1 level each, and cannot be used to split nodes.")
        }
        factors <- nlevels > 0L
    }
    
    # Check nsplitvar.
    if (!is.numeric(nsplitvar) ||
            abs(nsplitvar - round(nsplitvar)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument nsplitvar must be an integer.")
    }
    nsplitvar <- as.integer(round(nsplitvar))
    if (nsplitvar > length(varselect)) {
        stop("Argument nsplitvar cannot be greater than the number of ",
             "variables used for modelling.")
    }
    
    # Check maxeslevels.
    if (!is.numeric(maxeslevels) ||
            abs(maxeslevels - round(maxeslevels)) >=
            .Machine$double.eps ^ 0.5) {
        stop ("Argument maxeslevels must be an integer.")
    }
    maxeslevels <- as.integer(round(maxeslevels))
    if (maxeslevels < 2L) {
        stop("Argument maxeslevels must be at least 2.")
    }
    
    # Check nrandsplit.
    if (!is.numeric(nrandsplit) ||
            abs(nrandsplit - round(nrandsplit)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument nrandsplit must be an integer.")
    }
    nrandsplit <- as.integer(round(nrandsplit))
    if (nrandsplit < 1) {
        stop("Argument nrandsplit must be at least 1.")
    }
    
    # Check maxndsize.
    if (!is.numeric(maxndsize) ||
            abs(maxndsize - round(maxndsize)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument maxndsize must be an integer.")
    }
    maxndsize <- as.integer(round(maxndsize))
    if (maxndsize > nrow(x)) {
        stop("Argument maxndsize cannot be greater than the number of rows in ",
             "x.")
    }
    if (maxndsize < 1L) {
        stop("Argument maxndsize must be at least 1.")
    }
    
    # Check classweights.
    if (!is.null(y)) {
        if (is.null(classweights)) {
            classweights <- rep.int(1L, nclass)
        } else {
            if (!is.numeric(classweights)) {
                stop("Argument classweights must be a numeric vector.")
            }
            if (length(classweights) != nclass) {
                stop("Argument classweights must have the same number of ",
                     "elements as there are levels in y.")
            }
        }
    }
    classweights <- scale(classweights, center=FALSE)
    
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
    if (!is.character(cachepath)) {
        stop("Argument cachepath must be a character string.")
    }
    if (!file.exists(cachepath)) {
        if (!dir.create(cachepath)) {
            stop("Cannot create directory ", cachepath, ".")
        }
    }
    
    
    
    # Initialize parameters ----------------------------------------------------
    
    if (trace) message("Initializing parameters.")
    
    # Convert x to big.matrix, as C functions only support this at the moment.
    if (class(x) != "big.matrix") {
        x <- as.big.matrix(x)
    }
    
    forest <- new("bigcforest",
                  supervised=supervised,
                  factors=factors,
                  ylevels=ylevels,
                  nlevels=nlevels,
                  varselect=varselect,
                  nclass=ifelse(supervised, nclass, 2L),
                  classweights=classweights,
                  ntrees=0L,
                  nsplitvar=nsplitvar,
                  maxndsize=maxndsize,
                  maxeslevels=maxeslevels,
                  nrandsplit=nrandsplit,
                  trainconfusion=table(0),
                  printerrfreq=printerrfreq,
                  printclserr=printclserr,
                  cachepath=cachepath)
    rm(supervised, factors, nlevels, varselect, nclass, classweights, nsplitvar,
       maxndsize, maxeslevels, nrandsplit, printerrfreq, printclserr, cachepath)
    
    # The number of cases in the test set, or 0 if there is no test set.
    # ntest <- 0
    
    # 0 if the test set has no class labels, 1 if the test set has class labels. 
    # labelts <- 0
    
    # Number of samples to be used to build model.
    forest@nsample <- ifelse(forest@supervised, as.integer(nrow(x)),
                             2L * as.integer(nrow(x)))
    
    # Maximum number of nodes in any given tree.
    # Theoretically, 2*nsample - 1 should be enough
    forest@maxnodes <- 2L * forest@nsample + 1L
    
    # Sequence number of categorical variables. Used to index columns of a and
    # a.out later.
    forest@contvarseq <- integer(length(forest@factors))
    forest@contvarseq[!forest@factors] <- seq_len(sum(!forest@factors))
    
    # iftest <- ntest / (ntest - .1)
    
    # ntest0 <- (1L - iftest) + ntest
    
    # if (supervised) {
    #   wtx <- classweights[y]
    # } else {
    #   wtx <- rep.int(1, nsample)
    # }
    
    # Number of trees for which the each sample has been out-of-bag.
    forest@oobtimes <- integer(forest@nsample)
    forest@oobvotes <- matrix(0, forest@nsample, forest@nclass)
    forest@oobpred <- integer(forest@nsample)
    forest@trainclserr <- numeric(forest@nclass)
    # if (ntest > 0) {
    #    qts <- matrix(0, nclass, ntest)
    # }
    forest@avgini <- numeric(length(forest@varselect))
    
    
    
    # Grow forest --------------------------------------------------------------
    
    forest <- grow.bigcforest(forest, x, y, ntrees, trace=trace)
    
    return(forest)
}
