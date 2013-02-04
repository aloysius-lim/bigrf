bigrfc <- function(x,
                   y=NULL,
                   ntrees=500L,
                   supervised=!is.null(y),
                   varselect=NULL,
                   varnlevels=NULL,
                   nsplitvar=round(sqrt(ifelse(is.null(varselect), ncol(x),
                                               length(varselect)))),
                   maxeslevels=11L,
                   nrandsplit=1023L,
                   maxndsize=1L,
                   yclasswts=NULL,
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
    if (trace < 0L || trace > 2L) {
        stop("Argument trace must be 0, 1 or 2.")
    }
    
    if (trace >= 1L) message("Checking arguments in bigrfc.")
    
    # Check supervised.
    if (!is.logical(supervised)) {
        stop("Argument supervised must be a logical.")
    }
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix, matrix or data.frame.")
    }
    
    # Check y, and set ynclass, the number of classes in the response variable.
    # Also, if y is a factor vector, set ylevels, the original factor labels.
    if (supervised) {
        if (is.null(y)) {
            stop("Argument y must be specified for supervised learning.")
        }
        if (is.factor(y)) {
            ylevels <- levels(y)
            ynclass <- length(ylevels)
            y <- as.integer(y)
        } else if (is.integer(y)) {
            ylevels <- as.character(unique(y))
            ynclass <- length(unique(y))
        } else {
            stop("Argument y must be a factor or integer vector of class ",
                 "labels.")
        }
        if (length(y) != nrow(x)) {
            stop("Argument y must have as many elements as there are rows in ",
                 "x.")
        }
    } else {
        ylevels = c("original", "synthesized")
        ynclass = 2L;
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

    # Check varnlevels, and set factorvars, a logical vector indicating which
    # variables of x are factors.
    if (class(x) == "data.frame") {
        if (!is.null(varnlevels)) {
            warning("Argument varnlevels argument is ignored for data.frames. ",
                    "varnlevels will be inferred from columns of x instead.")
        }
        factorvars <- sapply(x[varselect], is.factor)
        varnlevels <- integer(length(varselect))
        varnlevels[factorvars] <- as.integer(sapply(x[varselect][factorvars],
                                                 function(n) length(levels(n))))
    } else if (is.null(varnlevels)) {
        warning("x is not a data.frame and argument varnlevels is not ",
                "specified. Number of levels cannot be inferred so all ",
                "variables will be treated as numeric.")
        factorvars <- logical(length(varselect))
        varnlevels <- numeric(length(varselect))
    } else {
        if (!is.integer(varnlevels)) {
            stop("Argument varnlevels must be an integer vector.")
        }
        if (length(varnlevels) != length(varselect)) {
            stop("Argument varnlevels must have the same number of elements as ",
                 "the number of variables being used.")
        }
        if (any(varnlevels == 1L)) {
            stop("Variable(s) ", paste(which(varnlevels == 1L), collapse=", "),
                 " have only 1 level each, and cannot be used to split nodes.")
        }
        factorvars <- varnlevels > 0L
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
    
    # Check yclasswts.
    if (!is.null(y)) {
        if (is.null(yclasswts)) {
            yclasswts <- rep.int(1L, ynclass)
        } else {
            if (!is.numeric(yclasswts)) {
                stop("Argument yclasswts must be a numeric vector.")
            }
            if (length(yclasswts) != ynclass) {
                stop("Argument yclasswts must have the same number of ",
                     "elements as there are levels in y.")
            }
        }
    }
    yclasswts <- scale(yclasswts, center=FALSE)
    
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
    
    
    
    # Initialize parameters ----------------------------------------------------
    
    if (trace >= 1L) message("Initializing parameters in bigrfc.")
    
    # Convert x to big.matrix, as C functions only support this at the moment.
    if (class(x) != "big.matrix") {
        if (is.null(cachepath)) {
            x <- as.big.matrix(x)
        } else {
            x <- as.big.matrix(x, backingfile="x", descriptorfile="x.desc",
                               backingpath=cachepath)
        }
    }
    
    forest <- new("bigcforest",
                  supervised=supervised,
                  varselect=varselect,
                  factorvars=factorvars,
                  varnlevels=varnlevels,
                  ylevels=ylevels,
                  ytable=table(0),
                  ynclass=ifelse(supervised, ynclass, 2L),
                  yclasswts=yclasswts,
                  ntrees=0L,
                  nsplitvar=nsplitvar,
                  maxndsize=maxndsize,
                  maxeslevels=maxeslevels,
                  nrandsplit=nrandsplit,
                  trainconfusion=table(0),
                  cachepath=cachepath)
    rm(supervised, factorvars, varnlevels, varselect, ynclass, yclasswts,
       nsplitvar, maxndsize, maxeslevels, nrandsplit, cachepath)
    
    # Number of examples to be used to build model.
    forest@nexamples <- ifelse(forest@supervised, as.integer(nrow(x)),
                             2L * as.integer(nrow(x)))
    
    # Maximum number of nodes in any given tree.
    # Theoretically, 2*nexamples - 1 should be enough
    forest@maxnodes <- 2L * forest@nexamples + 1L
    
    # Sequence number of categorical variables. Used to index columns of a and
    # a.out later.
    forest@contvarseq <- integer(length(forest@factorvars))
    forest@contvarseq[!forest@factorvars] <- seq_len(sum(!forest@factorvars))
    
    # if (supervised) {
    #   wtx <- yclasswts[y]
    # } else {
    #   wtx <- rep.int(1, nexamples)
    # }
    
    # Number of trees for which the each example has been out-of-bag.
    forest@oobtimes <- integer(forest@nexamples)
    forest@oobvotes <- matrix(0, forest@nexamples, forest@ynclass)
    forest@oobpred <- integer(forest@nexamples)
    forest@trainclserr <- numeric(forest@ynclass)
    forest@varginidec <- numeric(length(forest@varselect))
    
    
    
    # Grow forest --------------------------------------------------------------
    
    forest <- grow.bigcforest(forest, x, y, ntrees, printerrfreq=printerrfreq,
                              printclserr=printclserr, trace=trace)
    
    return(forest)
}
