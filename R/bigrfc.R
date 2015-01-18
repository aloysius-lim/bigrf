bigrfc <- function(x,
                   y,
                   ntrees=50L,
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
    
    # Check x.
    if (!(class(x) %in% c("big.matrix", "matrix", "data.frame"))) {
        stop("Argument x must be a big.matrix, matrix or data.frame.")
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
    ytable <- table(y, deparse.level=0)
    
    # Check ntrees.
    if (!is.numeric(ntrees) ||
            abs(ntrees - round(ntrees)) >= .Machine$double.eps ^ 0.5) {
        stop ("Argument ntrees must be an integer.")
    }
    ntrees <- as.integer(round(ntrees))
    if (ntrees < 0L) {
        stop("Argument ntrees must not be negative.")
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
    names(varselect) <- dimnames(x)[[2]][varselect]

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
        # Logicals should be considered factors with 2 levels.
        logicalvars <- sapply(x[varselect], is.logical)
        factorvars[logicalvars] <- TRUE
        varnlevels[logicalvars] <- 2L
        rm(logicalvars)
    } else if (is.null(varnlevels)) {
        warning("x is not a data.frame and argument varnlevels is not ",
                "specified. Number of levels cannot be inferred so all ",
                "variables will be treated as numeric.")
        factorvars <- logical(length(varselect))
        varnlevels <- integer(length(varselect))
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
    names(factorvars) <- names(varselect)
    names(varnlevels) <- names(varselect)
    
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
    if (is.null(yclasswts)) {
        yclasswts <- rep.int(1L, length(levels(y)))
    } else {
        if (!is.numeric(yclasswts)) {
            stop("Argument yclasswts must be a numeric vector.")
        }
        if (length(yclasswts) != length(levels(y))) {
            stop("Argument yclasswts must have the same number of ",
                 "elements as there are levels in y.")
        }
    }
    yclasswts <- scale(yclasswts, center=FALSE)
    dimnames(yclasswts) <- list(Class=levels(y), NULL)
    
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
    
    forest <- new("bigcforest",
                  nexamples=as.integer(nrow(x)),
                  varselect=varselect,
                  factorvars=factorvars,
                  varnlevels=varnlevels,
                  y=y,
                  ytable=ytable,
                  yclasswts=yclasswts,
                  ntrees=0L,
                  nsplitvar=nsplitvar,
                  maxndsize=maxndsize,
                  maxeslevels=maxeslevels,
                  nrandsplit=nrandsplit,
                  trainconfusion=table(0),
                  cachepath=cachepath)
    rm(factorvars, varnlevels, varselect, ytable, yclasswts, nsplitvar,
       maxndsize, maxeslevels, nrandsplit, cachepath)
    
    # Sequence number of continuous variables. Used to index columns of a and
    # a.out later.
    forest@contvarseq <- integer(length(forest@factorvars))
    forest@contvarseq[!forest@factorvars] <- seq_len(sum(!forest@factorvars))
    names(forest@contvarseq) <- names(forest@varselect)
    
    # Out-of-bag results and error estimates.
    forest@oobtimes <- integer(forest@nexamples)
    forest@oobvotes <- matrix(0, forest@nexamples, length(levels(y)),
                              dimnames=list(Example=NULL,
                              Class=levels(forest@y)))
    forest@oobpred <- integer(forest@nexamples)
    forest@trainerr <- numeric()
    forest@trainclserr <- matrix(0, 0, length(levels(y)),
                                 dimnames=list(NTrees=NULL,
                                               Class=levels(forest@y)))
    forest@varginidec <- numeric(length(forest@varselect))
    names(forest@varginidec) <- names(forest@varselect)
    
    
    
    # Build forest -------------------------------------------------------------
    
    if (ntrees > 0L) {
        forest <- grow(forest, x, ntrees, printerrfreq=printerrfreq,
                       printclserr=printclserr, trace=trace)
    }
    
    return(forest)
}
