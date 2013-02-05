\name{grow-methods}
\docType{methods}
\alias{grow}
\alias{grow-methods}
\alias{grow,bigcforest-method}
\title{Methods for Function \code{grow} in Package \pkg{bigrf}}
\description{
 Grows more trees in a random forest.
}
\usage{
  \S4method{grow}{bigcforest}(forest, x=NULL, y, ntrees=500L, printerrfreq=10L,
    printclserr=TRUE, reuse.cache=FALSE, trace=0L)
}
\arguments{
  \item{forest}{
    A random forest of class \code{"\linkS4class{bigcforest}"}.
  }
  \item{x}{
    A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} of predictor variables. The data must not have changed, otherwise unexpected modelling results may occur. If a \code{matrix} or \code{data.frame} is specified, it will be converted into a \code{big.matrix} for computation.
  }
  \item{y}{
    An integer or factor vector of response variables that was used in growing the forest. The data must not have changed, otherwise unexpected modelling results may occur.
  }
  \item{ntrees}{
    The number of additional trees to grow in the forest.
  }
  \item{printerrfreq}{
    An integer, specifying how often error estimates should be printed to the screen. Default: error estimates will be printed every 10 trees.
  }
  \item{printclserr}{
    \code{TRUE} for error estimates for individual classes to be printed, in addition to the overall error estimates. Default: \code{TRUE}.
  }
  \item{reuse.cache}{
    \code{TRUE} to reuse disk caches from previous tree-growing runs, which may significantly reduce initialization time for large data sets. Currently, only the \code{big.matrix}'s \code{asave}, \code{b} and \code{x} are reused. If \code{TRUE}, user must ensure that the files \file{asave}, \file{asave.desc}, \file{b}, \file{b.desc}, \file{x} and \file{x.desc} in \code{forest@cachepath} have not been modified or deleted since the previous run.
  }
  \item{trace}{
    \code{0} for no verbose output. \code{1} to print verbose output on tree building. \code{2} to print more verbose output on processing of individual nodes. Default: \code{0}. Due to the way \code{\link{\%dopar\%}} handles the output of the tree-building iterations, you may not see the verbose output in some GUIs like RStudio. For best results, run R from the command line in order to see all the verbose output.
  }
}
\value{
  The object supplied as the \code{forest} argument, with additional trees grown.
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest")}}{
      Grow more trees in the given classification forest.
    }
  }
}
\examples{
# Classify cars in the Cars93 data set by type (Compact, Large,
# Midsize, Small, Sporty, or Van).

# Load data.
data(Cars93, package="MASS")
x <- Cars93
y <- Cars93$Type

# Select variables with which to train model.
vars <- c(4:22)

# Run model, build 50 trees.
forest <- bigrfc(x, y, ntree=50L, varselect=vars)

# Grow 30 more trees.
forest <- grow(forest, x, y, ntree=30L)
}
\keyword{methods}
\keyword{models}
\keyword{tree}
