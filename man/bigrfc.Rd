\name{bigrfc}
\alias{bigrfc}
\title{
  Random Forest Classification Model
}
\description{
  Train a random forest classification model using Leo Breiman and Adele Cutler's algorithm, with enhancements for large data sets. This implementation uses the \pkg{\link[=bigmemory-package]{bigmemory}} package for disk-based caching during tree building, and the \pkg{\link[=foreach-package]{foreach}} package to parallelize the tree building process.
}
\usage{
  bigrfc(x, y = NULL, ntrees = 500L, supervised = !is.null(y),
    varselect = NULL, varnlevels = NULL,
    nsplitvar = round(sqrt(ifelse(is.null(varselect), ncol(x),
    length(varselect)))), maxeslevels = 11L, nrandsplit = 1023L,
    maxndsize = 1L, yclasswts = NULL, printerrfreq = 10L,
    printclserr = TRUE, cachepath = tempdir(), trace = 0L)
}
\arguments{
  \item{x}{
    A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} of predictor variables. If a \code{matrix} or \code{data.frame} is specified, it will be converted into a \code{big.matrix} for computation.
  }
  \item{y}{
    An integer or factor vector of response variables. Must be specified for supervised learning.
  }
  \item{ntrees}{
    The number of trees to be grown in the forest.
  }
  \item{supervised}{
    \code{TRUE} for supervised learning (training set has class labels), \code{FALSE} for unsupervised learning (training set has no class labels). Default: \code{TRUE} if \code{y} is specified, FALSE otherwise.
  }
  \item{varselect}{
    An integer vector specifying which columns in \code{x} to use. If not specified, all variables will be used.
  }
  \item{varnlevels}{
    An integer vector with elements specifying the number of levels in the corresponding variables in use, or 0 for numeric variables. Used only when \code{x} does not contain levels information (i.e. \code{x} is a \code{matrix} or \code{big.matrix}). If \code{x} is a \code{data.frame}, \code{varnlevels} will be inferred from \code{x}. If \code{x} is not a \code{data.frame} and \code{varnlevels} is \code{NULL}, all variables will be treated as numeric. If all columns of \code{x} are used, \code{varnlevels} should have as many elements as there are columns of \code{x}. But if varselect is specified, then \code{varnlevels} and \code{varselect} should be of the same length.
  }
  \item{nsplitvar}{
    The number of variables to split on at each node. Default: If \code{varselect} is specified, the square root of the number of variables specified; otherwise, the square root of the number of columns of \code{x}.
  }
  \item{maxeslevels}{
    Maximum number of levels for categorical variables for which exhaustive search of possible splits will be performed. Default: 11. This will amount to searching (2 ^ (11 - 1)) - 1 = 1,023 splits.
  }
  \item{nrandsplit}{
    Number of random splits to examine for categorical variables with more than maxeslevels levels. Default: 1,023.
  }
  \item{maxndsize}{
    Maximum number of examples in each node when growing the trees. Nodes will be split if they have more than this number of examples. Default: 1.
  }
  \item{yclasswts}{
    For supervised learning, a numeric vector of class weights, or \code{NULL} if all classes should be weighted equally. Not used for unsupervised learning.
  }
  \item{printerrfreq}{
    An integer, specifying how often error estimates should be printed to the screen. Default: error estimates will be printed every 10 trees.
  }
  \item{printclserr}{
    \code{TRUE} for error estimates for individual classes to be printed, in addition to the overall error estimates. Default: \code{TRUE}.
  }
  \item{cachepath}{
    Path to folder where data caches used in building the forest can be stored. If \code{NULL}, then the \code{big.matrix}'s will be created in memory, with no disk caching; this will be suitable for small data sets. Default: \code{tempdir()}.
  }
  \item{trace}{
    \code{0} for no verbose output. \code{1} to print verbose output on tree building. \code{2} to print more verbose output on processing of individual nodes. Default: \code{0}. Due to the way \code{\link{\%dopar\%}} handles the output of the tree-building iterations, you may not see the verbose output in some GUIs like RStudio. For best results, run R from the command line in order to see all the verbose output.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class \code{"\linkS4class{bigcforest}"} containing the specified number of trees, which are objects of class \code{"\linkS4class{bigctree}"}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[randomForest]{randomForest}}
\code{\link[party]{cforest}}
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
}
\keyword{models}
\keyword{tree}
