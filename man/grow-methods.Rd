\name{grow-methods}
\docType{methods}
\alias{grow}
\alias{grow-methods}
\alias{grow,bigcforest-method}
\title{Grow More Trees in a Random Forest}
\description{
  Grow more trees in a random forest, using the same parameters.
}
\usage{
\S4method{grow}{bigcforest}(forest, x=NULL, ntrees=50L, printerrfreq=10L,
    printclserr=TRUE, reuse.cache=FALSE, trace=0L)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{x}{A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} of predictor variables. The data must not have changed, otherwise unexpected modelling results may occur. If a \code{matrix} or \code{data.frame} is specified, it will be converted into a \code{big.matrix} for computation. Optional if \code{reuse.cache} is \code{TRUE}.}
  \item{ntrees}{The number of additional trees to grow in the forest. Default: 50.}
  \item{printerrfreq}{An integer, specifying how often error estimates should be printed to the screen. Default: error estimates will be printed every 10 trees.}
  \item{printclserr}{\code{TRUE} for error estimates for individual classes to be printed, in addition to the overall error estimates. Default: \code{TRUE}.}
  \item{reuse.cache}{\code{TRUE} to reuse disk caches of the \code{big.matrix}'s \code{x} and \code{asave} from the initial building of the random forest, which may significantly reduce initialization time for large data sets. If \code{TRUE}, the user must ensure that the files \file{x} and \file{x.desc} in \code{forest@cachepath} have not been modified or deleted. Files \file{asave} and \file{asave.desc} must also be present if there are any numeric variables.}
  \item{trace}{\code{0} for no verbose output. \code{1} to print verbose output on growing of trees, and a summary of the grown forest. \code{2} to print more verbose output on processing of individual nodes. Default: \code{0}. Due to the way \code{\link{\%dopar\%}} handles the output of the tree-growing iterations, you may not see the verbose output in some GUIs like RStudio. For best results, run R from the command line in order to see all the verbose output.}
}
\value{
  The object supplied as the \code{forest} argument, with additional trees grown.
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest")}}{Grow more trees in a classification random forest.}
  }
}
\references{
  Breiman, L. (2001). Random forests. \emph{Machine learning}, \bold{45(1)}, 5-32.

  Breiman, L. & Cutler, A. (n.d.). Random Forests. Retrieved from \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm}.
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

# Run model, grow 30 trees.
forest <- bigrfc(x, y, ntree=30L, varselect=vars, cachepath=NULL)

# Grow 10 more trees.
forest <- grow(forest, x, ntree=10L)
}
\keyword{methods}
\keyword{models}
\keyword{tree}
