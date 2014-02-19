\name{proximities-methods}
\docType{methods}
\alias{proximities}
\alias{proximities-methods}
\alias{proximities,bigcforest-method}
\title{Compute Proximity Matrix}
\description{
  Compute the proximity matrix for a random forest, for the \code{nnearest} most proximate examples to each training example.
}
\usage{
\S4method{proximities}{bigcforest}(forest, nnearest=forest@nexamples,
    cachepath=tempdir(), trace=0L)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{nnearest}{The number of most proximate examples for which to compute proximity measures for each training example. Setting this to a smaller number will speed up computation of scaling co-ordinates. Default: \code{forest@nexamples}.}
  \item{cachepath}{Path to folder where the proximity matrix can be stored. If \code{NULL}, then the \code{big.matrix}'s will be created in memory with no disk caching, which would be suitable for small data sets. If the user wishes to reuse the cached files, it is suggested that a folder other than \code{tempdir()} is used, as the operating system may automatically delete any cache files in \code{tempdir()}. Default: \code{tempdir()}.}
  \item{trace}{\code{0} for no verbose output. \code{1} to print verbose output. \code{2} to print even more verbose output on processing of each tree and example. Default: \code{0}.}
}
\value{
  An object of class \code{"\linkS4class{bigrfprox}"} containing the proximity matrix.
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest")}}{Compute the proximity matrix for a classification random forest.}
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

# Calculate proximity matrix.
prox <- proximities(forest, cachepath=NULL)
}
\keyword{methods}
