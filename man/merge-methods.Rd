\name{merge-methods}
\docType{methods}
\alias{merge}
\alias{merge-methods}
\alias{merge,bigcforest,bigcforest-method}
\title{Merge Two Random Forests}
\description{
  Merge two random forests into one. This is useful, for example, for building random forests in parallel on multiple machines, then combining them into a single forest.
}
\usage{
\S4method{merge}{bigcforest,bigcforest}(x, y)
}
\arguments{
  \item{x}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{y}{A random forest of class \code{"\linkS4class{bigcforest}"}, built using the same data sources and forest-building parameters as \code{x}.}
}
\details{
  These methods copy all the trees from \code{y} into \code{x}, and calculates the error statistics and confusion matrices of the merged forest.
}
\value{
  The merged forest of the same class as the input arguments.
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "bigcforest", y = "bigcforest")}}{Merges the classification random forests x and y.}
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

# Run model, grow 15 trees.
forest <- bigrfc(x, y, ntree=15L, varselect=vars, cachepath=NULL)

# Build a second forest.
forest2 <- bigrfc(x, y, ntree=10L, varselect=vars, cachepath=NULL)

# Merge the two forests.
big.forest <- merge(forest, forest2)
}
\keyword{methods}
