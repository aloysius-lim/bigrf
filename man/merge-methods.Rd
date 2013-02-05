\name{merge-methods}
\docType{methods}
\alias{merge}
\alias{merge-methods}
\alias{merge,bigcforest,bigcforest-method}
\title{Methods for Function \code{merge} in Package \pkg{bigrf}}
\description{Merges two random forests into one.}
\usage{
  \S4method{merge}{bigcforest,bigcforest}(x, y, class.labels)
}
\arguments{
  \item{x}{
    A random forest of class \code{"\linkS4class{bigcforest}"}.
  }
  \item{y}{
    A random forest of class \code{"\linkS4class{bigcforest}"}, built using the same data sources and forest-building parameters as \code{x}.
  }
  \item{class.labels}{
    The \code{factor} or \code{integer} vector of response variables (class labels). This is called \code{y} in most other functions in this package, but had to be renamed here to avoid clashing with the method signature for \code{merge}.
  }
}
\details{
  These methods copy all the trees from \code{y} into \code{x}, and calculates the error statistics and confusion matrices of the merged forest.
}
\value{
  The merged forest of the same class as the input arguments.
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "bigcforest", y = "bigcforest")}}{
  Merges the classification random forests x and y of class \code{"\linkS4class{bigcforest}"}.
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

# Build a second forest.
forest2 <- bigrfc(x, y, ntree=20L, varselect=vars)

# Merge the two forests.
big.forest <- merge(forest, forest2, y)
}
\keyword{methods}
