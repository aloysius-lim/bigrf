\name{fastimp-methods}
\docType{methods}
\alias{fastimp}
\alias{fastimp-methods}
\alias{fastimp,bigcforest-method}
\title{Compute Fast (Gini) Variable Importance}
\description{
  Calculates variable importance using a fast method, by adding up the decreases in Gini impurity for each variable over all trees. The results are often consistent with the full variable importance calculated with \code{\link[bigrf]{varimp}}.
}
\usage{
\S4method{fastimp}{bigcforest}(forest)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
}
\value{
A numeric vector containing the fast (Gini) variable importance measures for each variable.
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest")}}{Compute the fast (Gini) variable importance for a classification random forest.}
  }
}
\references{
  Breiman, L. (2001). Random forests. \emph{Machine learning}, \bold{45(1)}, 5-32.

  Breiman, L. & Cutler, A. (n.d.). Random Forests. Retrieved from \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm}.
}
\seealso{
  \code{\link[bigrf]{varimp}}
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

# Calculate variable importance, including those for each out-of-bag example.
fastimportance <- fastimp(forest)
}
\keyword{methods}
