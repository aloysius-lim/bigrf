\name{outliers-methods}
\docType{methods}
\alias{outliers}
\alias{outliers-methods}
\alias{outliers,bigcforest-method}
\title{Compute Outlier Scores}
\description{
  Compute outlier scores for each class of examples used to train a random forest. Outliers are defined as examples whose proximities to other examples in the same class are small.
}
\usage{
\S4method{outliers}{bigcforest}(forest, trace=0L)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{trace}{\code{0} for no verbose output. \code{1} to print verbose output. Default: \code{0}.}
}
\value{
  A numeric vector containing the outlier scores for each training example. Higher scores indicate greater dissimilarity from other training examples in the same class.
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest")}}{Compute outlier scores for a classification random forest.}
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

# Calculate proximity matrix and scaling co-ordinates, and plot
# them.
prox <- proximities(forest, cachepath=NULL)
scale <- scaling(prox)
plot(scale, col=as.integer(y) + 2, pch=as.integer(y) + 2)

# Calculate outlier scores, and circle the top 20\% percent of
# them in red.
outscores <- outliers(forest)
points(scale[outscores > quantile(outscores, probs=0.8), ],
    col=2, pch=1, cex=1.5)
}
\keyword{methods}
