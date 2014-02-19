\name{varimp-methods}
\docType{methods}
\alias{varimp}
\alias{varimp-methods}
\alias{varimp,bigcforest-method}
\title{Compute Variable Importance}
\description{
  Compute variable importance based on out-of-bag estimates. For each tree in the forest, the predictions of the out-of-bag examples are recorded. Then, a variable v is randomly permuted in the out-of-bag examples, and the tree is used to classify the out-of-bag examples again. The difference in votes for the correct class in the original data and the permuted data is used to calculate the variable importance for variable v. This process is then repeated for all variables.
}
\usage{
\S4method{varimp}{bigcforest}(forest, x=NULL, impbyexample=FALSE,
    reuse.cache=FALSE, trace=0L)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{x}{A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} of predictor variables. The data must not have changed, otherwise unexpected modelling results may occur. If a \code{matrix} or \code{data.frame} is specified, it will be converted into a \code{big.matrix} for computation. Optional if \code{reuse.cache} is \code{TRUE}.}
  \item{impbyexample}{A logical indicating whether to compute the variable importance for each out-of-bag example.}
  \item{reuse.cache}{\code{TRUE} to reuse disk caches of the \code{big.matrix} \code{x} from the initial building of the random forest, which may significantly reduce initialization time for large data sets. If \code{TRUE}, the user must ensure that the files \file{x} and \file{x.desc} in \code{forest@cachepath} have not been modified or deleted.}
  \item{trace}{\code{0} for no verbose output. \code{1} to print verbose output. Default: \code{0}.}
}
\value{
  A list with the following components:
  \describe{
    \item{\code{importance}:}{Importance of each variable, which is the number of votes for the correct class in the out-of-bag examples with variable v permuted subtracted from the number of votes for the correct class in the original out-of-bag examples, averaged over all trees in the forest.}
    \item{\code{importance.ex}:}{Importance of each variable for each out-of-bag example.}
    \item{\code{zscore}:}{Z-score of each variable, computed by dividing the raw variable importance score by the standard error.}
    \item{\code{significance}:}{Significance level of each variable importance, computed by applying the complementary error function on the z-score.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest")}}{Compute variable importance for a classification random forest.}
  }
}
\references{
  Breiman, L. (2001). Random forests. \emph{Machine learning}, \bold{45(1)}, 5-32.

  Breiman, L. & Cutler, A. (n.d.). Random Forests. Retrieved from \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm}.
}
\seealso{
  \code{\link[bigrf]{fastimp}}
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
importance <- varimp(forest, x, impbyexample=TRUE)
}
\keyword{methods}
