\name{prototypes-methods}
\docType{methods}
\alias{prototypes}
\alias{prototypes-methods}
\alias{prototypes,bigcforest,bigrfprox-method}
\title{Compute Class Prototypes}
\description{
  Compute the prototypes for each class in the training set, which provide a picture of how each variable relates to the classification. They are useful representations of a "typical" example of each class.
}
\usage{
\S4method{prototypes}{bigcforest,bigrfprox}(forest, prox, nprot=1L, x=NULL,
    reuse.cache=FALSE, trace=0L)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{prox}{A proximity matrix of class \code{"\linkS4class{bigrfprox}"}.}
  \item{nprot}{The number of prototypes to compute for each class. Default: \code{1}.}
  \item{x}{A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} of predictor variables. The data must not have changed, otherwise unexpected modelling results may occur. If a \code{matrix} or \code{data.frame} is specified, it will be converted into a \code{big.matrix} for computation. Optional if \code{reuse.cache} is \code{TRUE}.}
  \item{reuse.cache}{\code{TRUE} to reuse disk caches of the \code{big.matrix} \code{x} from the initial building of the random forest, which may significantly reduce initialization time for large data sets. If \code{TRUE}, the user must ensure that the files \file{x} and \file{x.desc} in \code{forest@cachepath} have not been modified or deleted.}
  \item{trace}{\code{0} for no verbose output. \code{1} to print verbose output. Default: \code{0}.}
}
\details{
  Prototypes are computed using proximities, as follows. For the first prototype for class c, find the example i with the largest number of class c examples among its k nearest neighbors. Among these examples, find the 25th, 50th and 75th percentiles of the numeric variables, and most frequent level of the categorical variables. For the second prototype, the procedure is repeated, considering only examples that are not among the k examples used to compute the first prototype, and so on.
}
\value{
  A list with the following components:
  \describe{
    \item{\code{nprotfound}:}{Number of prototypes found for each class.}
    \item{\code{clustersize}:}{\code{forest@ynclass} by \code{nprot} matrix indicating the number of examples used to compute each prototype.}
    \item{\code{prot}:}{\code{forest@ynclass} by \code{nprot} by \code{length(forest@varselect)} by \code{3} array containing the raw prototype values. For numeric variables, the prototypes are represented by the medians, with the 25th and 75th percentiles given as estimates of the prototype stability. For categorical variables, the values are the most frequent level.}
    \item{\code{prot.std}:}{\code{forest@ynclass} by \code{nprot} by \code{length(forest@varselect)} by \code{3} array containing standardized prototype values. Prototype values for numeric variables are subtracted by the 5th percentile, then divided by the difference between the 95th and 5th percentile. Prototype values for categorical variables are divided by the number of levels in that variable.}
    \item{\code{levelsfreq}:}{List of length \code{length(forest@varselect)} containing, for each categorical variable v, an \code{forest@ynclass} by \code{nprot} by \code{forest@varnlevels[v]} array that indicate the frequency of levels used to compute the prototype level. These are useful for estimating prototype stability for categorical variables.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{signature(forest = "bigcforest", prox = "bigrfprox")}}{Compute prototypes for a classification random forest.}
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

# Compute prototypes.
prot <- prototypes(forest, prox, x=x)

# Plot first prototypes, using one colour for each class.
plot(seq_along(vars), prot$prot.std[1, 1, , 2], type="l", col=1,
     ylim=c(min(prot$prot.std[, 1, , 2]), max(prot$prot.std[, 1, , 2])))
for (i in 2:length(levels(y))) {
    lines(seq_along(vars), prot$prot.std[i, 1, , 2], type="l", col=i)
}

# Plot first prototype for class 1, including quartile values for numeric
# variables.
plot(seq_along(vars), prot$prot.std[1, 1, , 1], type="l", col=1,
     ylim=c(min(prot$prot.std[1, 1, , ]), max(prot$prot.std[1, 1, , ])))
for (i in 2:3) {
    lines(seq_along(vars), prot$prot.std[1, 1, , i], type="l", col=i)
}
}
\keyword{methods}
