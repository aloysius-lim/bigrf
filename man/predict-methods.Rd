\name{predict-methods}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,bigcforest-method}
\title{Predict Classes of Test Examples}
\description{
  Predict the classes of a set of test examples using a random forest.
}
\usage{
\S4method{predict}{bigcforest}(object, x, y=NULL, printerrfreq=10L,
    printclserr=TRUE, cachepath=tempdir(), trace=0L)
}
\arguments{
  \item{object}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
  \item{x}{A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} of predictor variables. If a \code{matrix} or \code{data.frame} is specified, it will be converted into a \code{big.matrix} for computation.}
  \item{y}{An integer or factor vector of response variables. Test errors will be calculated only if \code{y} is specified.}
  \item{printerrfreq}{An integer, specifying how often error estimates should be printed to the screen. Default: error estimates will be printed every 10 trees.}
  \item{printclserr}{\code{TRUE} for error estimates for individual classes to be printed, in addition to the overall error estimates. Default: \code{TRUE}.}
  \item{cachepath}{Path to folder where data caches used in building the forest can be stored. If \code{NULL}, then the \code{big.matrix}'s will be created in memory, with no disk caching; this will be suitable for small data sets. Default: \code{tempdir()}.}
  \item{trace}{\code{0} for no verbose output. \code{1} to print verbose output on prediction by trees, and a summary of the predictions. Default: \code{0}. Due to the way \code{\link{\%dopar\%}} handles output, you may not see the verbose output in some GUIs like RStudio. For best results, run R from the command line in order to see all the verbose output.}
}
\details{
  These methods copy all the trees from \code{y} into \code{x}, and calculates the error statistics and confusion matrices of the merged forest.
}
\value{
  An object of class \code{"\linkS4class{bigcprediction}"} containing the prediction results.
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "bigcforest")}}{Predict classes of a set of test examples using a classification random forest.}
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

# Run model, grow 30 trees on the first 60 examples.
forest <- bigrfc(x[1:60, ], y[1:60], ntree=30L, varselect=vars, cachepath=NULL)

# Get predictions for the remaining examples.
predictions <- predict(forest, x[-(1:60), ], y[-(1:60)], cachepath=NULL)
}
\keyword{methods}
