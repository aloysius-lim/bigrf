\name{generateSyntheticClass}
\alias{generateSyntheticClass}
\title{Generate Synthetic Second Class for Unsupervised Learning}
\description{
  To use Random Forests for unsupervised learning, the training set \code{x} is treated as a single class. This function creates a synthetic second class for classification by sampling at random from the univariate distributions of the original data. This is useful, for example, for clustering.
}
\usage{
generateSyntheticClass(x, ...)
}
\arguments{
  \item{x}{A \code{\link{big.matrix}}, \code{matrix} or \code{data.frame} containing the predictor variables of the original training set.}
  \item{\dots}{If \code{x} is a \code{big.matrix}, these arguments will be passed on to \code{\link{big.matrix}} to control how the \code{big.matrix} for the two-class training set is created.}
}
\value{
  A list containing the following components:
  \item{x}{The two-class training set, comprising the original training set and the synthesized second class. It will be an object of the same type as the argument \code{x}.}
  \item{y}{A factor vector that labels the two classes in \code{x}.}
}
\references{
  Breiman, L. (2001). Random forests. \emph{Machine learning}, \bold{45(1)}, 5-32.

  Breiman, L. & Cutler, A. (n.d.). Random Forests. Retrieved from \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm}.
}
\examples{
# Perform unsupervised learning on the Cars93 data set.

# Load data.
data(Cars93, package="MASS")

# Create second synthetic class for unsupervised learning.
newdata <- generateSyntheticClass(Cars93)

# Select variables with which to train model.
vars <- c(4:22)

# Run model, grow 30 trees.
forest <- bigrfc(newdata$x, newdata$y, ntree=30L, varselect=vars,
                 cachepath=NULL)
}