\name{interactions-methods}
\docType{methods}
\alias{interactions}
\alias{interactions-methods}
\alias{interactions,bigcforest-method}
\title{Compute Variable Interactions}
\description{
  Compute the interactions between variables, using an experimental method described by Breiman and Cutler. Variables m and n interact if a split on variable m in a tree makes a split on n either systematically less possible or more possible.
}
\usage{
\S4method{interactions}{bigcforest}(forest)
}
\arguments{
  \item{forest}{A random forest of class \code{"\linkS4class{bigcforest}"}.}
}
\value{
  A symmetrical matrix of interactions between variables. A large positive number indicates that a split on one variable inhibits a split on the other variable, and conversely. This could indicate that the two variables are strongly correlated.
}
\section{Methods}{
\describe{
  \item{\code{signature(forest = "bigcforest")}}{Compute variable interactions for a classification random forest.}
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

# Calculate variable interactions.
inter <- interactions(forest)
}
\keyword{methods}
