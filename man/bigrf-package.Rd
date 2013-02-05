\name{bigrf-package}
\alias{bigrf-package}
\alias{bigrf}
\docType{package}
\title{Big Random Forests--Classification and Regression Forests for Large Data Sets}
\description{
This is an implementation of Leo Breiman and Adele Cutler's Random Forest algorithms for classification and regression, optimized for data sets that are too large to fit or be processed in memory. It utilizes \pkg{\link[bigmemory:bigmemory-package]{bigmemory}} disk-based memory for caching of intermediate computations, and \pkg{\link[foreach:foreach-package]{foreach}} to parallelize the tree-buliding process.
}
\details{
\tabular{ll}{
Package: \tab bigrf\cr
Version: \tab 0.1-1\cr
Date: \tab 2013-01-21\cr
Depends: \tab methods, bigmemory\cr
Imports: \tab foreach\cr
LinkingTo: \tab bigmemory\cr
Suggests: \tab MASS, doMC, doMPI, doRedis, doSNOW\cr
License: \tab GPL-3\cr
Built: \tab R 2.15.1; universal-apple-darwin9.8.0; 2013-01-31 04:38:12 UTC; unix\cr
}

Index:

\preformatted{
bigcforest-class        Class '"bigcforest"'
bigctree-class          Class '"bigctree"'
bigrf-package           Big Random Forests-Classification and
                        Regression Forests for Large Data Sets
bigrfc                  Random Forest Classification Model
grow-methods            Methods for Function 'grow' in Package 'bigrf'
merge-methods           Methods for Function 'merge' in Package 'bigrf'
show,bigcforest-method  Methods for Function 'show' in Package 'bigrf'
}

The main entry point for this package is \code{\link[bigrf]{bigrfc}}, which is used to grow a classification random forest based on the given data and forest-growing parameters. \code{bigrfc} returns an object of class \code{"\linkS4class{bigcforest}"} representing the grown forest and containing objects of class \code{"\linkS4class{bigctree}"}.

After a forest is grown, more trees can be grown in the same forest by passing the forest to \code{\link[bigrf]{grow}}.

Multiple forests grown with the same data and parameters can be merged together using \code{\link[bigrf]{merge}}. This is useful, for example, for building forests in parallel on multiple machines, then merging the results into one big forest.

To build trees in parallel, the appropriate parallel backend for \pkg{\link[foreach:foreach-package]{foreach}} must be registered. For example, if you wish to use \pkg{\link[doMC:doMC-package]{doMC}} utilizing all available cores on your machine, run the following before calling \code{bigrfc} or \code{grow} (see \pkg{\link[=foreach-package]{foreach}} for more details):

\preformatted{
library(doMC)
registerDoMC(cores=multicore:::detectCores(all.tests=TRUE))
}

Grown forests can be used to predict the results for new examples using the \code{\link[bigrf]{predict}} method.

Unsupervised learning (as described by Breiman and Cutler) is supported via the \code{\link[bigrf]{generateSyntheticClass}} to generate the second synthetic class. 

The \code{\link[bigrf]{varimp}} and \code{\link[bigrf]{fastimp}} offer two ways to calculate the relative importance of the training variables in growing the forest.

This package performs particularly well for large data sets that cause excessive virtual memory swapping by the OS, which often renders the system unusable for other tasks or unresponsive to user input. Other random forest algorithms (e.g. \code{\link[randomForest]{randomForest}} and \code{\link[party]{cforest}}) may achieve higher performance on smaller data sets.

Currently, only classification trees with limited functionality are implemented. More functions and regression trees will be added in the future.
}
\author{
Original Fortran77 code by Leo Breiman and Adele Cutler.

R port with disk caching and parallelization enhancements by Aloysius Lim.

Maintainer: Aloysius Lim <aloysius.lim@gmail.com>
}
\references{
Breiman, Leo. "Random forests." Machine learning 45, no. 1 (2001): 5-32.

Random Forests web page: \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm}
}
\keyword{package}
\keyword{models}
\keyword{tree}
\seealso{
\code{\link[randomForest]{randomForest}}
\code{\link[party]{cforest}}
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
