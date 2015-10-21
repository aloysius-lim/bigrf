\name{bigrf-package}
\alias{bigrf-package}
\alias{bigrf}
\docType{package}
\title{Big Random Forests: Classification and Regression Forests for Large Data Sets}
\description{
  This is an implementation of Leo Breiman's and Adele Cutler's Random Forest algorithms for classification and regression, with optimizations for performance and for handling of data sets that are too large to be processed in memory. Forests can be built in parallel at two levels. First, trees can be built in parallel on a single machine using \code{\link[foreach]{foreach}}. Second, multiple forests can be built in parallel on multiple machines, then merged into one. For large data sets, disk-based \code{\link[bigmemory]{big.matrix}}'s may be used for storing data and intermediate computations, to prevent excessive virtual memory swapping by the operating system. Currently, only classification forests with a subset of the functionality in Breiman and Cutler's original code are implemented. More functionality and regression trees will be added in the future. See file INSTALL-WINDOWS in the source package for Windows installation instructions.
}
\details{
  \tabular{ll}{
    Package: \tab bigrf\cr
    Version: \tab 0.1-12\cr
    Date: \tab 2015-10-21\cr
    OS_type: \tab unix\cr
    Depends: \tab R (>= 2.14), methods, bigmemory (>= 4.5.8)\cr
    Imports: \tab foreach\cr
    Suggests: \tab MASS, doParallel\cr
    LinkingTo: \tab bigmemory, BH\cr
    License: \tab GPL-3\cr
    Copyright: \tab 2013-2015 Aloysius Lim\cr
    URL: \tab \url{https://github.com/aloysius-lim/bigrf}\cr
    BugReports: \tab \url{https://github.com/aloysius-lim/bigrf/issues}\cr
  }

  Index:
  \preformatted{
    bigcforest-class        Classification Random Forests
    bigcprediction-class    Random Forest Predictions
    bigctree-class          Classification Trees in Random
                            Forests
    bigrf-package           Big Random Forests: Classification
                            and Regression Forests for Large Data
                            Sets
    bigrfc                  Build a Classification Random Forest
                            Model
    bigrfprox-class         Proximity Matrices
    fastimp-methods         Compute Fast (Gini) Variable
                            Importance
    generateSyntheticClass  Generate Synthetic Second Class for
                            Unsupervised Learning
    grow-methods            Grow More Trees in a Random Forest
    interactions-methods    Compute Variable Interactions
    merge-methods           Merge Two Random Forests
    outliers-methods        Compute Outlier Scores
    predict-methods         Predict Classes of Test Examples
    proximities-methods     Compute Proximity Matrix
    scaling-methods         Compute Metric Scaling Co-ordinates
    varimp-methods          Compute Variable Importance
  }

  The main entry point for this package is \code{\link[bigrf]{bigrfc}}, which is used to build a classification random forest on the given training data and forest-building parameters. \code{bigrfc} returns the forest as an object of class \code{"\linkS4class{bigcforest}"}, which contains the trees grown as objects of class \code{"\linkS4class{bigctree}"}. After a forest is built, more trees can be grown using \code{\link[bigrf]{grow}}.
}
\section{Performance Optimizations}{
  For better performance, trees may be grown in parallel by registering an appropriate parallel backend for \pkg{\link[foreach:foreach-package]{foreach}}. As an example, the following code uses the \pkg{\link[doParallel:doParallel-package]{doParallel}} package to enable tree-growing on all available cores on the machine. This code must be executed before calling \code{bigrfc} or \code{grow}. See \pkg{\link[foreach:foreach-package]{foreach}} for more details on supported parallel backends.
  \preformatted{
    library(doParallel)
    registerDoParallel(cores=detectCores(all.tests=TRUE))
  }
  Multiple random forests can also be built in parallel on multiple machines (using the same training data and parameters), then merged into one forest using \code{\link[bigrf]{merge}}.
  \cr\cr
  For large data sets, the training data, intermediate computations and some outputs (e.g. proximity matrices) may be cached on disk using \code{\link[bigmemory]{big.matrix}} objects. This enables random forests to be built on fairly large data sets without hitting RAM limits, which will cause excessive virtual memory swapping by the operating system.
  \cr\cr
  Disk caching may be turned off for optimal performance on smaller data sets by setting function / method argument \code{cachepath} to \code{NULL}, causing the \code{big.matrix}'s to be created in memory.
}
\author{
  Original Fortran77 code by Leo Breiman and Adele Cutler.

  R port with disk caching and parallelization enhancements by Aloysius Lim.

  Maintainer: Aloysius Lim
}
\references{
  Breiman, L. (2001). Random forests. \emph{Machine learning}, \bold{45(1)}, 5-32.

  Breiman, L. & Cutler, A. (n.d.). Random Forests. Retrieved from \url{http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm}.
}
\seealso{
  \code{\link[randomForest]{randomForest}}
  \code{\link[party]{cforest}}
}
\keyword{package}
\keyword{models}
\keyword{tree}
