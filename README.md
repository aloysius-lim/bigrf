bigrf
=====

This is an R implementation of Leo Breiman's and Adele Cutler's Random Forest algorithms for classification and regression, with optimizations for performance and for handling of data sets that are too large to be processed in memory. Forests can be built in parallel at two levels. First, trees can be grown in parallel on a single machine using [`foreach`](http://cran.r-project.org/web/packages/foreach/). Second, multiple forests can be built in parallel on multiple machines, then merged into one. For large data sets, disk-based [`big.matrix`](http://cran.r-project.org/web/packages/bigmemory/)'s may be used for storing data and intermediate computations, to prevent excessive virtual memory swapping by the operating system. Currently, only classification forests with a subset of the functionality in Breiman and Cutler's original code are implemented. More functionality and regression trees may be added in the future.

Capabilities and Usage
----------------------

  The main entry point for this package is `bigrfc`, which is used to build a classification random forest on the given training data and forest-building parameters. `bigrfc` returns the forest as an object of class `"bigcforest"`, which contains the trees grown as objects of class `"bigctree"`. After a forest is built, more trees can be grown using `grow`.

Performance Optimizations
-------------------------

For better performance, trees may be grown in parallel by registering an appropriate parallel backend for [`foreach`](http://cran.r-project.org/web/packages). As an example, the following code uses the [`doParallel`](http://cran.r-project.org/web/packages/doParallel/) backend to enable tree-growing on all available cores on the machine. This code must be executed before calling `bigrfc` or `grow`. See the documentation for [`foreach`](http://cran.r-project.org/web/packages) for more details on supported parallel backends.

    library(doParallel)
    registerDoParallel(cores=detectCores(all.tests=TRUE))

Multiple random forests can also be built in parallel on multiple machines (using the same training data and parameters), then merged into one forest using `merge`.

For large data sets, the training data, intermediate computations and some outputs (e.g. proximity matrices) may be cached on disk using [`"big.matrix"`](http://cran.r-project.org/web/packages/bigmemory/) objects. This enables random forests to be built on fairly large data sets without hitting RAM limits, which will cause excessive virtual memory swapping by the operating system.

Disk caching may be turned off for optimal performance on smaller data sets by setting function / method argument `cachepath` to `NULL`, causing the `big.matrix`'s to be created in memory.

Authors
-------

Original Fortran77 code by Leo Breiman and Adele Cutler.
  
R port with disk caching and parallelization enhancements by Aloysius Lim.

Licences
--------

`bigrf` is licensed under the [GNU General Public License Version 3 (GPLv3)](http://www.gnu.org/licenses/gpl.html).

References
----------

Breiman, L. (2001). Random forests. *Machine learning*, 45(1), 5-32.
  
Breiman, L. & Cutler, A. (n.d.). Random Forests. Retrieved from [http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm](http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm).

Usage Examples
--------------
    # Classify cars in the Cars93 data set by type (Compact, Large, Midsize,
    # Small, Sporty, or Van).
    
    # Load data.
    data(Cars93, package="MASS")
    x <- Cars93
    y <- Cars93$Type
    
    # Select variables with which to train model.
    vars <- c(4:22)
    
    # Run model, grow 50 trees on the first 60 examples.
    forest1 <- bigrfc(x[1:60, ], y[1:60], ntree=50L, varselect=vars)
    
    # Build a second forest.
    forest2 <- bigrfc(x[1:60, ], y[1:60], ntree=50L, varselect=vars)
    
    # Merge the 2 forests.
    forest <- merge(forest1, forest2)
	  
    # Grow even more trees.
    forest <- grow(forest, x[1:60, ], ntree=30L)
	  
    # Get predictions for the remaining examples.
    predictions <- predict(forest, x[-(1:60), ], y[-(1:60)])
    
    # Calculate variable importance, including those for each out-of-bag
    # example.
    importance <- varimp(forest, x[1:60, ], impbyexample=TRUE)
    
    # Calculate variable importance using the fast (Gini) method.
    fastimportance <- fastimp(forest)
    
    # Calculate variable interactions.
    inter <- interactions(forest)
    
    # Calculate proximity matrix and first two scaling co-ordinates.
    prox <- proximities(forest)
    scale <- scaling(prox)
    
    # Plot the 1st vs 2nd scaling co-ordinates.
    plot(scale, col=as.integer(y) + 2, pch=as.integer(y) + 2)

    # Calculate outlier scores, and circle the top 20% percent of them in red.
    outscores <- outliers(forest)
    points(scale[outscores > quantile(outscores, probs=0.8), ], col=2, pch=1,
           cex=1.5)
    
    # Compute prototypes.
    prot <- prototypes(forest, prox, x=x[1:60, ])
    
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
    
