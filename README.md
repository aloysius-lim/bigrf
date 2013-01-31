bigrf
=====

This is an implementation of Leo Breiman and Adele Cutler's Random Forest algorithms for classification and regression, optimized for data sets that are too large to fit or be processed in memory. It utilizes [`bigmemory`](http://cran.r-project.org/web/packages/bigmemory/) disk-based memory for caching of intermediate computations, and [`foreach`](http://cran.r-project.org/web/packages/foreach/) to parallelize the tree-buliding process. This package performs particularly well for large data sets that cause excessive virtual memory swapping by the OS. Currently, only classification trees with limited functionality are implemented. More functions and regression trees will be added in the future.

Capabilities and Usage
----------------------

The main entry point for this package is `bigrfc`, which is used to grow a classification random forest based on the given data and forest-growing parameters. `bigrfc` returns an object of class `"bigcforest"` representing the grown forest and containing objects of class `"bigctree"`.

After a forest is grown, more trees can be grown in the same forest by passing the forest to `grow`.

Multiple forests grown with the same data and parameters can be merged together using `merge`. This is useful, for example, for building forests in parallel on multiple machines, then merging the results into one big forest.

To build trees in parallel, the appropriate parallel backend for [`foreach`](http://cran.r-project.org/web/packages/foreach/) must be registered. For example, if you wish to use [`doMC`](http://cran.r-project.org/web/packages/doMC/) utilizing all available cores on your machine, run the following before calling `bigrfc` or `grow` (see the documentation for `foreach` for more details):

    library(doMC)
    registerDoMC(cores=multicore:::detectCores(all.tests=TRUE))

This package performs particularly well for large data sets that cause excessive virtual memory swapping by the OS, which often renders the system unusable for other tasks or unresponsive to user input. Other random forest algorithms (e.g. [`randomForest`](http://cran.r-project.org/web/packages/randomForest/) and [`cforest`](http://cran.r-project.org/web/packages/party/)) may achieve higher performance on smaller data sets.

Currently, only classification trees with limited functionality are implemented. More functions and regression trees will be added in the future.

Authors
-------

Original Fortran77 code by Leo Breiman and Adele Cutler.

R port with disk caching and parallelization enhancements by Aloysius Lim.

Licences
--------

`bigrf` is licensed under the [GNU General Public License Version 3 (GPLv3)](http://www.gnu.org/licenses/gpl.html).

License information from original code:

> Copyright 2002-2003  Leo Breiman and Adele Cutler.
> 
> This is free open source software but its use, in part or in whole, in any commercial product that is sold for profit is prohibited without the written consent of Leo Breiman and Adele Cutler.
> 
> We very much appreciate bug notices and suggested improvements.
> 
> leo@stat.berkeley.edu   adele@math.usu.edu

References
----------

Breiman, Leo. "Random forests." *Machine learning* 45, no. 1 (2001): 5-32.

Usage Examples
--------------

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
