# rocTree 1.1.1
  * Resubmission; previous version was archived 
# rocTree 1.1.0
  * Resolve CRAN error: "calling math functions such as log pow sqrt with integer arguments is not portable."
  * Rewrite tree/ensemble algorithms with Rcpp
	* Significantly speed improvement

# rocTree Version 1.0.0
	* Update DESCRIPTION as suggested by CRAN

# rocTree 0.99-9
	* Modified dCON
	* Added forest
	* Slight speed improvement.
	* Rebuild with roxygen2

# rocTree 0.99-5
	* Modified the print funciton; tree now print with data.tree package, but the old way can still be called.
	* Added the plot function (from DiagrammeR); tree can be plotted in seveal ways.
	* Added plot saving feature.

# rocTree 0.99-4
	* Added the print function; tree can be printed in two ways determined by the logical value of top2bottom

# rocTree 0.99-3
	* Fixed folding in cross validation; now the folding are generated based on percentile as in caret::createFolds

# rocTree 0.99-2
	* Added plot.rocTree to plot hazard functions for the final terminal nodes.
	* Added and cleaned dfFinal and ndFinal.
	* Added an additional smoothing parameter for plots.

# rocTree 0.99-1
	* Added parallel computing for VC.
	* Added rocTree.prune.
	* Drafted the help pages.
	* Improved speed in data preparation.
