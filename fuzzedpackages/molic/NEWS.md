# molic 2.0.0

The package has undergone a major make-over. A slight, but breakable, change in the api of `fit_outlier`. The documentation of `fit_outlier` has been updated and now includes more and better examples of how and when to use the function. The `fit_graph` function is no longer a part of **molic**. It now lives in its own package at [ess](https://github.com/mlindsk/ess) and **molic** is now dependend on **ess**. It is therefore now required to run `include(ess)` to have access to `fit_graph`. 

The readme file has also undergone a major change - the former example using `cars` data has been removed; it was never really a good example showing how to do outlier detection with **molic**.

# molic 1.2.2

 * A new data set, `derma` has been included and a new vignette using this data has been added. 
 * The `tgp_dat` data has now been compressed to save disk space.
 * The `plot.gengraph` function applied to an object (`gengraph`) returned from one of the graph fitting functions (`fit_graph`, `fit_components` etc.) now takes an input that let the user specify the color of the nodes.

# molic 1.2.1

 * `subgraph` function is now provided. 
 * `sapply'`s are now converted to `vapply'`s for safety and potentially more speed when fitting graphs.

# molic 1.2.0

 * `pmf` no longer plots the density of the deviances of a `outlier_model` object. Use `plot` for this instead; this is now consistent with the other related functions like `fit_outlier`. Instead `pmf` is used to construct the probability mass function of a decomposable graphical model which can be used to obtain probabilities of observing specific cells/observations/configurations.

# molic 1.1.0

**Development Model**

From this release we adopt the branching model introduced by Vincent Driessen

 * [Git branching model](https://nvie.com/posts/a-successful-git-branching-model/)

This means, that there are now two branches: the **master branch** is always the current stable version, and the **develop branch** is the develop version.

**New API**

 * Functions like `fit_outlier` that depends on an adjacency list no accept `gengraph` objects returned from `fit_graph` - i.e. no need to use `adj_lst()` first.

**New functions**

 * `generate_multiple_models`
     + Given a class variable with $1,2\ldots, l$ levels and a new observation $y$, this function is a convenient wrapper around `fit_graph` and `fit_outlier` that conducts all the hypothesis $H_k:$ $y$ has level $k$ for $k = 1,2,\ldots, l$.
 * `plot.multiple_models`
     + Given an object returned from `fit_multiple_models` this function is used to visualize all the hypothesis tests for a single observation simultaneously. It is a `ggplot2` object
 * `plot.outlier`
     + Given an object returned from `fit_outlier` this function is used to visualize the approximated density of the deviance under the null hypothesis. It is a `ggplot2` object.
 * `components`
     + Return a list with all components of a graph
 * `fit_components` 

**Misc**
 * All deviances are now non-negative as they should be! Before, a constant was neglected which could potentially confuse the users since a deviance is per definition non-negative.
 
# molic 1.0.0

 * First release.
