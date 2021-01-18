# intervalaverage 0.8.0 (Release)

1. No major changes from 0.0.0.9010 other than changing the readme to reflect that this may be soon downloadable from CRAN.



# intervalaverage 0.0.0.9010 (Development)

1. Slight reversal to the breaking change before: I had forgotten I was already allowing different group_vars in x and y via providing named group_vars for intervalintersect so this functionality is turned back on. Naming either interval_vars or group_vars in the intervalintersect function is still explicitly disallowed since the functionality of allowing different join column names in x and y for that function is still not implemented.

2. Rewrote some internals of the intervalintersect function and added a few additional tests. This, combined with the major rewrite of the intervalaverage function in 9009, should hopefully resolve the memory issues that were occurring i386 windows during R CMD CHECK.

# intervalaverage 0.0.0.9009 (Development)

1. POTENTIALLY BREAKING CHANGE: interval_vars and group_vars arguments cannot have names for the intervalaverage function. This is to allow for a future release where x and y could have differently named interval and/or grouping columns. Additionally, in the intervalintersect function group_vars cannot have names for the same reason. intervalintersect already allows for different interval_vars names in x and y via this exact scheme, so named interval_vars arguments to intervalintersect are obviously still allowed.

2. Turned off some tests when building on CRAN that are likely never going to pass on i386 due to memory issues with the slower non-exported testing/validation function (interval_average_slow_f).


# intervalaverage 0.0.0.9008 (Development)

1. Optimized interval average for memory usage.


# intervalaverage 0.0.0.9007 (Development)

1. Performance check is no longer tested on cran.

# intervalaverage 0.0.0.9006 (Development)

1. A new unit test/regression test has been added to ensure that the function is using the optimal algorithm based on data.table performance.

# intervalaverage 0.0.0.9005 (Development)

1. fixed a bug where the y argument of intervalaverage was not being restored to its original state 
(row order and removing ordering column rowindex) if there were overlapping intervals in y.

# intervalaverage 0.0.0.9004 (Development)

1. More improvements to the vignette. The final example is now complete 
(previously it did not actually average the clipped exposures).

# intervalaverage 0.0.0.9003 (Development)

1. More minor improvements in documentation (edits to intro vignette and help files)
2. Improvements in code to avoid side-effects (ie changes in row order) that were originally implemented in 0.0.0.9001, based on suggestions from https://github.com/Rdatatable/data.table/issues/4575#issuecomment-650656185  
3. intervalintersect now returns columns in a more sensible order.

# intervalaverage 0.0.0.9002 (Development)

## BUG FIXES
1 . Minor improvements in documentation (readme,intro vignette)

# intervalaverage 0.0.0.9001 (Development)

## POTENTIALLY BREAKING CHANGES

1. Interval columns used by `intervalaverage`, `intervalintersect`,
and `isolateoverlaps` functions (i.e. those columns of `x` and `y`
specified by the `interval_vars` argument in all of those functions)
must now be of either class `integer`
or `IDate`. This is a change from previous where columns of class `Date` were allowed. 
`Date` objects could have hidden decimals which might cause problems in this function
which is written to deal with only discrete time. Changing Dates to IDates should fix existing code.


2. `intervalaverage` now restores the keys of x and y to their state prior to the function call if x and y are not already keyed on c(group_vars,interval_vars) and prints a message recommending that the keys be set first (to avoid unsetting the key at the end of the function if you don't care about the order). This is a breaking change only if code relied on `intervalaverage` changing the key.




## NEW FEATURES

## BUG FIXES

1. `isolateoverlaps` no longer sets the key of x internal to the function as this was never necessary for the foverlaps call to work. This should reduce the possibility of side-effects for this function even though x is still passed by reference. This fixes the problem that the key of x might changed if the function returns an error (the function was already written to restore the key of x if it successfully completed).

2. CJ.dt no longer passes-by-references the data.tables provided as arguments in `...`; While this may add some more memory usage and processing time to copy the table, it really shouldn't matter since if you don't have memory to copy a table you don't have memory to cartesian join it with another table.

3. `is.overlapping` now restores the key and column order of x to its original state even if the function returns an error.
