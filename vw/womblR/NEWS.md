# womblR 1.0.4 (2018-12-04)

* The Makevars file was updated. 

# womblR 1.0.3 (2018-06-19)

* The function `PlotSensitivity` has been updated so the plotted values correspond to the user provided legend limits.

* Vignette has been updated to reflect package improvements.

# womblR 1.0.2 (2017-09-12)

* Fixed a bug related to the use of the truncated normal from the `msm` package.

* Updated the vignette.

# womblR 1.0.1 (2017-07-17)

* New `STBDwDM()` functionality, related to the adjacency weights that can be used. A new option is included, `Weights` allows for the `binary` weights of the Lee and Mitchell (2011) specification.

* New `PlotVFTimeSeries()` functionality, related to the location specific regression line (including `line.col`, `line.reg`, and `line.type`). 

* Default bounds for lower and upper bounds for `Phi` are now specific to the temporal correlation structure specification (i.e., if `TemporalStructure` = `"ar1"` the bounds are now appropriate).

* Adjusted the truncated normal random sampler to allow it to use the inverse CDF method or accept-reject method of the `msm` package when necessary. 

* Tidying of the help pages.

# womblR 1.0.0 (2017-06-15)

* First release.
