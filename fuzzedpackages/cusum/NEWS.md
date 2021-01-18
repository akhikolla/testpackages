# cusum 0.4.1 (2019/10/02)

### Bug fixes
* improved passing of custom weights
* `alpha` calculation works correctly for process improvements


# cusum 0.4.0 (2019/08/28)

### New features
* new function to calculate group-sequential CUSUM (GSCUSUM) for risk-adjusted and non-risk-adjusted processes
* vignette describing gscusum-functions
* custom weights can be passed to cusum and racusum functions

# cusum 0.3.0 (2019/06/17)

### New features
* new function to calculate CUSUM control limits exactly for very small sample sizes
* function to calculate DPCL following Zhang and Woodall (2016) for RA-CUSUM
* detect process improvements with cusum() and racusum()

# cusum 0.2.1 (2019/03/19)

### Bug fixes
* fixed typo in cusum_limit_sim and cusum_alpha_sim


# cusum 0.2.0 (2019/03/15)

### New features
* S3 method for `plot.cusum()`
* `cusum` class
