gastempt: Fitting gastric emptying curves
====================================

Dieter Menne  
Menne Biomed Consulting Tübingen, Germany  
http://www.menne-biomed.de  
dieter.menne@menne-biomed.de   



[![](https://travis-ci.org/dmenne/gastempt.svg?branch=master)](https://travis-ci.org/dmenne/gastempt) [![](https://coveralls.io/repos/github/dmenne/gastempt/badge.svg?branch=master)](https://coveralls.io/github/dmenne/gastempt?branch=master) [![](https://cranlogs.r-pkg.org/badges/grand-total/gastempt)](https://cran.r-project.org/package=gastempt) [![](https://www.r-pkg.org/badges/last-release/gastempt)](https://CRAN.R-project.org/package=gastempt)


A package and a [Shiny](https://shiny.rstudio.com/) web application to create simulated gastric emptying data, and to analyze gastric emptying from clinical studies using a population fit with R and package `nlme`. In addition,Bayesian fits with [Stan](http://mc-stan.org/) to handle critical cases are implemented.

Part of the work has been supported by section GI MRT, Klinik für Gastroenterologie und Hepatologie, Universitätsspital Zürich; thanks to Werner Schwizer and Andreas Steingötter for their contributions.

### Download
The package is available from [CRAN](https://CRAN.R-project.org/package=gastempt) and github ([source](https://github.com/dmenne/gastempt), [documentation](https://dmenne.github.io/gastempt/)). To install, use:

```
devtools::install_github("dmenne/gastempt")
```

Compilation of the Stan models needs several minutes. 

### Shiny online interface

The web interface can be installed on your computer, or run as [web app](  
https://apps.menne-biomed.de/gastempt/).

Two __models__ are implemented in the web interface

* `linexp, vol = v0 * (1 + kappa * t / tempt) * exp(-t / tempt):`Recommended for gastric emptying curves with an initial volume overshoot from secretion. With parameter kappa > 1, there is a maximum after t=0.  When all emptying curves start with a steep drop, this model can be difficult to fit.
* `powexp, vol = v0 * exp(-(t / tempt) ^ beta):` The power exponential function introduced by Elashof et. al. to fit scintigraphic emptying data; this type of data does not have an initial overshoot by design. Compared to the `linexp` model, fitting `powexp` is more reliable and rarely fails to converge in the presence of noise and outliers. The power exponential can be useful with MRI data when there is an unusual late phase in emptying.

### Methods with variants

* Population fits based on function `nlme` in package R `nlme`.
* [Stan-based fits](http://menne-biomed.de/blog/tag:Stan), both [without](https://github.com/dmenne/gastempt/blob/master/inst/stan/linexp_gastro_1b.stan) and [with](https://github.com/dmenne/gastempt/blob/master/inst/stan/linexp_gastro_2b.stan) covariance estimation. Thanks to priors, fitting with Bayesian methods also works for single records, even if stability strongly improves with more data sets available; see  [stan_gastempt](https://dmenne.github.io/gastempt/reference/stan_gastempt.html). Some details can be found in [my blog](https://menne-biomed.de/blog/ballot-and-bazaar). The rationale for using Stan to fit non-linear curves is discussed [here](http://menne-biomed.de/blog/breath-test-stan) for <sup>13</sup>C breath test data, but is equally valid for gastric emptying data. 

### Data entry:
* Data can be entered directly from the clipboard copied from Excel, or can be simulated using a Shiny app.
* Several preset simulations are provided in the Shiny app, with different amounts of noise and outliers 
* Robustness of models can be tested by manipulating noise quality and between-subject variance. 
* Fits are displayed with data.
* The coefficients of the analysis including t50 and the slope in t50 can be downloaded in .csv format.


### Example 

Program with simulated data (needs about 40 seconds till plot shows):

```
library(gastempt)
dd = simulate_gastempt(n_records = 6, seed = 471)
d = dd$data
ret = stan_gastempt(d)
print(ret$coef)
print(ret$plot)
```

![Screenshot](tools/readme/screenshot.png)


## Docker image

The image cannot be compiled on the Docker hub because the build runs out of memory in the standard configuration.

### Installing Docker 
- For Windows 10, you can get the installer from the [Docker store](https://store.docker.com/editions/community/docker-ce-desktop-windows). For installation details, see [here](https://docs.docker.com/docker-for-windows/install/).  
- Linux users know how to install Docker anyway. 
- Docker should have at least 2 GB of memory; on Windows, use Settings from the Docker tray icon. If you want to build the Docker image, you need at least 4 GB and 2 cores; confusing error messages are being emitted when memory is low.

### Installing `gastempt` 

- From the command line, enter the following to start the container

```
docker run --name gastempt  --restart unless-stopped -p 3838:3838 -d dmenne/gastempt
```
- The first startup needs some time because 1 GB has to be downloaded. Subsequent startups require only a few seconds.
- Connect to the app with your browser via `localhost:3838`.





