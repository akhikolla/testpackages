# Latent Order Logistic (LOLOG) Graph Models

LOLOG is a general framework for generative statistical modeling of graph datasets motivated by
the principle of network growth. This class of models is fully general and terms modeling
different important network features can be mixed and matched to provide a rich generative
description of complex networks.


## Resources

* The mathematical details are outlined in a **[technical paper](http://arxiv.org/abs/1804.04583)**.
* For a more detailed description of what can be done with the ``lolog`` package, **[see the introductory vignette](inst/doc/lolog-introduction.pdf)**.
* An application of LOLOG modeling to a UK Faculty data set with comparisons to an ERGM fit can be found **[here](inst/doc/lolog-ergm.pdf)**.

## Installation

### The Easy Way

To install the latest release from CRAN run:

```
install.packages("lolog")
```

### The Slightly Less Easy Way

To install the latest development version from the github repo run:
```
# If devtools is not installed:
# install.packages("devtools")

devtools::install_github("statnet/lolog")
```
If this is your first R source package that you have installed, you’ll also need a set of development tools. On Windows, download and install [Rtools]( https://cran.r-project.org/bin/windows/Rtools/), and ``devtools`` takes care of the rest. On a Mac, install the [Xcode command line tools]( https://developer.apple.com/downloads). On Linux, install the R development package, usually called ``r-devel`` or ``r-base-dev``. For details see [Package Development Prerequisites](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites).

## Using The Package

```
library(lolog)
library(network)
data(ukFaculty)

# Delete 2 vertices missing group
delete.vertices(ukFaculty, which(is.na(ukFaculty %v% "Group")))

# A dyad independent model
fitind <- lolog(ukFaculty ~ edges() + nodeMatch("GroupC") + nodeCov("GroupC"))
summary(fitind)
```


## Development

[Development Practices and Policies for Contributers](../../wiki/How-to-Contribute:-Git-Practices)

### Using Eclipse

This package is set up as an Eclipse project, and the C++ code can be compiled and run without reinstalling the package. To set up in your eclipse IDE, select import project -> General -> Existing Projects into Workspace and select the lolog directory.

This project was set up following the methods outlined in:

<http://blog.fellstat.com/?p=170>



