
<!-- README.md is generated from README.Rmd. Please edit that file -->

# uFTIR

**THE BRANCH hcl HAS MODIFICATIONS TO ADAPT THE PACKAGE TO THE
HIGH-PERFORMANCE CLUSTER**

An R package to read and process Agilent Cary 620 FTIR microscope
-hyperspectral- images. As it is, the package performs image
classification to look for microplastics in environmental samples.

When I wrote the package I took some ideas from:

1.  The implementation of the Spectal Angle Mapper algorithm in
    RcppArmadillo was copied from
    [RStoolbox](https://bleutner.github.io/RStoolbox/). I modified the
    function to do the two for-loops in C++, instead of only one as in
    Leutner, B. et al.’s version.
2.  I’m using the [R raster
    package](https://cran.r-project.org/package=raster) to plot. This is
    not efficient -I need to coerce every object to a RasterLayer-, and
    it makes the package huge -I need to import the raster methods-. But
    it is enough for now.
3.  The spectral library included as reference comes from [Primpke et al
    (2018)](https://doi.org/10.1007/s00216-018-1156-x).
4.  The function read\_tile is based on the Octave/Matlab code of [Alex
    Henderson](https://bitbucket.org/AlexHenderson/agilent-file-formats/src/master/),
    published [here](https://doi.org/10.5281/zenodo.399238). He solved
    the Agilent-provides-binary-files problem.

## Installation

Install the current development from github via:

``` r
remotes::install_github("fcorra/uFTIR")
```
