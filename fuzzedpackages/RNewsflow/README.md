RNewsflow: Tools for analyzing content homogeneity and news diffusion using computational text analysis
=======================================================================================================

Given the sheer amount of news sources in the digital age (e.g., newspapers, blogs, social media) it has become difficult to determine where news is first introduced and how it diffuses across sources. RNewsflow provides tools for analyzing content homogeneity and diffusion patterns using computational text analysis. The content of news messages is compared using techniques from the field of information retrieval, similar to plagiarism detection. By using a sliding window approach to only compare messages within a given time distance, many sources can be compared over long periods of time. Furthermore, the package introduces an approach for analyzing the news similarity data as a network, and includes various functions to analyze and visualize this network.

Installation
============

You can install the development version of RNewsflow directly from github:

``` r
library(devtools)
install_github("kasperwelbers/RNewsflow")
```

Vignette
========

The vignette containing a step-by-step tutorial for using RNewsflow can be called from within R.

``` r
library(RNewsflow)
vignette('RNewsflow')
```
