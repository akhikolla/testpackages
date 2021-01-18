---
title: "About"
author: "Anders Ellern Bilgrau"
date: "2019-03-25"
output: html_document
---



## About

This is the [**shiny**](https://shiny.rstudio.com/) web-application accompanying the [**GMCM**](https://cran.r-project.com/package=GMCM) package[1] for [**R**](https://r-project.com). The application is developed and maintained by [Anders Ellern Bilgrau](https://bilgrau.com). You can read more about the underlying methods employed for unsupervised clustering [here](http://doi.org/10.18637/jss.v070.i02) and [here](https://github.com/AEBilgrau/GMCM#GMCM) and the references therein.

### Availability 
The application is available online at 

> https://gmcm.shinyapps.io 

The application can also be available as a local instance on your own machine via the **GMCM** package. To run the application locally, make sure **GMCM** is installed and run


```r
GMCM::runGMCM()
```

to start a local shiny server. More information on installing **GMCM** and the prerequisites is found at [CRAN](https://cran.r-project.org/web/packages/GMCM/readme/README.html) or [the repository at GitHub](https://github.com/AEBilgrau/GMCM#installation).

This instance is running **GMCM version 1.3.1**.


### Citation
If you use the web-application or **GMCM** package itself, please cite the reference below.

1. Bilgrau AE, Eriksen PS, Rasmussen JG, Johnsen HE, Dybkaer K,
Boegsted M (2016). "GMCM: Unsupervised Clustering and
Meta-Analysis Using Gaussian Mixture Copula Models." _Journal of
Statistical Software_, *70*(2), 1-23. doi: 10.18637/jss.v070.i02
(URL: http://doi.org/10.18637/jss.v070.i02).

In **R** you can run `citation("GMCM")` to get a copy-friendly bibtex format.

---


<div align="center">
  Anders Ellern Bilgrau
  &mdash;
  <a href="https://github.com/AEBilgrau/GMCM" padding="10px">
    <i class="fab fa-github"></i>
  </a>
  <a href="https://stackoverflow.com/users/1568306/anders-ellern-bilgrau">
    <i class="fab fa-stack-overflow"></i>
  </a>
  <a href="https://scholar.google.dk/citations?user=zQNl61YAAAAJ&amp;hl=en">
    <i class="fab fa-google"></i>
  </a>
  <a href="https://www.linkedin.com/in/aebilgrau/">
    <i class="fab fa-linkedin"></i>
  </a>
  &mdash;
  2019
</div>

