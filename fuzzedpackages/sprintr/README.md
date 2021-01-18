# sprintr

The `sprintr` package contains implementation of a computationally efficient method to fit large-scale interaction models based on a reluctant interaction selection principle.
The details of the method can be found in 
[Yu, Bien, and Tibshirani (2019) *Reluctant Interaction Modeling*](https://arxiv.org/abs/1907.08414).

To install `sprintr` from [github](http://github.com), type in R console
```R
devtools::install_github("hugogogo/sprintr", build_vignettes = TRUE)
```
Note that the installation above requires using R package [devtools](https://CRAN.R-project.org/package=devtools)
(which can be installed using `install.packages("devtools")`).

Please check the accompanying vignette on how to use the `sprintr` package. To read vignette, after installing the package, type in R console
```R
browseVignettes("sprintr")
```

