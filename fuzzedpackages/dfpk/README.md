 [![CRAN Version](https://www.r-pkg.org/badges/version/dfpk)](https://cran.r-project.org/package=dfpk)
 ![](https://cranlogs.r-pkg.org/badges/grand-total/dfpk)
  
# dfpk

The **dfpk** R package provides an interface to fit Bayesian generalized (non-)linear mixed models using Stan, which is a C++ package for obtaining Bayesian inference using the No-U-turn sampler (see http://mc-stan.org/). 

### Description

dfpk package includes methods involving PK measures in the dose allocation process during a Phase I clinical trials. These methods enter PK in the dose finding designs in different ways, including covariates models, dependent variable or hierarchical models. This package provides functions to generate scenarios, and to run simulations which their objective is to determine the maximum tolerated dose (MTD).

#### Installation 

### Establish Version  

A latest version of the package **dfpk** is available on CRAN and can be loaded via 

```{r} 
install.packages("dfpk")
library(dfpk) 
```  

### Development Version 
To install the **dfpk** package from GitHub, first make sure that you can install the **rstan** package and C++ toolchain by following these [instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). The program Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) comes with a C++ compiler for Windows. On OS-X, you should install Xcode. Once **rstan** is successfully installed, you can install **dfpk** from GitHub using the **devtools** package by executing the following in R:

```{r}
if (!require(devtools)){
  install.packages("devtools") 
  library(devtools) 
}

install_github("artemis-toumazi/dfpk")
```

If installation fails, please let us know by [filing an issue](https://github.com/artemis-toumazi/dfpk/issues). 

Details on formula syntax, families and link functions, as well as prior distributions can be found on the help page of the dfpk function:
```{r help.dfpk, eval=FALSE}
help(dfpk) 
```

#### FAQ

### Can I avoid compiling models? 

Unfortunately, fitting your model with **dfpk**, there is currently no way to avoid the compilation. 

### What is the best way to ask a question or propose a new feature? 

You can either open an issue on [github](https://github.com/artemis-toumazi/dfpk) or write me an email to (artemis.toumazi@gmail.com). 
