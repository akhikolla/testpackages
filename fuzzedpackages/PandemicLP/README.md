# Installing PandemicLP

First, make sure you have the `devtools` package installed by typing 

```R
install.packages("devtools")
```

Then run code

```R
devtools::install_github("CovidLP/PandemicLP", build_vignettes = TRUE)
```

The option `build_vignettes = TRUE` will make the installation take longer, but it will make the example vignette available for viewing.

## If installation fails on Mac OS

Installation on Mac requires 3.6.2, so make sure it is up to date.

If you are using the Catalina version of the Mac Operating System, go to [here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-a-Mac) and follow the same instructions for installing rstan from source.

## If installation fails on Windows

In RStudio (preferably) or otherwise in R, execute once

pkgbuild::has_build_tools(debug = TRUE)

to check your C++ toolchain using the pkgbuild package that gets installed when you install RStan.

Then,

1. If you are using Rstudio, a pop-up will appear asking if you want to install Rtools, if you do not already have it. Click Yes and wait about a minute until the installation is finished.
2. If you are not using RStudio, a message will appear in the R console telling you to install Rtools.

In some cases, Windows requires this block of code to be run in R for the installation to be successful. However, do not do this unless the installation has failed.

```R
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
if (!file.exists(M)) file.create(M)
cat(if( grepl("^darwin", R.version$os)) "\nCXX14FLAGS=-O3 -march=native -mtune=native -arch x86_64 -ftemplate-depth-256" else 
    if (.Platform$OS.type == "windows") "\nCXX14FLAGS=-O3 -mtune=native -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2" else
    "CXX14FLAGS = -fPIC",
    file = M, sep = "\n", append = TRUE)
```

## If installation still fails

If the steps above were tried for your OS and it still didn't work, please e-mail us at covidlp.team@gmail.com
