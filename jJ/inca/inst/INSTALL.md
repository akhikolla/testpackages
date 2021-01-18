# Installation requirements and instructions

## Requirements
Most recent versions of the following software are always preferred, however the minimal requirements are also specified. The common software required by three main-stream operative systems are provided below, and any specific requirement is treated separately:

 * [R software](http://www.r-project.org/): the minimal requirement is R version 3.5.1.

 * [R packages](http://cran.r-project.org/) (it is recommended if the most recent versions of the following packages are installed):
     * `Rcpp`, for interfacing C++ code with R;
     * `RcppArmadillo`, for interfacing C++ algorithms for linear algebra with R;
     * `devtools`, for developing utilities and interfaces with on-line repositories;

To install the required packages, the following R code should be exectued before the installation:
```R
pkgnames <- c("Rcpp", "RcppArmadillo", "devtools")
pkgnames <- pkgnames[!pkgnames %in% .packages(TRUE)]
if (length(pkgnames)) install.packages(pkgnames)
```
To update all the installed packages to the last version available, the following command line should be typed into an R console:
```R
update.packages(ask = FALSE)
```
The additional software must be installed before to update the required packages.

### Additional software on Linux

 * [GNU compilers](https://gcc.gnu.org/): the minimal requirement is gcc version 4.9.3.

### Additional software on Windows

 * [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is a library which provides the GNU compilers. The minimal requirement is `Rtools33.exe`, or other versions according to the R version already installed. Successively, the environmental variable `PATH` must be edited, because it must include the folders containing the GNU compilers provided.

### Additional software on (Mac) OS X

 * [X-code](https://developer.apple.com/xcode/download/) is required to install the GNU compilers. The minimal requirement is gcc version 4.9.3, which allows for C++11 syntax.

## Package management
### Installing the stable release of the inca package

The use of the following R command is highly suggested to install the **inca** package:
```R
install.packages("inca")
```

The other alternative to install an R package is from its source-code compressed as a tarball archive. This can be done by entering the following command into a terminal session on **Linux** and **(Mac) OS X**  :
```bash
R CMD INSTALL inca_0.0.4.tar.gz
```

On **Windows**, by opening the command prompt (`cmd.exe`), it is possible to point to the proper directory with `cd`, and then install the package via `Rcmd.exe` with the following command:

```bash
Rcmd.exe INSTALL inca_0.0.4.tar.gz
```

More details can be found on the "Installing packages" section of the [R-admin](https://cran.r-project.org/doc/manuals/R-admin.html) manual.

### Installing the current development version of the package
```R
devtools::install_github("drwolf85/inca")
```

### Updating the inca package
To update the **inca** package, it is necessary to type the following code from the R console:
```R
update.packages("inca")
```

### Removing the inca package
To remove **inca** from the list of R packages, it is necessary to type the following code from the R console:
```R
remove.packages("inca")
```
