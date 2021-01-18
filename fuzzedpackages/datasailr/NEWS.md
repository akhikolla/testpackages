# HISTORY

## Version 0.8.7

* Thanksgiving Day Release (11/26/2020)
* In configure.ac, Onigmo's configure script is changed to be called with --disable-maintainer-mode, by which users do not re-run autoconf. (In some situations where timestamp information is not correctly preserved, autoconf tries to run again on user machines, which may cause aclocal missing error.)
    + To enable/disable maintainer mode, Onigmo's configure.ac needs to have 'AM_MAINTAINER_MODE([enable])'. DataSailr packaging toolchain modifies Onigmo's configure.ac to have this line.
* README and documents now use the term 'DataSailr script'.
* Options for configure script, --enable-libsailr-debug and --enable-datasailr-debug, are now available.
    + R CMD INSTALL can pass these options to configure with '--configure-args=' option.
    + By setting these options, -DDEBUG options are passed to compilers and DEBUG macro variable is defined.
* Fixes
    + When assignment does not happen in script, an error occured. However, there are cases where only discard!() and push!() are executed, and in such cases processed dataset should be returned if fullData = T. This is fixed.


## Version 0.8.6

* Update to follow libsailr's external function availability.
* push!() function is implemented as a DataSailr extended function.
    + It is implemented as an external function of libsailr, because it does more than just row level calculation. Libsailr cannot handle this. At DataSailr level, dataframe structure can be modified, and new rows can be added to result dataframe.
    + The function name, push!, comes from two reasons.
        + It is important to clarify that this function is different from Sailr built-in functions. This function works for more than each row level. It changes dataframe structures, and it has ! at the end.
        + This function does not stop but suspends libsailr vm execution, and adds or pushes rows. Therefore it is named as push.
    + The common usage of this function is converting data of wide format into long format.
* discard!() function is implemented as a DataSailr extended function.
    + discard!() function drops the current row from the result dataframe.
    + A common usage of this function is to filter out specific rows using if condition.


## Version 0.8.5

* Minor update
    + Warning with gcc link time optimization (on Fedora) is resolved.

## Version 0.8.4

* Cross-platform compilation is improved.
* The core processing library, libsailr has become able to report runtime errors with correspoding script location, and DataSailr supports it.
* Fixes
    + Variables that appear only on RHS of script crashed the program, which is fixed and the program stops safely.
    + String objects and regular expression objects were not freed in some cases. Memory leaks are fixed. 
    + NA in StringVector is treated as NULL, which results in passing an empty string to libsailr.


## Version 0.8.3

* Memory leaks fixed at various points that are pointed out after CRAN submission. (Apr.5-13 2020)
* Check whether parsing succeded or not. Function stops if the Sailr script has a wrong syntax. (Apr.14 2020)


## Version 0.8.2

* CRAN release candidate (Mar. 12 2020)
    + Biarch: true in DESCRIPTION.
        + In windows environment, if you have configure.win file, cross-compilation is not executed. This Biarch option solves it.
        + (ref) https://community.rstudio.com/t/configure-win-and-cran-submission/24684/4?fbclid=IwAR3RIMNzABwyBicz6hSPG6m1hmurtEzLGPjgMTFsjXJVkHEqBdW6wYFPd-Q
    + Package C++ main file (data_sailr_cpp_main.cpp) is updated for minor warnings.
        + Rcpp vectors' size() function returns singed int. Chage my code to be compatible with it.
    + Pass R defined makefile variables to Onigmo and libsailr makefiles
        + Makevars.win is updated.

* Minor updates (Mar. 15 2020)
    + Typos are corrected in README.md
    + Also, the following files are updated to resolve CRAN submission problems.
    + library() functions are removed.
        + These are used especially in test functions.
    + Title and description in DESCRIPTION are updated.

* Minor updates (Mar. 18 2020)
    + DESCRIPTION file is updated
    + printf() functions are enclosed within IF_DEBUG() macro in src/data_sailr_cpp_main.cpp
        + Including stdout/stderr related functions is not allowed in CRAN package.
    + bin directory is renamed to exec

* Minor updates (Mar. 22 2020)
    + DESCRIPTION file is updated. (Encoding and SystemRequirements)
    + Example script is added to package document.

* Factor class vector in dataframe is supported. They are now dealt as CharactorVector. (Apr. 1 2020)
    + In this implementation, factor levels is supposed to be CharacterVector.


## Version 0.8.1

* Beta release

## Version 0.8

* First commit
* Thanksgiving Day Release (11/28/2019)
* Christmas Update (12/25/2019) to release to the public
* Updating to submit to CRAN (01/11/2020)
    + The variables updated by this package have character types, not factors.
    + New argument, rowname, is added to sail() function.
    + libsailr is updated.
* Change following the libsailr API update. (01/19/2020)
* cleanup script is updated (2/4/2020)
* Resolve warnings for submitting CRAN (2/5/2020)
* Resolve windows compilation failures. (avoid autotools dependency during installation) (2/8/2020)
* configure.ac updates to cope with autotools nonexistent system. (2/9/2020)
* Makevars.win is updated. To build onigmo, win32/Makefile.mingw is now used. 
    + This leads to successfully build on rhub windows environment. (3/8/2020)


## Version 0.0 (Birth)

* Package skeleton was created. (11/15/2018)
   + The original name was RCppCalc, which was intended only for arithmetic calculation.

