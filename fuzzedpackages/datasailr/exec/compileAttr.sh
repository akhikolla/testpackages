#!/bin/sh
#
# Note
# This script is discussed in https://github.com/RcppCore/Rcpp/issues/136
# To preare src/RcppExports.cpp and R/RcppExports.R, you need to run Rcpp::compileAttributes at the package root directory, and this script does this work.
# In this package, this script is called from configure.ac and cleanup script for the following reasons.
# 
# * configure.ac
#   + Whether configure script is executed or not depends on the derived files from configure script. Sometimes configure script does not run and RcppExports files are not created.
#   + I need to run autoconf to reflect configure.ac change to configure.
#   + 
# * cleanup script
#   + This alsway runs as long as I use option --preclean for R CMD INSTALL.
#   +
# 
#
#Rscript -e 'library(Rcpp); compileAttributes(".")'
#
