## -*- mode: R -*-
###########################################################################
## This file is part of isqg: R package for in silico quantitative genetics
##
##              Copyright (C) 2018 Fernando H. Toledo CIMMYT
##              
## * Filename: ISQG.R
## 
## * Description: Package Documentations
## 
## * Author: Fernando H. Toledo
## 
## * Maintainer: Fernando H. Toledo
## 
## * Created: Fr Mar 09 2018
## 
##   This program is free software; you can redistribute it and/or modify 
##   it under the terms of the GNU General Public License as published by 
##   the Free Software Foundation; either version 2 of the License, or 
##  (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful, but 
##   WITHOUT ANY WARRANTY; without even the implied warranty of 
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
##   General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software Foundation, 
##   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##                                                        
##   `` Far better an approximate answer to the right question, which is 
##   often vague, than the exact answer to the wrong question, which can
##   always be made precise ''
##                          --John Tukey, Ann. Math. Stat. 33(1):13 1962
##
###########################################################################

##' @importFrom Rcpp evalCpp
##' @importFrom R6 R6Class
##' @importFrom Rdpack reprompt
##' @useDynLib isqg
NULL

##' @title isqg: A package to perform \strong{in silico} quantitative genetics
##'
##' @name isqg-package
##'
##' @description \code{isqg} provides R6/C++ classes for in silico quantitative 
##'     genetics. Mimic the meiosis recombination. Allows user defined extensions 
##'     which provides great flexibility. Details of the implementation are described 
##'     in \insertCite{g3;textual}{isqg}.
##'
##' @references
##'     \insertAllCited{}
##' 
##' @author Fernando H. Toledo \email{f.toledo@cgiar.org}
##'
##' @aliases isqg-package isqg
##'
##' @keywords package
##'
##' @rdname ISQG
##'
##' @docType package
NULL

##' @title Toy Example of a Map for \emph{in silico} Quantitative Genetics
##'
##' @name ToyMap
##'
##' @description This data comprise 2 chromossomes with 2cM each and 21 
##'     monitored loci in each chromosome
##'
##' @format a data.frame, 42 rows and 3 collumns (snp, chr, pos).
##'
##' @usage data(ToyMap)
##'
##' @rdname Data
##'
##' @docType data
NULL

## \EOF
###########################################################################
