## -*- mode: R -*-
###########################################################################
## This file is part of isqg: R package for in silico quantitative genetics
##
##              Copyright (C) 2018 Fernando H. Toledo CIMMYT
##              
## * Filename: Mating.R
## 
## * Description: Simple Mating designs definitions used by isqg R package
## 
## * Author: Fernando H. Toledo
## 
## * Maintainer: Fernando H. Toledo
## 
## * Created: Fr Nov 10 2017
## 
## * Updated: -
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

##' @title Breed Simulated Individuals According to Basic Mating Schemes
##'
##' @name mating
##'
##' @description Performs the simple mating schemes bi-parental cross, 
##'   self-cross and haploid duplication, respectively through the functions 
##'   \code{cross}, \code{selfcross} and \code{dh} and return the respective size 
##'   \emph{n} progeny involving the parental individuals belonging to the 
##'   same specie.
##'
##' @return a size \emph{n} list with instances of the class Specimen that 
##'   represent new individuals belonging to the progeny of the respective mating 
##'   scheme.
##'
##' @param n a length-one integer vector with the size of the progeny.
##' @param p1,p2,gid are instances of the class specimen which will be used 
##'     as the parents.
##'
##' @details Basically this family of functions take simulated individuals 
##'   belonging to the same simulated specie, performs the meiosis that 
##'   generates individual's gametes. According to the scheme applied the 
##'   gametes are merged into new simulated individuals. These are wrap 
##'   functions to the C++ class that mimic the meiosis recombination process.
##' 
##' @examples
##' data(ToyMap)
##' spc <- set_specie(ToyMap)
##' AA <- founder(spc, "AA")
##' aa <- founder(spc, "aa")
##' 
##' ## Mather Design
##' F1  <- cross(n = 1, AA, aa)
##' BC1 <- cross(n = 5, F1, AA)
##' BC2 <- F1$cross(n = 5, aa)   # using R6 methods
##' F2  <- selfcross(n = 10, F1) 
##' RIL <- dh(n = 10, F1) 
##' ## chainable R6 methods
##' F3  <- F1$selfcross(n = 1, replace = TRUE)$selfcross(n = 1, replace = TRUE)
##'
##' @rdname Mating
NULL

##' @export
##' @rdname Mating
'cross' <- function(n = 1, p1, p2) {
  progeny <- .Cpp_cross_ctor(n, p1, p2)
  ## names(progeny) <- paste0(substitute(p1), "_x_", substitute(p2), "_n_", 1:n)
  if(n == 1) progeny <- progeny[[1]]
  return(progeny)
}

##' @export
##' @rdname Mating
'selfcross' <- function(n = 1, gid) {
  progeny <- .Cpp_selfcross_ctor(n, gid)
  ## names(progeny) <- paste0(substitute(gid), "_s_n_", 1:n)
  if(n == 1) progeny <- progeny[[1]]
  return(progeny)
}

##' @export
##' @rdname Mating
'dh' <- function(n = 1, gid) {
  progeny <- .Cpp_doublehaploid_ctor(n, gid)
  ## names(progeny) <- paste0(substitute(gid), "_d_n_", 1:n)
  if(n == 1) progeny <- progeny[[1]]
  return(progeny)
}

## \EOF
###########################################################################
