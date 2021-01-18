## -*- mode: R -*-
###########################################################################
## This file is part of isqg: R package for in silico quantitative genetics
##
##              Copyright (C) 2018 Fernando H. Toledo CIMMYT
##              
## * Filename: Hooks.R
## 
## * Description: Exposes Rcpp Module used by isqg R package
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

## replaces: loadModule(module = "ISQG", what = TRUE)
.onUnload <- function (libpath) {
  gc()
  library.dynam.unload("isqg", libpath)
}

## \EOF
###########################################################################
