## -*- mode: R -*-
###########################################################################
## This file is part of isqg: R package for in silico quantitative genetics
##
##              Copyright (C) 2018 Fernando H. Toledo CIMMYT
##              
## * Filename: Trait.R
## 
## * Description: Defines Trait functionalities (standalone methods)
## 
## * Author: Fernando H. Toledo
## 
## * Maintainer: Fernando H. Toledo
## 
## * Created: Tu Apr 10 2018
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

##' @title Simulated Trait for Individuals According to Basic Models
##'
##' @name fitness
##'
##' @description Constructor of instances of the Trait class given the focal 
##'     specie and the parameters to define infinitesimal or quantitative fitness.
##'
##' @return Objects of R6 class with methods to mimic in silico Traits.
##'
##' @param specie an instance of the R6 class Specie with the genome's parameters.
##' @param m,a,d  a length-one numeric vector with respectively the mean, the 
##'     additive and the dominant effects.
##' @param genes  a character vector with the putative genes.
##' @param data   a data frame with the genes (snp) and their additive and dominant 
##'     effects
##'
##' @details Infinitesimal traits need the mean, the additive and the dominant 
##'     effect and optionally the vector of the putative genes. Quantitative 
##'     traits are defined given the mean and a data frame with the putative 
##'     genes and their additive and dominant effects.
##' 
##' @examples
##' data(ToyMap)
##' spc <- set_specie(ToyMap)
##' AA <- founder(spc, "AA")
##' aa <- spc$founder("aa")
##' 
##' F1  <- cross(n = 1, AA, aa) # the hybrid
##' 
##' ## set a infinitesimal & a quantitative fitness
##' infty <- set_infty(spc, m = 0, a = 1, d = .5) # partial dominance
##' genes <- data.frame(snp = sample(ToyMap$snp, 10), add = rnorm(10), dom = rnorm(10))
##' quant <- set_quant(spc, m = 0, data = genes)
##' 
##' ## evaluating the breeding value 
##' infty$alpha(AA)
##' quant$alpha(F1)
##'
##' @rdname fitness
NULL

##' @export
##' @rdname fitness
"set_infty" <- function(specie, m = 0, a = 1, d = 0, genes = NULL) {
  if (all(c(is.null(a), is.null(d))))
    stop( "Trait without genetic variation" )
  if (all(c(a == 0, d == 0)))
    stop( "Trait without genetic variation" )
  map <- specie$map()
  if (is.null(genes))
    genes <- map$snp
  if(!all(genes %in% map$snp))
    stop( "Genes were wrongly specified" )
  data_joint <- merge(map, data.frame(snp = genes, add = a, dom = d), 
                      by = 'snp', all = TRUE)
  data_order <- data_joint[order(data_joint$chr, data_joint$pos),]
  loci <- with(data_order, ifelse(is.na(add) & is.na(dom), 0, 1))
  vec_loci <- unname(sapply(split(loci, data_order$chr),
                            function(x) paste(x, collapse = '')))
  return(.Cpp_trait_infty_ctor(specie, vec_loci, m, a, d))
}

##' @export
##' @rdname fitness
"set_quant" <- function(specie, m, data) {
  if (!all(c('snp', 'add', 'dom') %in% names(data)))
    stop( "data must have 'snp', 'add' and 'dom' collumns" )
  map <- specie$map()
  if (!all(data$snp %in% map$snp))
    stop( "Genes were wrongly specified" )
  data_joint <- merge(map, data, by = 'snp', all = TRUE)
  data_order <- data_joint[order(data_joint$chr, data_joint$pos),]
  loci <- with(data_order, ifelse(is.na(add) | is.na(dom), 0, 1))
  add <- with(data_order, ifelse(is.na(add), 0, add))
  dom <- with(data_order, ifelse(is.na(dom), 0, dom))
  vec_loci <- unname(sapply(split(loci, data_order$chr),
                            function(x) paste(x, collapse = '')))
  vec_add <- unname(lapply(split(add, data_order$chr), rev)) # needs reverse order 
  vec_dom <- unname(lapply(split(dom, data_order$chr), rev)) # idem above
  return(.Cpp_trait_quant_ctor(specie, vec_loci, m, vec_add, vec_dom))
}

## \EOF
###########################################################################
