## -*- mode: R -*-
###########################################################################
## This file is part of isqg: R package for in silico quantitative genetics
##
##              Copyright (C) 2018 Fernando H. Toledo CIMMYT
##              
## * Filename: R6Classes.R
## 
## * Description: Defines seamless R6 classes (pointer holders)
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

##' @title Class providing object with methods to mimic in silico Genomes
##'
##' @name Specie
##'
##' @description Mean to mimic a in silico Genomes. It is the machine instances of
##'     the simulator.
##' 
##' @return Objects of R6 class with methods to mimic in silico Genomes.
##'
##' @details Object of R6 class that points to C++ objetcs. 
##'
##' @rdname Specie
##'
##' @docType class
NULL 

##' @rdname Specie
.R_Specie_ctor <-
  R6::R6Class(
    classname = "Specie",
    inherit = NULL,
    portable = TRUE,
    public = list(
      ##' @field .ptr External pointer to the instance of the C++ class Specie. 
      .ptr = NULL,
      ##' @description Create an instance of a Specie.
      ##' @param ptr an Smart pointer to an instance of a Specie C++ class.
      ##' @return A new `Specie` object.
      initialize = function(ptr) {
        self$.ptr <- ptr
      },
      ##' @description Print/Show an instance of the Specie class.
      ##' @param ... further arguments to be passed to print.
      print = function(...) {
        cat("<isqg Specie id> ", capture.output(self$.ptr), "\n", sep = "")
        invisible(self)
      },
      ##' @description Constructor of a Founder Instances of the Specimen Class
      ##' @param code   a length one character vector with one of the genotype codes:
      ##'     "AA", "Aa", "aA" or "aa".
      founder = function(code) {
        return(.Cpp_founder_ctor(self, code))
      }, 
      ##' @description Generate gamete prototypes for this specie.
      ##' @param n an integer number with the number of gametes prototype to be generated.
      ##' @return a vector of strings with the gamete prototypes.
      gamete = function(n) {
        return(.Cpp_Gamete_ctor(n, self))
      },
      ##' @description Retrieve the map of the current specie.
      ##' @return a data.frame with the map of the specie.
      map = function() {
        df <- data.frame(snp = .Cpp_spc_snps(self), 
                         chr = .Cpp_spc_chrs(self) + 1, ## return indices to R 
                         pos = .Cpp_spc_loci(self),
                         stringsAsFactors = FALSE)      ## snps as characters
        return(df)
      }
    )
  )

##' @title Class providing object with methods to mimic in silico Specimens
##'
##' @name Specimen
##'
##' @description Mean to mimic a in silico Specimens. It is the working instances of
##'     the simulator.
##' 
##' @return Objects of R6 class with methods to mimic in silico Specimens.
##'
##' @details Object of R6 class that points to C++ objetcs. 
##'
##' @rdname Specimen
##'
##' @docType class
NULL

##' @rdname Specimen
.R_Specimen_ctor <- 
  R6::R6Class(
    classname = 'Specimen',
    inherit = NULL,
    portable = TRUE,
    public = list(
      ##' @field .ptr External pointer to the instance of the C++ class Specie. 
      .ptr = NULL,
      ##' @description Create an instance of a Specimen.
      ##' @param ptr an Smart pointer to an instance of a Specimen C++ class.
      ##' @return A new `Specimen` object.
      initialize = function(ptr) {
        self$.ptr <- ptr
      },
      ##' @description Print/Show an instance of the Specimen class.
      ##' @param ... further arguments to be passed to print.
      print = function(...) {
        cat("<isqg Specimen id> ", capture.output(self$.ptr), "\n", sep = "")
        invisible(self)
      },
      ##' @description Evaluates the breeding value.
      ##' @param trait an instance of the class Trait.
      ##' @return the breeding value of the specimen for the given trait.
      alpha = function(trait) {
        return(.Cpp_trait_alpha_eval(trait, self))
      },
      ##' @description Codify Specimen's Genotypes.
      ##' @param phase logical should the codes keep the phase. 
      ##' @return A numeric or character vector with the codified Specimen's genotypes.
      genotype = function(phase = FALSE) {
        if (phase) {
          codes <- .Cpp_Genotype_cod(self)
          names(codes) <- .Cpp_get_snps(self)
          return(codes)
        } else {
          codes <- .Cpp_Genotype_num(self)
          names(codes) <- .Cpp_get_snps(self)
          return(codes)
        }
      },
      ##' @description Performs the simple bi-parental cross.
      ##' @param n a length-one integer vector with the size of the progeny.
      ##' @param gid instance of the class specimen which will be used to mate.
      ##' @return a size \emph{n} list with instances of the class Specimen that 
      ##'     represent new individuals belonging to the progeny of the respective mating 
      ##'     scheme.
      cross = function(n = 1, gid) {
        return(cross(n, self, gid))
      },
      ##' @description Performs the selfcross.
      ##' @param n a length-one integer vector with the size of the progeny.
      ##' @param replace logical scalar indicating if the outcome of the function will
      ##'      replace the current instance of the Specimen
      ##' @return a size \emph{n} list with instances of the class Specimen that 
      ##'     represent new individuals belonging to the progeny of the respective mating 
      ##'     scheme.
      selfcross = function(n = 1, replace = FALSE) {
        if (n == 1 && replace) {
          self <<- selfcross(n, self)
          return(invisible(self))
        } else {
          return(selfcross(n, self))
        }
      },
      ##' @description Performs the double-haploid duplication
      ##' @param n a length-one integer vector with the size of the progeny.
      ##' @param replace logical scalar indicating if the outcome of the function will 
      ##'     replace the current instance of the Specimen
      ##' @return a size \emph{n} list with instances of the class Specimen that 
      ##'     represent new individuals belonging to the progeny of the respective mating 
      ##'     scheme.
      dh = function(n = 1, replace = FALSE) {
        if (n == 1 && replace) {
          self <<- dh(n, self)
          return(invisible(self))
        } else {
          return(dh(n, self))
        }
      },
      ##' @description Generates a `mirrored` specimen.
      ##' @return an instance of the Specimen class with all loci mirrored.
      mirror = function() {
        return(.Cpp_specimen_mirror(self))
      },
      ##' @description Acess specific locus' value from specimen.
      ##' @param snp an character string with the name of the locus to lookup.
      ##' @param phase logical should the codes keep the phase. 
      ##' @return the genotype of the given locus.
      look = function(snp, phase = FALSE) {
        if (phase) {
          return(.Cpp_look_cod(self, snp))
        } else {
          return(.Cpp_look_num(self, snp))
        }
      }
    )
  )

##' @title Class providing object with methods to mimic in silico Traits
##'
##' @name Trait
##'
##' @description Mean to mimic a in silico Trait. It is the working instances of
##'     the simulator.
##' 
##' @return Objects of R6 class with methods to mimic in silico Traits.
##'
##' @details Object of R6 class that points to C++ objetcs. 
##'
##' @rdname Trait
##'
##' @docType class
NULL

##' @rdname Trait
.R_Trait_ctor <- 
  R6::R6Class(
    classname = 'Trait',
    inherit = NULL,
    portable = TRUE,
    public = list(
      ##' @field .ptr External pointer to the instance of the C++ class Specie. 
      .ptr = NULL,
      ##' @description Create an instance of a Trait.
      ##' @param ptr an Smart pointer to an instance of a Trait C++ class.
      ##' @return A new `Trait` object.
      initialize = function(ptr) {
        self$.ptr <- ptr
      },
      ##' @description Print/Show an instance of the Trait class.
      ##' @param ... further arguments to be passed to print.
      print = function(...) {
        cat("<isqg Trait id> ", capture.output(self$.ptr), "\n", sep = "")
        invisible(self)
      },
      ##' @description Evaluates the breeding value.
      ##' @param gid an instance of the class Specimen.
      ##' @return the breeding value of the given specimen for the trait.
      alpha = function(gid) { 
        return(.Cpp_trait_alpha_eval(self, gid)) 
      }
    )
  )
                           
## \EOF
###########################################################################
