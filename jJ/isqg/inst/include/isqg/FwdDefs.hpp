// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 This file is part of isqg, a R package for in silico quantitative genetics

              Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
 * Filename: FwdDefs.hpp
 
 * Description: C++ headers to be used by isqg R package
 
 * Author: Fernando H. Toledo
 
 * Maintainer: Fernando H. Toledo
 
 * Created: Fr Jan 12 2018
 
 * Updated: -
 
  This program is free software; you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation; either version 2 of the License, or 
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but 
  WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software Foundation, 
  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
                                                        
  `` Far better an approximate answer to the right question, which is 
  often vague, than the exact answer to the wrong question, which can
  always be made precise ''
                         --John Tukey, Ann. Math. Stat. 33(1):13 1962
*/
///////////////////////////////////////////////////////////////////////////

# ifndef _FORWARD_ALIAS_HPP_
# define _FORWARD_ALIAS_HPP_

// types' alias and forward declarations:
 
typedef std::vector<double>                Map ;

typedef std::vector<Map>                   Maps ;

typedef std::vector<bool>                  Guides ;

typedef Map::iterator                      Edge ;

// equal to Strand -- kept up just to track/debug
typedef boost::dynamic_bitset<>            Index ;

typedef std::vector<Index>                 Gamete ;

class                                      Chromosome ;

// runtime polymorphism device to enable user defined meiosis
class                                      Meiosis ; // abstraction for user defined meiosis

// function pointer for user defined function
typedef Map (* FPtrM) (Chromosome *) ;  

typedef Rcpp::XPtr<FPtrM>                  MPtr ; // external/smart pointer for meiosis

class                                      Extended ; // w/ user defined C++ function

typedef std::vector<Chromosome>            Chip ;     

// Alias within Catalog
typedef std::vector<std::string>           Names ;

typedef std::vector<int>                   Spots ;

typedef std::string                        Code ;

typedef std::tuple<Code, int, double, int> Position ;

class                                      Catalog ;

class                                      Genome ;

typedef Rcpp::XPtr<Genome>                 GPtr ; // external/smart pointer for Genome

class                                      Specie ; // slot for GPtr

// equal to Index -- kept up just to track/debug
typedef boost::dynamic_bitset<>            Strand ;

// equal to Gamete -- kept up just to track/debug
typedef std::vector<Strand>                Tape ;

class                                      DNA ;

typedef std::vector<DNA>                   Karyotype ;

class                                      Specimen ;

typedef std::vector<Specimen>              Population ;

typedef std::vector<Code>                  Codes ;

typedef std::vector<int>                   Genotype ;

// equal to Strand -- kept up just to track/debug
typedef std::vector<Index>                 Switcher ;

class                                      Trait ;

// runtime polymorphism device to enable user defined breeding value
class                                      Alpha ; // abstraction for user defined alpha

class                                      Infinitesimal ;

class                                      Quantitative ;

typedef double (* FPtrA) (Specimen) ;      // function pointer for user defined function

typedef Rcpp::XPtr<FPtrA>                  APtr ; // external/smart pointer for alpha

class                                      Custom ; // [[TODO]]

# endif // _FORWARD_ALIAS_HPP_

// \EOF
///////////////////////////////////////////////////////////////////////////
