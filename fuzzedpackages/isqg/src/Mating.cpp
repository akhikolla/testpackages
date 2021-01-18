// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*

 This file is part of isqg, a R package for in silico quantitative genetics

              Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
 * Filename: Mating.cpp
 
 * Description: C++ implementations to be used by isqg R package
 
 * Author: Fernando H. Toledo
 
 * Maintainer: Fernando H. Toledo
 
 * Created: We Apr 04 2018
 
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
//////////////////////////////////////////////////////////////////////////

// --Rcpp attributes--
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]

# include <isqg.h>

// --Initializations--

// [[Rcpp::export(name = .Cpp_founder_ctor)]]
Specimen founder(isqg::seamless::Trap<Specie> origin, Code code) {

  Karyotype information(origin->slot->ensemble.size()) ;

  for ( auto it = 0; it < origin->slot->ensemble.size(); it++ )
    information.at(it) = DNA(origin->slot->ensemble.at(it), code) ;

  Specimen individual(origin->slot, information) ;

  return individual ;

}

// [[Rcpp::export(name = .Cpp_import_ctor)]]
Specimen import(isqg::seamless::Trap<Specie> origin, Code cis, Code trans) {

  Karyotype information(origin->slot->ensemble.size()) ;
  
  Codes seq_cis  (origin->split(cis)  ) ;
  Codes seq_trans(origin->split(trans)) ;
  
  for ( auto it = 0; it < origin->slot->ensemble.size(); it++ )
    information.at(it) = DNA(seq_cis.at(it), seq_trans.at(it)) ;
  
  Specimen individual(origin->slot, information) ;
  
  return individual ;

}

// --Mating Designs--

// [[Rcpp::export(name = .Cpp_cross_ctor)]]
Population cross(int number, isqg::seamless::Trap<Specimen> female, isqg::seamless::Trap<Specimen> male) {
  
  // equality comparison
  if ( female->root != male->root ) 
    Rcpp::stop( "Provided Specimens belong to different Species" ) ;

  Population progeny(number, * female) ;

  for ( auto it = 0; it < number; it++ ) {
    female->meiosis() ; 
    male->meiosis() ;
    progeny.at(it) = Specimen(female->root, hybridization(female->nucleous, male->nucleous)) ;
  }
  
  return progeny ;

}

// [[Rcpp::export(name = .Cpp_selfcross_ctor)]]
Population self(int number, isqg::seamless::Trap<Specimen> individual) {

  Specimen clone(* individual) ;

  // quality control
  if ( individual->root != clone.root ) 
    Rcpp::stop( "Reference Specie was lost" ) ;

  Population progeny(number, * individual) ;
  
  for ( auto it = 0; it < number; it++ ) {
    individual->meiosis() ;
    clone.meiosis() ;
    progeny.at(it) = Specimen(clone.root, hybridization(individual->nucleous, clone.nucleous)) ;
  }

  return progeny ;

}

// [[Rcpp::export(name = .Cpp_doublehaploid_ctor)]]
Population dh(int number, isqg::seamless::Trap<Specimen> individual) {

  Population progeny(number, * individual) ;
  
  for ( auto it = 0; it < number; it++ ) {
    individual->meiosis() ;
    progeny.at(it) = Specimen(individual->root, duplication(individual->nucleous)) ;
  }

  return progeny ;

}

// --Breeding Designs--

// \EOF
///////////////////////////////////////////////////////////////////////////
