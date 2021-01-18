// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 This file is part of isqg, a R package for in silico quantitative genetics

              Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
 * Filename: Trait.cpp
 
 * Description: C++ implementations to be used by isqg R package
 
 * Author: Fernando H. Toledo
 
 * Maintainer: Fernando H. Toledo
 
 * Created: Mo Apr 09 2018
 
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

// --Trait--
Trait::Trait(GPtr origin, Codes loci, double mu, double add, double dom) :

  root    (origin),  
  trigger (new Infinitesimal(loci, mu, add, dom)) // regular pointer

  { }
  
Trait::Trait(GPtr origin, Codes loci, double mu, Maps add, Maps dom) :

  root    (origin),  
  trigger (new Quantitative(loci, mu, add, dom)) // regular pointer

  { }
  
Trait::Trait(GPtr origin, APtr extension, SEXP auxiliar) :
  
  root    (origin),
  trigger (new Custom(extension, auxiliar)) // regular pointer

  { }
  
double Trait::alpha(isqg::seamless::Trap<Specimen> gid) { return trigger->value(gid) ; }

// --Infinitesimal Trait--
Infinitesimal::Infinitesimal(Codes loci, double mu, double add, double dom) : 

  mean(mu),
  additive(add), 
  dominant(dom),
  genes(parser(loci))
  
  { }

double Infinitesimal::value(isqg::seamless::Trap<Specimen> gid) {

  double alpha(mean) ;
  
  for (auto it = 0; it < gid->nucleous.size(); it++) {
    alpha += (genes.at(it) & gid->nucleous.at(it).dom()).count() *   additive ;
    alpha += (genes.at(it) & gid->nucleous.at(it).het()).count() *   dominant ;
    alpha += (genes.at(it) & gid->nucleous.at(it).rec()).count() * - additive ;  
  }
  
  return alpha ;

}

Switcher Infinitesimal::parser(Codes loci) {

  Switcher index(loci.size()) ;
  
  for (auto it = 0; it < loci.size(); it++)
    index.at(it) = Index(loci.at(it)) ;
    
  return index ;

} 

// --Quantitative Trait--
Quantitative::Quantitative(Codes loci, double mu, Maps adds, Maps doms) : 

  mean(mu), 
  additives(adds), 
  dominants(doms),
  genes(parser(loci))  
  
  { }

double Quantitative::value(isqg::seamless::Trap<Specimen> gid) {

  double alpha(mean) ;

  for (auto it = 0; it < gid->nucleous.size(); it++) {

    Strand genotype(genes.at(it) & gid->nucleous.at(it).dom()) ;
    auto index(genotype.find_first()) ;
    while (index != genotype.npos) {
      alpha += additives.at(it).at(index) ;
      index  = genotype.find_next(index) ;
    }
    
    genotype = genes.at(it) & gid->nucleous.at(it).het() ;
    index    = genotype.find_first() ;
    while (index != genotype.npos) {
      alpha += dominants.at(it).at(index) ;
      index  = genotype.find_next(index) ;
    }
    
    genotype = genes.at(it) & gid->nucleous.at(it).rec() ;
    index    = genotype.find_first() ;
    while (index != genotype.npos) {
      alpha += -additives.at(it).at(index) ;
      index  =  genotype.find_next(index) ;
    }
  
  }
  
  return alpha ;

}

Switcher Quantitative::parser(Codes loci) {

  Switcher index(loci.size()) ;
  
  for (auto it = 0; it < loci.size(); it++)
    index.at(it) = Index(loci.at(it)) ;
    
  return index ;

}

// --Extended Meiosis--
double Custom::value(isqg::seamless::Trap<Specimen> gid) { return eval(* gid) ; }

// [[Rcpp::export(name = .Cpp_trait_infty_ctor)]]
Trait trait_infty_ctor(isqg::seamless::Trap<Specie> origin, Codes loci, double mu, double add, double dom) {

  return Trait(origin->slot, loci,  mu, add, dom) ;

}

// [[Rcpp::export(name = .Cpp_trait_quant_ctor)]]
Trait trait_quant_ctor(isqg::seamless::Trap<Specie> origin, Codes loci, double mu, Maps add, Maps dom) {

  return Trait(origin->slot, loci, mu, add, dom) ;

}

// [[Rcpp::export(name = .Cpp_trait_custm_ctor)]]
Trait trait_custm_ctor(isqg::seamless::Trap<Specie> origin, APtr extension, SEXP auxiliar) {

  return Trait(origin->slot, extension, auxiliar) ;

}

// [[Rcpp::export(name = .Cpp_trait_alpha_eval)]]
double trait_alpha_eval(isqg::seamless::Trap<Trait> trait, isqg::seamless::Trap<Specimen> gid) {

  // equality comparison
  if ( trait->root != gid->root ) 
    Rcpp::stop( "Provided Trait and Specimen belong to different Species" ) ;

  return trait->alpha(gid) ;
  
}

// \EOF
///////////////////////////////////////////////////////////////////////////
