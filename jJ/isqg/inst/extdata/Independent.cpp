// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 Sample extention for isqg: a R package for in silico quantitative genetics

              Copyright (C) 2019 Fernando H. Toledo CIMMYT
              
 * Filename: Independent.cpp
 
 * Description: extension of meiosis recombination for independent events
 
 * Author: Fernando H. Toledo
 
 * Created: Tue Dec 3 2019
 
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

// [[Rcpp::depends(isqg)]] // RCPP attribute

// loading headers
# include <isqg.h>    // from the "isqg" package
# include <vector>    // support for vectors from standard template library
# include <algorithm> // support for algorithms from standard template library

// function to define the custom recombination process
Map indep(Chromosome * group) {

  // retrieve the map of the chromosome  
  Map map(group->get_map()) ;

  // loop over all positions of the map and... 
  // raffle (Bernoulli draw) if each position will recombine or not
  for (auto it = 0; it < map.size(); it++)
    if (static_cast<bool>(R::rbinom(1., .5))) map.at(it) = 2. + 1. ;

  // remove any contiguous positions where recombination happened
  map.erase(std::remove(map.begin(), map.end(), 2. + 1.), map.end()) ;

  // return the recombination positions
  return map ; 
  
}

// wrap the function as a "smart" external pointer
// [[Rcpp::export]] // RCPP attribute
MPtr indepp() { return MPtr(new FPtrM(& indep), true) ; }

// \EOF
///////////////////////////////////////////////////////////////////////////
