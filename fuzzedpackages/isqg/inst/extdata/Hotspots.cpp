// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 Sample extention for isqg: a R package for in silico quantitative genetics

              Copyright (C) 2019 Fernando H. Toledo CIMMYT
              
 * Filename: Hotspots.cpp
 
 * Description: extension of meiosis recombination for events in hotspots
 
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

// auxiliar function to sample from a "shifted beta"" distribution
Rcpp::NumericVector shiftbeta(int n, double x) {

  // shift a beta(10, 2.5) draw by x
  return x * Rcpp::rbeta(n, 10.0, 2.5) ;

}

// function to define the custom recombination process
Map hotspot(Chromosome * group) {

  // retrive the length of the chromosome
  double length(group->get_length()) ;

  // sample how many recombination events will occur
  int event(static_cast<int>(R::rpois(length))) ;

  // conditional whenever recombination happens
  if ( event == 0 ) { // case: no recombinations

    // return an empty vector
    return Map() ;

  } else {            // case: with recombinations

    // draw from shifted beta distribution
    Map chiasmata(Rcpp::as<Map>(shiftbeta(event, length))) ;

    // sort events of recombination
    std::sort(chiasmata.begin(), chiasmata.end()) ;

    // return the recombination positions
    return chiasmata ;

  }

}

// wrap the function as a "smart" external pointer
// [[Rcpp::export]] // RCPP attribute
MPtr hotspotp() { return MPtr(new FPtrM(& hotspot), true) ; }

// \EOF
///////////////////////////////////////////////////////////////////////////
