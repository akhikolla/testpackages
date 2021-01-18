// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 This file is part of isqg, a R package for in silico quantitative genetics

              Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
 * Filename: FwdFuncs.hpp
 
 * Description: C++ headers to be used by isqg R package
 
 * Author: Fernando H. Toledo
 
 * Maintainer: Fernando H. Toledo
 
 * Created: Fr Mar 09 2018
 
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

# ifndef _FORWARD_FUNCS_HPP_
# define _FORWARD_FUNCS_HPP_

// --Internals--

bool       operator!=(Catalog, Catalog) ;

bool       operator!=(GPtr, GPtr) ;

Karyotype  hybridization(Karyotype, Karyotype) ;

Karyotype  duplication(Karyotype) ;

// --Initializations--

Specimen   founder(isqg::seamless::Trap<Specie>, Code) ;

Specimen   import(isqg::seamless::Trap<Specie>, Code, Code) ;

// --Mating Designs--

Population cross(int, isqg::seamless::Trap<Specimen>, isqg::seamless::Trap<Specimen>) ;

Population self(int, isqg::seamless::Trap<Specimen>) ;

Population dh(int, isqg::seamless::Trap<Specimen>) ;

// --Breeding Desings--

// --Externals--

MPtr       standard_meiosis(void) ;

Codes      gamete_ctor(int, isqg::seamless::Trap<Specie>) ;

// Specie     specie_std_ctor(Maps, Names, Spots, Map, Spots, Spots, Spots) ;

Specie     specie_cus_ctor(Maps, Names, Spots, Map, Spots, Spots, Spots, MPtr) ;

Names      specie_get_snps(isqg::seamless::Trap<Specie>) ;

Spots      specie_get_chrs(isqg::seamless::Trap<Specie>) ;

Map        specie_get_loci(isqg::seamless::Trap<Specie>) ;

Names      specimen_get_snps(isqg::seamless::Trap<Specimen>) ;

int        specimen_look_num(isqg::seamless::Trap<Specimen>, Code) ;

Code       specimen_look_cod(isqg::seamless::Trap<Specimen>, Code) ;

Specimen   specimen_mirror(isqg::seamless::Trap<Specimen>) ;

Genotype   genotype_num_ctor(isqg::seamless::Trap<Specimen>) ;

Codes      genotype_cod_ctor(isqg::seamless::Trap<Specimen>) ;

Trait      trait_infty_ctor(isqg::seamless::Trap<Specie>, Codes, double, double, double) ;

Trait      trait_quant_ctor(isqg::seamless::Trap<Specie>, Codes, double, Maps, Maps) ;

Trait      trait_custm_ctor(isqg::seamless::Trap<Specie>, APtr, SEXP) ;

double     trait_alpha_eval(isqg::seamless::Trap<Trait>, isqg::seamless::Trap<Specimen>) ;

# endif // _FORWARD_FUNCS_HPP_

// \EOF
///////////////////////////////////////////////////////////////////////////
