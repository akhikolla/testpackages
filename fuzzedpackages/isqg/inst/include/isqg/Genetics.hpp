///////////////////////////////////////////////////////////////////////////
// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 This file is part of isqg, a R package for in silico quantitative genetics

              Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
 * Filename: Genetics.hpp
 
 * Description: C++ headers to be used by isqg R package
 
 * Author: Fernando H. Toledo
 
 * Maintainer: Fernando H. Toledo
 
 * Created: Mo Jan 22 2018
 
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

# ifndef _QUANTGENETICS_HPP_
# define _QUANTGENETICS_HPP_

# include <isqg/FwdDefs.hpp>
# include <isqg/FwdFuncs.hpp>

class Chromosome {

  // polymorphic meiosis:
  friend class Extended ;
  friend class DNA ;  
  friend class Genome ;
  friend class Specimen ;

public:
  
  // constructors
  Chromosome(void) { } // needs default !
  Chromosome(Map, MPtr) ;

  // retrievers:
  double get_length(void) const { return length ; }
  double get_centromere(void) const { return centromere ; }
  Map    get_map(void) const { return map ; }  
  
  // meiosis intra chromosome
  void  meiosis(void) ;
  Index gamete(void) { return prototype ; }
  Map   pseudo_gamete(void) ;
  Index lazy_gamete(Map, bool) ;
  
private:

  double    length, centromere ; // length in cM & "mid point" position
  Map       map ;                // positions
  Edge      _5p, _3p ;           // iterate over range (extremes)
  Index     prototype ;          // gamete -- independent from the genotype
  Meiosis * trigger ;            // dispatcher for meiosis process
  
} ; // Chromosome

class Meiosis {

public:

  virtual Map meiosis(Chromosome *) = 0 ; // pure abstract :)

} ; // Meiosis

class Extended : public Meiosis {

public:

  Extended(MPtr extension) : process(*(extension)) { }

  Map meiosis(Chromosome *) ; // custom specialization

private:

  FPtrM process ; // holds the user defined process for meiosis

} ; // Extended

class Catalog {

  // comparisons
  friend bool operator==(Catalog, Catalog) ;

public:

  Catalog(Names, Spots, Map, Spots, Spots, Spots) ;

  // features
  Position search(Code) ;
  Codes    split(Code) ;
  
  // iterators
  Names::iterator snps_it(bool start) { return start ? names.begin() : names.end() ; }
  Spots::iterator chrs_it(bool start) { return start ? group.begin() : group.end() ; }
  Map::iterator   loci_it(bool start) { return start ? point.begin() : point.end() ; }
  Spots::iterator posi_it(bool start) { return start ? posis.begin() : posis.end() ; }
  
  // accessors/retrievers
  Names snps(void) const { return names ; }
  Spots chrs(void) const { return group ; }
  Map   loci(void) const { return point ; }
  Spots posi(void) const { return posis ; }

private:

  // genome identification
  Names names  ;
  Spots group  ;
  Map   point  ;
  Spots posis ;
  
  // chromosome limits
  Spots lower ;
  Spots upper ;

} ; // Catalog

class Genome {

  friend class Specimen ;
  friend class Specie ;
  friend class Proxy ;
  
  // interfaced constructors
  friend Specimen founder (isqg::seamless::Trap<Specie>, Code) ;
  friend Specimen import  (isqg::seamless::Trap<Specie>, Code, Code) ;
  
  // equality comparison
  friend bool operator!=(GPtr, GPtr) ;
  
  // accessors/retrievers
  friend Names specie_get_snps   (isqg::seamless::Trap<Specie>) ;
  friend Spots specie_get_chrs   (isqg::seamless::Trap<Specie>) ;
  friend Map   specie_get_loci   (isqg::seamless::Trap<Specie>) ;  
  friend Names specimen_get_snps (isqg::seamless::Trap<Specimen>) ;

public:

  Genome(Maps, Names, Spots, Map, Spots, Spots, Spots, MPtr) ; // extended
  Genome(const Genome &) ;                                     // copy constructor

 // meiosis inter chromosome -- iterate over chromosomes
 void     meiosis(void) ;
 Gamete   gamete(void)  ;
 Maps     pseudo_gamete(void) ;
 Gamete   lazy_gamete(Maps, Guides) ; // sifter
 
 Position search(Code snp) { return directory.search(snp) ; }

private:

  Chip parser_cus(Maps, MPtr) ;
    
  Chip    ensemble ; // would be public
  Catalog directory ;

} ; // Genome

class Specie {

public:

  Specie(Maps, Names, Spots, Map, Spots, Spots, Spots, MPtr) ;
  
  Codes gamete(int) ;
 
  Codes split(Code seq)  { return slot->directory.split(seq) ; }

  GPtr  slot ;
    
private:

  void meiosis(void) { slot->meiosis() ; }

} ; // Specie

class DNA {

  friend class Specimen ;

public:

  DNA(void) { }           // needs default !
  DNA(Chromosome, Code) ; // starting from code: AA, Aa, aA and/or aa
  DNA(DNA, DNA) ;         // starting from parents DNA and self
  DNA(Strand) ;           // starting from single parent dh
  DNA(Code, Code) ;       // starting from string when importing

  // genetic process
  void     meiosis(const Index &) ;
  Strand   recombination(void) ;
  Genotype genotype_num() ;
  Codes    genotype_cod() ;

  // auxiliars 
  void     flip(void)      { cis.flip(); trans.flip() ; } // when mirror
  Strand   dom(void)       { return  cis &  trans ; }     // AA    (dominant loci)
  Strand   het(void)       { return  cis ^  trans ; }     // Aa/aA (heterozigous loci)
  Strand   rec(void)       { return ~cis & ~trans ; }     // aa    (recessive loci)
  Strand   get_cis(void)   { return cis ; }               // get cis
  Strand   get_trans(void) { return cis ; }               // get trans
    
private:

  // double helix
  Strand cis, trans ;
    
  // recombination reference
  mutable Index  arrow ;

} ; // DNA

class Specimen {

  friend class DNA ;
  friend class Trait ;
  friend class Infinitesimal ;
  friend class Quantitative ;
  friend class Custom ;

  // fertilization  
  friend Karyotype  hybridization(Karyotype, Karyotype) ; // Gemome pair2Karyotype
  friend Karyotype  duplication(Karyotype) ;              // Genome2dh

  // interfaced constructors
  friend Specimen   founder(isqg::seamless::Trap<Specie>, Code) ;
  friend Specimen   import(isqg::seamless::Trap<Specie>, Code, Code) ;

  // mating designs
  friend Population cross(int, isqg::seamless::Trap<Specimen>, isqg::seamless::Trap<Specimen>) ;
  friend Population self(int, isqg::seamless::Trap<Specimen>) ;
  friend Population dh(int, isqg::seamless::Trap<Specimen>) ;
  
  // trait
  friend double trait_alpha_eval(isqg::seamless::Trap<Trait>, isqg::seamless::Trap<Specimen>) ;

  // equality comparison
  friend bool operator!=(GPtr, GPtr) ;

  // accessors/retrievers
  friend Names specimen_get_snps(isqg::seamless::Trap<Specimen>) ;
  friend int   specimen_look_num(isqg::seamless::Trap<Specimen>, Code) ;
  friend Code  specimen_look_cod(isqg::seamless::Trap<Specimen>, Code) ;

public:

  // factory
  Specimen(const Specimen &) ;     // copy constructor
  Specimen(GPtr, Karyotype) ;      // default constructor

  // meiosis within individual
  void meiosis(void) ;
 
  // export/codify interface
  Genotype genotype_num(void) ;
  Codes    genotype_cod(void) ;
  
  // auxiliars
  Specimen mirror(void) ;
  Position search(Code snp) { return root->search(snp) ; }
  int      look_num(Code) ;
  Code     look_cod(Code) ;
  Tape     get_cis(void) ;
  Tape     get_trans(void) ;

 private:

  Karyotype parser(Chip, Code) ;   // Genome2Karyotype

  GPtr      root ;                 // pointer to genome
  Karyotype nucleous ;             // stores DNA

} ; // Specimen

class Trait {

  // polymorphic meiosis:
  friend class Infinitesimal ;
  friend class Quantitative ;
  friend class Custom ;
  
  // trait
  friend double trait_alpha_eval(isqg::seamless::Trap<Trait>, isqg::seamless::Trap<Specimen>) ;
  
  // equality comparison
  friend bool operator!=(GPtr, GPtr) ;

public:
  
  // constructors
  Trait(GPtr, Codes, double, double, double) ; // when infinitesimal
  Trait(GPtr, Codes, double, Maps, Maps) ;       // when quantitative
  Trait(GPtr, APtr, SEXP) ;               // when custom
  
  // dispatch trigger and retrieve breeding value
  double  alpha(isqg::seamless::Trap<Specimen>) ;
  
private:

  GPtr    root ;
  Alpha * trigger ; // dispatcher for breeding value
  
} ; // Trait

class Alpha {

public:

  virtual double value(isqg::seamless::Trap<Specimen>) = 0 ; // pure abstract :)

} ; // Alpha

class Infinitesimal : public Alpha {

public:

  Infinitesimal(Codes, double, double, double) ;

  // return Fisher's infinitesimal Specimen's breeding value given "m", "a" and "d"
  double   value(isqg::seamless::Trap<Specimen>) ; 

private:

  Switcher parser(Codes) ;
  
  double   mean, additive, dominant ;
  
  Switcher genes ;

} ; // Infinitesimal

class Quantitative : public Alpha {

public:

  Quantitative(Codes, double, Maps, Maps) ;

  // return quantitative Specimen's breeding value given "m", "a"s and "d"s
  double value(isqg::seamless::Trap<Specimen>) ; 

private:

  Switcher parser(Codes) ;

  double   mean ;  
  
  Maps     additives, dominants ;
  
  Switcher genes ;

} ; // Quantitative

class Custom : public Alpha {

public:

  Custom(APtr extension, SEXP auxiliar) : eval(*(extension)), extra(auxiliar) { }

  double value(isqg::seamless::Trap<Specimen>) ; // return user defined Specimen's breeding value

private:

  FPtrA eval ;  // holds the user defined breeding value
  
  SEXP  extra ; // holds any extra information to generate breeding values

} ; // Custom

# endif // _GENETICS_HPP_

// \EOF
///////////////////////////////////////////////////////////////////////////
