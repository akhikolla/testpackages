// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
  This file is part of isqg, a R package for in silico quantitative genetics

  Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
  * Filename: Genetics.cpp
 
  * Description: C++ implementations to be used by isqg R package
 
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
//////////////////////////////////////////////////////////////////////////

# include <isqg.h>

// --Chromosome--  
Chromosome::Chromosome(Map input, MPtr custom) :
  
  length     (input.back()),
  centromere (length),
  map        (input),
  _5p        (map.begin()),
  _3p        (map.end()),
  prototype  (map.size()),
  trigger    (new Extended(custom)) // regular pointer

{ }

void Chromosome::meiosis(void) {

  prototype.reset() ;
  
  Map chiasmata(trigger->meiosis(this)) ;

  if ( chiasmata.size() > 0 ) { // only if crossing happens

    int breaks(0) ; // holds the breaks' distance from _5p

    // right shift XOR chain -- find interval
    for ( auto && chiasma : chiasmata ) {
  
      breaks     = std::distance(_5p, std::upper_bound(_5p, _3p, chiasma)) ;
      prototype ^= (Index(map.size()).set() >> breaks) ;

    }
  
  }
  
  if (static_cast<bool>(R::rbinom(1.0, 0.5))) // raffle w.r.t. reference strand
    prototype.flip() ;

}

// ref: Karlin & Liberman (1978) [Proc. Natl. Acad. Sci. 75(12):6332--6336]
Map count_location(Chromosome * chromosome) {

  double length(chromosome->get_length()) ;

  int event(static_cast<int>(R::rpois(length))) ;
  
  if ( event == 0 ) {
  
    return Map() ;
  
  } else {
  
    // sugar :)
    Map chiasmata(Rcpp::as<Map>(Rcpp::runif(event, 0.0, length))) ;     
    
    std::sort(chiasmata.begin(), chiasmata.end()) ;
  
    return chiasmata ;
    
  }
  
}

Map Chromosome::pseudo_gamete(void) { return trigger->meiosis(this) ; }

Index Chromosome::lazy_gamete(Map chiasmata, bool tape) {

  prototype.reset() ;
  
  if ( chiasmata.size() > 0 ) { // only if crossing happens
  
    int breaks(0) ; // holds the breaks' distance from _5p

    // right shift XOR chain -- find interval
    for (auto && chiasma : chiasmata) {
    
      breaks     = std::distance(_5p, std::upper_bound(_5p, _3p, chiasma)) ;
      prototype ^= (Index(map.size()).set() >> breaks) ;
    
    }
  
  }
  
  if (tape) // raffle w.r.t. reference strand
    prototype.flip() ;
    
  return prototype ;

}

// [[Rcpp::export(name = .Cpp_meiosis_standard)]]
MPtr standard_meiosis() { return MPtr(new FPtrM(& count_location), true) ; }

// --Extended Meiosis--
Map Extended::meiosis(Chromosome * chromosome) { 

  return process(chromosome) ;

}

// --Catalog--
Catalog::Catalog(Names snps, Spots chrs, Map loci, Spots index, Spots lwr, Spots upr) :

  names (snps),
  group (chrs),
  point (loci),
  posis (index),
  lower (lwr),
  upper (upr)

{ }
  
Position Catalog::search(Code snp) {

  int index(std::distance(names.begin(), std::find(names.begin(), names.end(), snp))) ;
  
  if (index >= names.size())
    Rcpp::stop( "Provided 'snp' doesn't found" ) ;

  return std::make_tuple(names.at(index), group.at(index), point.at(index), posis.at(index)) ;

}
  
Codes Catalog::split(Code seq) {

  Codes vec_seq(lower.size()) ;

  for ( auto it = 0; it < lower.size(); it++ ) {
    vec_seq.at(it) = seq.substr(lower.at(it), upper.at(it) - lower.at(it) + 1) ;
  }

  return vec_seq ;

}

Genome::Genome(Maps input, Names snps, Spots chrs, Map loci, Spots index, Spots lwr, Spots upr, MPtr custom) : 

  ensemble  (parser_cus(input, custom)),
  directory (snps, chrs, loci, index, lwr, upr) 
  
{ }
  
Genome::Genome(const Genome & original) : 

  ensemble  (original.ensemble),
  directory (original.directory)
  
{ }

Chip Genome::parser_cus(Maps input, MPtr custom) {

  int size(input.size()) ;
  
  Chip information(size) ;
  
  for ( auto it = 0; it < size; it++ )
    information.at(it) = Chromosome(input.at(it), custom) ;
  
  return information ;

}

void Genome::meiosis(void) { for ( auto && chr : ensemble ) { chr.meiosis() ; } }

Gamete Genome::gamete(void) {

  Gamete prototypes(ensemble.size()) ;
  
  for (auto it = 0; it < ensemble.size(); it++) 
    prototypes.at(it) = ensemble.at(it).gamete() ;
    
  return prototypes ;

}

Maps Genome::pseudo_gamete(void) {

  Maps chiasmatas(ensemble.size()) ;
  
  for (auto it = 0; it < ensemble.size(); it++)
    chiasmatas.at(it) = ensemble.at(it).pseudo_gamete() ;

  return chiasmatas ;

}

Gamete Genome::lazy_gamete(Maps points, Guides tapes) {

  Gamete prototypes(ensemble.size()) ;
  
  for (auto it = 0; it < ensemble.size(); it++)
    prototypes.at(it) = ensemble.at(it).lazy_gamete(points.at(it), tapes.at(it)) ;
    
  return prototypes ;

}

// --Specie--
Specie::Specie(Maps input, Names snps, Spots chrs, Map loci, Spots index, Spots lwr, Spots upr, MPtr custom) : 

  slot(new Genome(input, snps, chrs, loci, index, lwr, upr, custom), true) 
  
{ }

Codes Specie::gamete(int number) {

  Codes sample(number) ;
  
  for (auto nt = 0; nt < number; nt++) {

    this->meiosis() ;

    Gamete gametes(this->slot->gamete()) ;
    Index reference(gametes.at(0)) ;
 
    for (auto it = 1; it < gametes.size(); it++) {
    
      reference.resize(reference.size() + gametes.at(it).size()) ;
      reference <<= gametes.at(it).size() ;
      Index vante(gametes.at(it)) ;
      vante.resize(reference.size()) ;
      reference |=  vante ;
      
    }

    Code guide ;
    boost::to_string(reference, guide) ;
    
    sample.at(nt) = guide ;
    
  }

  return sample ;

}

// [[Rcpp::export(name = .Cpp_Gamete_ctor)]]
Codes gamete_ctor(int number, isqg::seamless::Trap<Specie> spc) { return spc->gamete(number) ; }

// [[Rcpp::export(name = .Cpp_Specie_cus_ctor)]]
Specie specie_cus_ctor(Maps input, Names snps, Spots chrs, Map loci, Spots index, Spots lwr, Spots upr, MPtr custom) { 

  return Specie(input, snps, chrs, loci, index, lwr, upr, custom) ; 
  
}

// [[Rcpp::export(name = .Cpp_spc_snps)]]
Names specie_get_snps(isqg::seamless::Trap<Specie> spc) { return spc->slot->directory.snps() ; }

// [[Rcpp::export(name = .Cpp_spc_chrs)]]
Spots specie_get_chrs(isqg::seamless::Trap<Specie> spc) { return spc->slot->directory.chrs() ; }

// [[Rcpp::export(name = .Cpp_spc_loci)]]
Map   specie_get_loci(isqg::seamless::Trap<Specie> spc) { return spc->slot->directory.loci() ; }

// --DNA--
DNA::DNA(Chromosome chromosome, Code genotype) : 

  arrow(Index(chromosome.map.size())) 
  
{

  if (genotype == "AA") {

    cis   = Strand(chromosome.map.size()).set() ;
    trans = Strand(chromosome.map.size()).set() ;

  } else if (genotype == "Aa") {

    cis   = Strand(chromosome.map.size()).set() ;
    trans = Strand(chromosome.map.size()) ;

  } else if (genotype == "aA") {

    cis   = Strand(chromosome.map.size()) ;
    trans = Strand(chromosome.map.size()).set() ;

  } else if (genotype == "aa") {

    cis   = Strand(chromosome.map.size()) ;
    trans = Strand(chromosome.map.size()) ;

  } else {

    Rcpp::stop("Unable to initialize genotype with the provided code") ;

  } 
    
}

DNA::DNA(DNA female, DNA male) :

  cis   (female.recombination()),
  trans (male.recombination()),
  arrow (female.cis.size())

{ }
  
DNA::DNA(Strand haploid) :

  cis   (haploid),
  trans (haploid),
  arrow (haploid.size())
  
{ }
  
DNA::DNA(Code cis_seq, Code trans_seq) :

  cis   (Index(cis_seq)),
  trans (Index(trans_seq)),
  arrow (Index(cis_seq.size()))
  
{ }

void DNA::meiosis(const Index & prototype) { arrow = prototype ; } // [change 04/17/2018]

Strand DNA::recombination(void) { return (arrow & cis) | (~arrow & trans) ; }

Genotype DNA::genotype_num() {

  Genotype genotype(cis.size()) ;
  
  for (auto it = 0; it < cis.size(); it++)
    genotype[it] = cis[it] & trans[it] ? 1 : (cis[it] ^ trans[it] ? 0 : -1) ;

  std::reverse(genotype.begin(), genotype.end()) ;

  return genotype ;

}

Codes DNA::genotype_cod() {

  Codes genotype(cis.size()) ;
  
  for (auto it = 0; it < cis.size(); it++)
    genotype[it] = cis[it] & trans[it] ? "1 1" : (~cis[it] & ~trans[it] ? "2 2" : (cis[it] & 1 ? "1 2" : "2 1")) ;
    
  std::reverse(genotype.begin(), genotype.end()) ;
    
  return genotype ;

}

// --Specimen--
Specimen::Specimen(const Specimen & original) :

  root     (original.root),
  nucleous (original.nucleous)
  
{ }

Specimen::Specimen(GPtr origin, Karyotype information) :

  root     (origin),
  nucleous (information)

{ }
  
void Specimen::meiosis(void) {

  root->meiosis() ;
  
  for ( auto it = 0; it < root -> ensemble.size(); it++ )
    nucleous.at(it).meiosis(root -> ensemble.at(it).prototype) ;

}

Genotype Specimen::genotype_num(void) {

  Genotype genotypes ;
  
  for (auto it = 0; it < nucleous.size(); it++) {
    Genotype re(nucleous.at(it).genotype_num()) ;
    genotypes.insert(genotypes.end(), re.begin(), re.end()) ;
  }

  return genotypes ;

}

Codes Specimen::genotype_cod(void) {

  Codes genotypes ;
  
  for (auto it = 0; it < nucleous.size(); it++) {
    Codes re(nucleous.at(it).genotype_cod()) ;
    genotypes.insert(genotypes.end(), re.begin(), re.end()) ;
  }
  
  return genotypes ;

}

Specimen Specimen::mirror(void) {

  Specimen clone(*this) ;
  
  for (auto it = 0; it < nucleous.size(); it++)
    clone.nucleous.at(it).flip() ;
    
  return clone ;

}

int Specimen::look_num(Code snp) {

  Position id(this->search(snp)) ;
      
  bool cis(nucleous.at(std::get<1>(id)).cis[std::get<3>(id)]) ;
  bool trans(nucleous.at(std::get<1>(id)).trans[std::get<3>(id)]) ;
  
  return cis && trans ? 1 : (cis != trans ? 0 : -1) ;

}

Code Specimen::look_cod(Code snp) {

  Position id(this->search(snp)) ;
      
  bool cis(nucleous.at(std::get<1>(id)).cis[std::get<3>(id)]) ;
  bool trans(nucleous.at(std::get<1>(id)).trans[std::get<3>(id)]) ;
    
  return cis && trans ? "1 1" : (!cis && !trans ? "2 2" : (cis ? "1 2" : "2 1")) ;

}

Tape Specimen::get_cis(void) {

  Tape cis(nucleous.size()) ;
  
  for (auto it = 0; it < nucleous.size(); it++)
    cis.at(it) = nucleous.at(it).get_cis() ;
    
  return cis ;

}

Tape Specimen::get_trans(void) {

  Tape trans(nucleous.size()) ;
  
  for (auto it = 0; it < nucleous.size(); it++)
    trans.at(it) = nucleous.at(it).get_trans() ;
    
  return trans ;

}

// [[Rcpp::export(name = .Cpp_Genotype_num)]]
Genotype genotype_num_ctor(isqg::seamless::Trap<Specimen> gid) { return gid->genotype_num() ; }

// [[Rcpp::export(name = .Cpp_Genotype_cod)]]
Codes genotype_cod_ctor(isqg::seamless::Trap<Specimen> gid) { return gid->genotype_cod() ; }

// [[Rcpp::export(name = .Cpp_get_snps)]]
Names specimen_get_snps(isqg::seamless::Trap<Specimen> gid) { return gid->root->directory.snps() ; }

// [[Rcpp::export(name = .Cpp_look_num)]]
int   specimen_look_num(isqg::seamless::Trap<Specimen> gid, Code snp) { return gid->look_num(snp) ; }

// [[Rcpp::export(name = .Cpp_look_cod)]]
Code  specimen_look_cod(isqg::seamless::Trap<Specimen> gid, Code snp) { return gid->look_cod(snp) ; }

// [[Rcpp::export(name = .Cpp_specimen_mirror)]]
Specimen   specimen_mirror(isqg::seamless::Trap<Specimen> gid) { return gid->mirror() ; }

// --Internals--
bool operator!=(Catalog A, Catalog B) {

  bool snps(std::equal(A.snps_it(1), A.snps_it(0), B.snps_it(1))) ;
  bool chrs(std::equal(A.chrs_it(1), A.chrs_it(0), B.chrs_it(1))) ;
  bool posi(std::equal(A.posi_it(1), A.posi_it(0), B.posi_it(1))) ;
  bool loci(std::equal(A.loci_it(1), A.loci_it(0), B.loci_it(1))) ;

  return !(snps && chrs && loci && posi) ;

}

bool operator!=(GPtr A, GPtr B) { return A->directory != B->directory ; }

Karyotype hybridization(Karyotype female, Karyotype male) {

  Karyotype information(female.size()) ;
  
  for ( auto it = 0; it < female.size(); it++ )
    information.at(it) = DNA(female.at(it), male.at(it)) ;
  
  return information ;

}

Karyotype duplication(Karyotype individual) {

  Karyotype information(individual.size()) ;
  
  for ( auto it = 0; it < individual.size(); it++ )
    information.at(it) = DNA(individual.at(it).recombination()) ;
  
  return information ;

}

// \EOF
///////////////////////////////////////////////////////////////////////////
