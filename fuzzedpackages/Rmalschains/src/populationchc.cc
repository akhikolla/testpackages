/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "populationchc.h"
#include <cmath>

using namespace realea;

PopulationRealCHC::PopulationRealCHC(Random *random, unsigned int max, unsigned int pob) : 
    PopulationReal(random,max, pob) {
}

void PopulationRealCHC::restart(DomainRealPtr domain) {
    reset(domain, getBest());
}

tIndividualReal *PopulationRealCHC::getInstance(tChromosomeReal &crom) {
    return new tIndividualRealCHC(crom);
}

tIndividualReal *PopulationRealCHC::getInstance(tChromosomeReal &crom, double fitness) {
    return new tIndividualRealCHC(crom, fitness);
}



tIndividualRealCHC::tIndividualRealCHC(const tChromosomeReal &com) :
    tIndividualReal(com), 
    m_codbin(NULL), m_is_codbin(false) {
}

tIndividualRealCHC::tIndividualRealCHC(const tChromosomeReal &com, double fitness) :
    tIndividualReal(com, fitness), 
    m_codbin(NULL), m_is_codbin(false) {
}


void StringRep(tGen *vect, char *Cad_sal, int length, tGen *max, tGen *min);

void tIndividualRealCHC::calculateBin(DomainRealPtr domain) {
    if (m_is_codbin) {
       return;
    }

    tChromosomeReal sol = this->sol();
    m_codbin_size = sol.size()*BITS_GEN;
    m_codbin = new char[m_codbin_size+1];

    tReal min, max;

    domain->getValues(0, &min, &max);
    // Uso la función anterior
    StringRep(&sol[0], m_codbin, sol.size(), &min, &max);

    for (unsigned i=0; i < sol.size(); ++i) {
        string str(m_codbin + i*BITS_GEN, BITS_GEN);
	bitset<BITS_GEN> genbin(str);
	m_codbin_opt.push_back(genbin);
    }
    m_is_codbin = true;
}

tIndividualRealCHC::~tIndividualRealCHC(void) {
    if (m_codbin) 
       delete[] m_codbin;
}


/**********************************************************/
/* Itoc and Ctoi translate ints to strings and vice versa */
/**********************************************************/
unsigned long int Ctoi(char *Cad_ent, int length)
{
  int i;		
	unsigned long n;	
	
	n = (unsigned long) 0;
	for (i=0; i<length; i++)
	  {
	    n <<= 1;
	    n += (*Cad_ent++ - (int) '0');
	  }
	return(n);
}

void Itoc(unsigned long int n, char *Cad_sal, int length)
{
  int i;		
  
  for (i=length-1; i>=0; i--)
    {
      Cad_sal[i] = '0' + (n & 1);
      n >>= 1;
    }
}


/*****************************************************************/
/* Translations between fixed point ints and reflected Gray code */
/*****************************************************************/

void Gray(char *Cad_ent, char *Cad_sal, int length)
{
  int i;
  char last;

  last = '0';
  for (i=0; i<length; i++)
    {
      Cad_sal[i] = '0' + (Cad_ent[i] != last);
      last = Cad_ent[i];
    }
}

/*************************************************************************/
/* Translations between string representation and floating point vectors */
/*************************************************************************/
void StringRep(tGen *vect, char *Cad_sal, int length, tGen *max, tGen *min)
{
  int i;		
  unsigned long int n;	
  int pos;	
  char tmpstring[BITS_GEN];	
  tGen INCREMENTO;	
  
  pos = 0;
  for (i=0; i < length; i++)
    {
      INCREMENTO=(max[i]-min[i])/(pow(2.0, (double) BITS_GEN) - 1.0);
      
      n = (int) ((vect[i] - min[i]) / INCREMENTO + 0.5);
      
      Itoc(n, tmpstring, BITS_GEN);
      Gray(tmpstring, &Cad_sal[pos], BITS_GEN);
      
      pos += BITS_GEN;
    }
  Cad_sal[pos] = '\0';
}

char *tIndividualRealCHC::getBin(void) {
    if (m_is_codbin)
      return m_codbin;
    else
	throw string("codbin has not been calculated");
}

unsigned tIndividualRealCHC::getBinSize(void) {
    return m_codbin_size;
}

unsigned tIndividualRealCHC::distHamming(tIndividualRealCHC *other) {
    tIndividualRealCHC *mom = this;
    tIndividualRealCHC *dad = other;

   char *momBin = mom->getBin();
    char *dadBin = dad->getBin();
    unsigned dist = 0;

    for (unsigned i = 0; momBin[i] != '\0'; i++) {
	if (momBin[i] != dadBin[i]) {
	    dist += 1;
	}
    }

    return dist;
}

unsigned tIndividualRealCHC::distHammingOpt(tIndividualRealCHC *other) {
    tIndividualRealCHC *mom = this;
    tIndividualRealCHC *dad = other;
    
    vector< bitset<BITS_GEN> > momBin = mom->m_codbin_opt;
    vector< bitset<BITS_GEN> > dadBin = dad->m_codbin_opt;
    unsigned dist = 0;

    for (unsigned i = 0; i < momBin.size(); i++) {
	bitset<BITS_GEN> together = momBin[i] ^dadBin[i];
	dist += together.count();
    }

    return dist;
}
