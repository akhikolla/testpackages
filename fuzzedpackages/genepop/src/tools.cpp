/***************************************************************************
@ F. Rousset 2005-2006
@ J. Lopez 2016-2017

francois.rousset@umontpellier.fr
jimmy.lopez@umontpellier.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/

#include "tools.h"

#include <unistd.h>

#include <sstream>
#include <cstdio> // remove
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h> 
static bool Grestrict = false;
#else
#include <iostream> // cerr
#endif



using namespace std;

int genepop_exit(int error, const char *message){ // error=-1 for error, O for normal
#ifdef COMPATIBILITYRCPP
  throw(Rcpp::exception(message,"tools.cpp",4));
#else
  cerr<<message;
  return(error); // http://stackoverflow.com/questions/30250934/how-to-end-c-code
#endif
}

// [[Rcpp::export]]
void Rinterrupt_genepop() {
  genepop_exit(-1,"Interruption of genepop"); 
}

void clean_temp_file(int locus, int population) {
  remove("poploc");
  for(int jfi=0;jfi<locus;jfi++){
    stringstream locstst;
    locstst<<"locc"<<jfi+1;
    remove(locstst.str().c_str());
  }
  for(int ifi=0;ifi<population;ifi++){
    stringstream locstst;
    locstst<<"popc"<<ifi+1;
    remove(locstst.str().c_str());
  }
  //     remove("P*_L*."); ne marche pas...
  for(int ifi=0;ifi<population;ifi++){
    for(int jfi=0;jfi<locus;jfi++){
      stringstream locstst;
      locstst<<"P"<<ifi+1<<"_L"<<jfi+1;
      remove(locstst.str().c_str());
    }
  }
}

#ifdef COMPATIBILITYRCPP

// [[Rcpp::export]]
void Rset_restriction(bool set) {
  Grestrict = set;
}

void check_restriction(int locus, int population) {
  static int max_locus = 300;
  static int max_population = 300;
  if(Grestrict){
    bool breed = false;
    ostringstream stream;
    stream.clear();
    stream << "The maximum number of loci is " << max_locus << " but the input file contains " << locus << " loci.";
    string str = stream.str();
    if(locus > max_locus) {
      throw(Rcpp::exception(str.c_str()));
      breed = true;
    }
    ostringstream stream2;
    stream2.clear();
    stream2 << "The maximum number of populations is " << max_population << " but the input file contains " << population << " populations.";
    string str2 = stream.str();
    if(population > max_population) {
      throw(Rcpp::exception(str2.c_str()));
      breed = true;
    }
    
    if(breed) {
      clean_temp_file(locus, population);
    }
  }
  
}
#endif
