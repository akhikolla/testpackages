/***************************************************************************
@ F. Rousset 2005-2007

francois.rousset@umontpellier.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

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
#include <iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<limits>
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
#endif
#include "genepop.h"
#include "GenepopS.h"
#include "F_est.h"
#include "multimig.h"
#include "tools.h"



using namespace std;

namespace multimig {
  // pas facile a retirer car deviendrait argument d'une function accedee via pointer dans le bootstrap:
   static vector<vector<double> >alllocusStats; // each line= one locus data for all pairs of pops
   size_t nb_loc_migf; // number of loci to be read from the multi matrix file
}

using namespace multimig;


void readMultiMigFile(const char nom_fich_mig[]) {
using namespace datamatrix;
alllocusStats.resize(0); // each line= one locus data for all pairs of pops
vector<double> perlocusStats;
double tmp;
long int Nmiss=0;
char ch;
  ifstream ipt(nom_fich_mig);
string prefiltre;
stringstream filtre(stringstream::in | stringstream::out);;
  while(!ipt.is_open()) {

        noR_cout<<"\n Cannot open file "<<nom_fich_mig<<". Give another input file again: ";

	 string nfs;
	 cin>>nfs;
	 cin.ignore();
	 ipt.clear();
	 ipt.open(nfs.c_str());
	 }   	/*  Data file  */
    ipt.get(ch);
    if (ipt.eof()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"\n This file exists but is empty. I must exit.";
            ipt.clear();
        #else
            cerr<<"\n This file exists but is empty. I must exit.";
            ipt.clear();
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(1, "This file exists but is empty.");
    } else while (ipt.get()!= '\n') continue; //reach the end of last read line
    ipt>>nb_sam_migf;
    data.resize(nb_sam_migf);
    for (vector<vector<long double> >::iterator ii=data.begin();ii!=data.end();ii++) (*ii).resize(nb_sam_migf);
    while (ipt.get()!= '\n') continue; //reach the end of last read line
    ipt>>nb_loc_migf;
    //cout<<nb_sam_migf<<" "<<nb_loc_migf;getchar();
    for (size_t iloc=0;iloc<nb_loc_migf;iloc++) {
         while (ipt.get()!= '\n') continue; //reach the end of last read line
         while (ipt.get()!= '\n') continue; //skips the commentary on the next line
         perlocusStats.resize(0);
    	  for (size_t s1=1; s1<nb_sam_migf; s1++) {
    		for (size_t s2=0; s2<s1; s2++) {
//cout<<"\r"<<s1<<" "<<s2;
// C++ does not seem to read a quiet_NaN even if it has written a quiet_NaN. So we write the string "NaN" rather than a quiet_NaN and read a string
                ipt>>prefiltre;
                if(prefiltre=="NaN") {
     				Nmiss++;

                        if (Nmiss<4) noR_cout<<"\a\n No genetic information for pair "<<s1+1<<" and "<<s2+1;
     				    else if (Nmiss==4) noR_cout<<"\a\n (more pairs without genetic information)";

                    tmp=numeric_limits<double>::quiet_NaN();
                } else {
                    filtre.str(prefiltre);
                    filtre>>tmp;
                    filtre.str("");
                    filtre.clear();
                }
                perlocusStats.push_back(tmp);
//cout<<"\n"<<s1<<" "<<s2<<" "<<tmp;
     		}
//getchar();
    	  }

              if (Nmiss>0) {
          		if (Nmiss>3) noR_cout<<"\a\n  "<<Nmiss<<" pairs without genetic information"<<endl;
          		noR_cout<<"\n\n For pairs of individuals, this will typically occur when such pairs"<<endl;
          		noR_cout<<"have no genotyped locus in common that are polymorphic in the population."<<endl;
          		noR_cout<<"\n The analysis can nevertheless proceed."<<endl;
                  if (pauseGP) { noR_cout<<"\n\n(Return) to continue"<<endl; getchar();}
          	  }


         alllocusStats.push_back(perlocusStats);
//getchar();
        } //end loop on locis
// puis lecture demi matrice geo
//    {ipt.ignore(std::numeric_limits<std::streamsize>::max(), '\n');}
    while (ipt.get()!= '\n') continue; //skipln
//    {ipt.ignore(std::numeric_limits<std::streamsize>::max(), '\n');}
    while (ipt.get()!= '\n') continue; //skipln
    for (size_t s1=1; s1<nb_sam_migf; s1++)
        for (size_t s2=0; s2<s1; s2++) {
            ipt>>data[s1][s2];
//cout<<"\n"<<s1<<" "<<s2<<" "<<data[s1][s2];
            if (ipt.fail()) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::Rcerr<<"\a\n Incomplete distance matrix! Check input file.\n...the program must terminate...\n";
                    ipt.clear();
                #else
                    cerr<<"\a\n Incomplete distance matrix! Check input file.\n...the program must terminate...\n";
                    if (cinGetOnError) cin.get();
                    ipt.clear();
                #endif
                genepop_exit(1, "Incomplete distance matrix!");
            }
        }
    ipt.close();
//cout<<"\aend readMultiMigFile";getchar();
} /*end, readMultiMigFile */

vector<double> ersatz(vector<double> ABCwei) {
/**  wrapper for CIs from several perlocus genetic distance matrices (rather than from num's and denom's)
**/
  using namespace datamatrix;
  using namespace NS_F_est;
   vector<double>t0(3);
   double num,denom;
  /** make the 'data' matrix of (weighted genetic) and geographic distances **/
  size_t pairit=0;
  for (size_t s1=0; s1<nb_sam_migf; s1++)
        for (size_t s2=0; s2<s1; s2++) {
                num=0;denom=0;
                for (size_t it=0;it<nb_loc_migf;it++) {
                    num+=alllocusStats[it][pairit]*ABCwei[it];
                    denom+=ABCwei[it];
                }
                data[s2][s1]=num/denom;
                pairit++;
        } //else cout<<"("<<s1<<","<<s2<<","<<mindist<<") ";
   conversionFst(); //conversion after each reweighting
   // geo distances were read but not converted when the multimig file was read
   if (_first_of_repl) conversionGeo(); // this is called only for the first replicate; note that this operates typeSelection on population (habitats) types
  /** compute a regression **/
  t0=calcwritecorw();
  /** return regression results**/
return t0;
}

double slope(vector<double> ABCwei) {
    vector<double> t0=ersatz(ABCwei);
return(t0[1]);
}

double intercept(vector<double> ABCwei) {
    vector<double> t0=ersatz(ABCwei);
return(t0[0]);
}

double mean(vector<double> ABCwei) {
  vector<double> t0=ersatz(ABCwei);
  return(t0[2]);
}

void initializeMultimig() {
  //nom_fich_multimig.clear();
  alllocusStats.clear();
  //size_t nb_loc_migf = 0;
}

void cleanMultimig() {
  //nom_fich_multimig.clear();
  alllocusStats.clear();
}
