// Contingency table tests and private alleles
/***************************************************************************
@ F. Rousset 2005-2006

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
#include <iomanip>
#include "genepop.h" //template enligne
#include "CT_tests.h"
#include "GenepopS.h"
#include <ctime>
#include <algorithm>
#include "MersenneTwister.h"
#include "tools.h"
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
#endif



//pointer to member function: cf Primer+ p. 1071
void (Cctable::*switchFnPtr)();
static std::vector<double> (Cctable::*testFnPtr)();


using namespace std;

namespace NS_GG { // for genic-genotypic tests

   static size_t allmax;
   static double (Cctable::*statFnPtr)();
}

namespace NS_pairs_for_GeneDivTest {
   static size_t pop1,pop2;
   static bool pairwiseBool;
}

Cctable::~Cctable() {}


Cctable::Cctable(vector<vector<unsigned long int> > table) {
	nb_lig=table.size();
	ctable.resize(nb_lig);
	nb_col=table[0].size();
    for (size_t ii=0;ii<nb_lig;ii++) {
        ctable[ii].resize(nb_col);
        copy(table[ii].begin(),table[ii].end(),ctable[ii].begin());
    }
}

size_t Cctable::maxCellCount() {
    size_t maxcount=0;
     for (unsigned int ii=0;ii<nb_lig;ii++) {
         for (unsigned int jj=0;jj<nb_col;jj++) {
             if (ctable[ii][jj]>maxcount) maxcount=ctable[ii][jj];
         }
     }
return maxcount;
}

void Cctable::filltypesGG(CGenotypes* allGenotypes,char coding){
  //noR_cerr<<"begin filltypesGG "<<coding;

  //genic-genotypic (only complete diploid genotypes
  typesGG.clear();
  typesGG.resize(2);
  allGenotypes->resetIterator();
  ssize_t genotype;  // must be signed type !!
  //int atuer=0;
  while ((genotype = allGenotypes->getNext()) >= 0) { //sur les genos complets seulement
    // 2 lignes suivantes code correct mais le but est d'éviter de refaire ces lectures
    //                          table[0].push_back((*pp1)->getEffective(genotype));
    //                          table[1].push_back((*pp2)->getEffective(genotype));
    typesGG[0].push_back(int(minAllele(genotype,coding)));
    typesGG[1].push_back(int(maxAllele(genotype,coding)));
    // RnoR_cerr<<genotype<<" "<<typesGG[0]<<" "<<typesGG[0];
    // atuer++;
    // if (atuer) genepop_exit(-1,"FIXME")
  }
  //RnoR_cerr<<"end filltypesGG";
}

bool Cctable::purgeZeros(bool genobool) {
using namespace NS_GG;
// vire les lignes et col vides et verif que le résultat a >1 ligs et >1 cols
vector<vector<unsigned long int> >::iterator q;
vector<unsigned long int>::iterator qq;
vector<int>::iterator qqq;
ligmarg.resize(nb_lig,0);colmarg.resize(nb_col,0);total=0;
// vire lignes vides
     for (unsigned int ii=0;ii<nb_lig;ii++) {
         for (unsigned int jj=0;jj<nb_col;jj++) {
             ligmarg[ii]+=ctable[ii][jj]; //(fully informative) genotype counts!
             colmarg[jj]+=ctable[ii][jj];
         }
     }
     q=ctable.begin();
     qq=ligmarg.begin();
     for (unsigned int ii=0;ii<nb_lig;ii++) {
         if (ligmarg[ii]==0) {
            q=ctable.erase(q); //ne bouge pas
            nb_lig--;
         }
         else q++; // ne bouge le pointer que s'il a une ligne non vide car sinon le pointeur est sur lacolonne déplacée
     }
     if (nb_lig<2) return false; // une seule ligne a des données
     //else
// vire cols vides
     for (unsigned int jj=0;jj<nb_col;) {
         if (colmarg[jj]==0) {
//cout<<jj<<" "<<colmarg[jj]<<" "<<nb_col<<" "<<ctable[1].size()<<" "<<types[0].size()<<" "<<types[1].size();getchar();
             for (unsigned int ii=0;ii<nb_lig;ii++) {
                 qq=ctable[ii].begin()+jj;
                 ctable[ii].erase(qq);
             } //for ii<nb_lig
             qq=colmarg.begin()+jj;
             colmarg.erase(qq);
             if (genobool) {
                 qqq=typesGG[0].begin()+jj;
                 typesGG[0].erase(qqq);
                 qqq=typesGG[1].begin()+jj;
                 typesGG[1].erase(qqq);
             } //genobool
             nb_col--;
         } //if colmarg==0
         else jj++;
     } //for jj<nb_col
     if (nb_col<2) return false; // une seule col a des données
// on recalcule les marges de la table nettoyée
     ligmarg.resize(0,0);ligmarg.resize(nb_lig,0);
     colmarg.resize(0,0);colmarg.resize(nb_col,0);
     for (unsigned int ii=0;ii<nb_lig;ii++) {
         for (unsigned int jj=0;jj<nb_col;jj++) {
             ligmarg[ii]+=ctable[ii][jj];
             colmarg[jj]+=ctable[ii][jj];
         }
         total+=ligmarg[ii];
     }
     return true;
}

void Cctable::cumul(double& privfreq,long int& nbpriv,vector<double>& ssize) {
// pour private alleles
    for (unsigned int lig=0;lig<nb_lig;lig++) {
        for (unsigned int col=0;col<nb_col;col++) {
            if (ctable[lig][col]==colmarg[col]) {
//cout<<lig<<" "<<col<<" "<<ligmarg[lig]<<" "<<colmarg[col]<<" "<<ctable[lig][col];getchar();
                privfreq+= double(colmarg[col])/ligmarg[lig];
                nbpriv++;
            }
        }
    }
    ssize[0]+=total;
    ssize[1]+=nb_lig;
}

bool Cctable::verifInfo() {
// finds whether at least one line and at least one column have >1 marginals
int nonzero=0;
bool info=false;
vector<unsigned long int>locligmarg(nb_lig,0),loccolmarg(nb_col,0);
for (size_t ii=0;ii<nb_lig;ii++) {
         for (size_t jj=0;jj<nb_col;jj++) {
             locligmarg[ii]+=ctable[ii][jj];
             loccolmarg[jj]+=ctable[ii][jj];
         }
     }
     for (size_t ii=0;ii<nb_lig;ii++) {
         if (locligmarg[ii]>1) info=true;
         if (locligmarg[ii]>0) nonzero++;
     }
     if (!info) return info;
     if (nonzero<2) return false; // une seule ligne a des données
     //else
     info=false;
     nonzero=0;
     for (size_t jj=0;jj<nb_col;jj++) {
         if (loccolmarg[jj]>1) info=true;
         if (loccolmarg[jj]>0) nonzero++;
     }
     if (nonzero<2) return false; // une seule col a des données
     return info;
}

void Cctable::choix() {
  bool locdebug=false;
   lig1 = size_t(alea.randExc(nb_lig));
  if (locdebug) RnoR_cerr<<lig1<<" ";
   do {
       lig2 = size_t(alea.randExc(nb_lig));
     if (locdebug)  RnoR_cerr<<lig2<<" ";
   } while (lig1==lig2);

   col1 = size_t(alea.randExc(nb_col));
   if (locdebug) RnoR_cerr<<col1<<" ";
   do {
        col2 = size_t(alea.randExc(nb_col));
     if (locdebug) RnoR_cerr<<col2<<" ";
   } while (col1==col2);
   if (locdebug) genepop_exit(-1,"ICI");
}

void Cctable::switchSP(){
     ctable[lig1][col1]-=1;
     ctable[lig2][col2]-=1;
     ctable[lig1][col2]+=1;
     ctable[lig2][col1]+=1;
} // fin de switchSP

void Cctable::switchSP_GG(){
using namespace NS_GG;
     ctable[lig1][col1]-=1;
     ctable[lig2][col2]-=1;
     ctable[lig1][col2]+=1;
     ctable[lig2][col1]+=1;
     atable[lig1][size_t(typesGG[0][col1])]-=1;
     atable[lig1][size_t(typesGG[1][col1])]-=1;
     atable[lig2][size_t(typesGG[0][col2])]-=1;
     atable[lig2][size_t(typesGG[1][col2])]-=1;
     atable[lig1][size_t(typesGG[0][col2])]+=1;
     atable[lig1][size_t(typesGG[1][col2])]+=1;
     atable[lig2][size_t(typesGG[0][col1])]+=1;
     atable[lig2][size_t(typesGG[1][col1])]+=1;
} // fin de switchSP_GG

int Cctable::print(ostream& ost){
     int champ=int(log(0.0001+total)/log(10))+2;
     ost<<setw(6+champ*int(this->nb_col))<<" "<<"  Total\n";
     for (unsigned int ii=0;ii<ctable.size();ii++) {
         ost<<"       ";
         enligne(ctable[ii],ost,champ);
         ost<<"  "<<ligmarg[ii]<<endl;
     }
     ost<<endl;
     ost<<"Total  ";
     enligne(colmarg,ost,champ);
     ost<<"  "<<total<<endl;
return 0;
} // fin de print

vector<double> Cctable::Proba_test() { //Fisher's probability test on the table
unsigned long int prob1,prob2,cnmoy=0;
double prob,rho=0.0,hM = 0.0,hECA = 0.0;
time_t _StartTime,_EndTime;
float _secs;
bool affichIntermBool=false;
long int switches;

    time(&_StartTime);
    for (unsigned long int k=0;k<dem;k++) {
      choix();//choix de 2 lignes et 2 colonnes
      prob1 = ctable[lig1][col1] * ctable[lig2][col2];
      if (prob1 != 0) {
        prob2 = (ctable[lig2][col1] + 1) * (ctable[lig1][col2] + 1);
        prob = double(prob1) / prob2;
        if (!((prob < 1) && (alea()>prob))) { // not (cas ou on ne bouge pas)
           (this->*switchFnPtr)();
           rho += log(prob);
           if (fabs(rho) < DRIFT_PREC) rho = 0;
        }
      }
    }
    switches=0;
  //Realisation des batches
    for (unsigned long int k=0;k<batchnbr;k++) {
       for (unsigned long int cl = 0;cl<batchlgth;cl++) {
           choix();//choix de 2 pops et 2 alleles pour swicher
           prob1 = ctable[lig1][col1] * ctable[lig2][col2];
           if (prob1 != 0) {
            prob2 = (ctable[lig2][col1] + 1) * (ctable[lig1][col2] + 1);
            prob = double(prob1) / prob2;
            if (!((prob < 1) && (alea()>prob))) { // NOT (cas ou on ne bouge pas)
               (this->*switchFnPtr)();
               rho += log(prob);
               switches++;
               if (fabs(rho) < DRIFT_PREC) rho = 0;
            }
           }
           if (rho <= 0) cnmoy += 1;
       }
       hM += double(cnmoy) / batchlgth;
       hECA +=pow(double(cnmoy) / batchlgth,2);
       cnmoy = 0;
       time(&_EndTime);
       _secs=float(_EndTime- _StartTime);
       if (_secs>1) {
          affichIntermBool=true;
          _gotoxy(0,18);
          noR_cout<<"Already "<<k<<" batches done out of "<<batchnbr<<"  "<<endl;
          _StartTime=_EndTime;
       }
    }
    if (affichIntermBool) {
        _gotoxy(0,18);
        noR_cout<<"Already "<<batchnbr<<" batches done out of "<<batchnbr<<"  "<<endl;
    }

    hECA = sqrt((hECA - (hM * hM) / batchnbr) / batchnbr / (batchnbr - 1.0));
    hM = hM / batchnbr;
    vector<double>zut(0);
    zut.push_back(hM);
    zut.push_back(hECA);
    zut.push_back(switches);
return zut;
} // fin Proba_test

double Cctable::calc_Gobs() {
     vector<unsigned long int>locligmarg(nb_lig,0),loccolmarg(nb_col,0);
  unsigned long int loctotal=0,obs;
     double rG=0;
     for (size_t ii=0;ii<nb_lig;ii++) {
         for (size_t jj=0;jj<nb_col;jj++) {
             locligmarg[ii]+=ctable[ii][jj];
             loccolmarg[jj]+=ctable[ii][jj];
         }
     loctotal+=locligmarg[ii];
     }
     expected.resize(nb_lig);
     for (size_t ii=0;ii<nb_lig;ii++) {
         expected[ii].resize(0);
         for (size_t jj=0;jj<nb_col;jj++)
             expected[ii].push_back(double(locligmarg[ii])*loccolmarg[jj]/loctotal);
     }
     for (size_t ii=0;ii<nb_lig;ii++)
         for (size_t jj=0;jj<nb_col;jj++)
             if ((obs=ctable[ii][jj])>0) rG+=log(double(obs)/expected[ii][jj])*obs;
return rG;
}

double Cctable::calc_G_geneDiv_trend() {
    // the return value must be low if div decreases with increasing rank
    // if called with default value cerr, no output to struc_file is produced
    // a cleaner (but longer) code would use a 'nullstream' with appropriate bit sets (ask google)
       double obs,rG=0,tmp;
       vector<double>obsQ(0);

       #ifdef COMPATIBILITYRCPP
           // Rcpp::Rcerr<<"\nGene diversities per population:\n";
       #else
            cerr<<"\nGene diversities per population:\n";
       #endif

        for (size_t ii=0;ii<nb_lig;ii++) {
            obs=0;
            for (size_t jj=0;jj<nb_col;jj++) {
                tmp=ctable[ii][jj]; //ctable vs atable et nb_col vs allmax diffs avec _GG_
                obs+=tmp*(tmp-1.);
            }
            // actually computes Q, not 1-Q !!
            obs/=ligmarg[ii]*(ligmarg[ii]-1); // observed IDENTITY in pop i // _GG_ a des facteurs 2

            #ifdef COMPATIBILITYRCPP
                // Rcpp::Rcerr<<"Pop "<<ii<<": "<<1.-obs<<endl;
            #else
                cerr<<"Pop "<<ii<<": "<<1.-obs<<endl;
            #endif

 //cout<<ii<<" "<<obs;getchar();
 //           if (ii>0) { // compares number of alleles in successive samples: order of samples has a meaning !
 //                // return value 'G' will be compared as G>=Gobs => increases P value
 //               rG+=(obs-preobs);
 //           }
 //           preobs=obs; // stores value for sample, then considers next sample
            obsQ.push_back(obs);
        }
        if (NS_pairs_for_GeneDivTest::pairwiseBool)
               rG=(obsQ[2]-obsQ[1])*
                         (sequenceGeneDivRanks[NS_pairs_for_GeneDivTest::pop2]-sequenceGeneDivRanks[NS_pairs_for_GeneDivTest::pop1]);
        else for (size_t ii=0;ii<nb_lig;ii++) {
           for (size_t jj=ii+1;jj<nb_lig;jj++) {
               rG+=(obsQ[jj]-obsQ[ii])*(sequenceGeneDivRanks[jj]-sequenceGeneDivRanks[ii]);
               /*In order to test that sample 2 has lower div than sample 1, one should have declared
                 sequenceGeneDivRanks[2]>sequenceGeneDivRanks[1] (second of 'higher' rank in diversity = lower diversity)
               */
           }
        }
 //cout<<rG;getchar();
 return rG;
}

double Cctable::calc_G_geneDiv_trend(ostream& ost) { //private option for CA
   // the return value must be low if div decreases with increasing rank
   // if called with default value cerr, no output to struc_file is produced
   // a cleaner (but longer) code would use a 'nullstream' with appropriate bit sets (ask google)
      double obs,rG=0,tmp;
      vector<double>obsQ(0);
      ost<<"\nGene diversities per population:\n";

       for (size_t ii=0;ii<nb_lig;ii++) {
           obs=0;
           for (size_t jj=0;jj<nb_col;jj++) {
               tmp=ctable[ii][jj]; //ctable vs atable et nb_col vs allmax diffs avec _GG_
               obs+=tmp*(tmp-1.);
           }
           // actually computes Q, not 1-Q !!
           obs/=ligmarg[ii]*(ligmarg[ii]-1); // observed IDENTITY in pop i // _GG_ a des facteurs 2
           ost<<"Pop "<<ii<<": "<<1.-obs<<endl;
//cout<<ii<<" "<<obs;getchar();
//           if (ii>0) { // compares number of alleles in successive samples: order of samples has a meaning !
//                // return value 'G' will be compared as G>=Gobs => increases P value
//               rG+=(obs-preobs);
//           }
//           preobs=obs; // stores value for sample, then considers next sample
           obsQ.push_back(obs);
       }
       if (NS_pairs_for_GeneDivTest::pairwiseBool)
              rG=(obsQ[2]-obsQ[1])*
                        (sequenceGeneDivRanks[NS_pairs_for_GeneDivTest::pop2]-sequenceGeneDivRanks[NS_pairs_for_GeneDivTest::pop1]);
       else for (size_t ii=0;ii<nb_lig;ii++) {
          for (size_t jj=ii+1;jj<nb_lig;jj++) {
              rG+=(obsQ[jj]-obsQ[ii])*(sequenceGeneDivRanks[jj]-sequenceGeneDivRanks[ii]);
              /*In order to test that sample 2 has lower div than sample 1, one should have declared
                sequenceGeneDivRanks[2]>sequenceGeneDivRanks[1] (second of 'higher' rank in diversity = lower diversity)
              */
          }
       }
//cout<<rG;getchar();
return rG;
}


double Cctable::calc_GGobs() {
using namespace NS_GG;
     vector<unsigned long int>locligmarg(nb_lig,0);//,colmarg(nb_col,0);
     allmax=size_t(*max_element(typesGG[1].begin(),typesGG[1].end()));
     vector<unsigned long int>all_colmarg(allmax+1,0);
     atable.resize(nb_lig);
     unsigned long int loctotal=0;
     for (size_t ii=0;ii<nb_lig;ii++) {
         atable[ii].resize(0);//remplit *la suite* de zero
         atable[ii].resize(allmax+1); //taille pas economique mais bof
         for (size_t jj=0;jj<nb_col;jj++) {
             locligmarg[ii]+=ctable[ii][jj];
             atable[ii][size_t(typesGG[0][jj])]+=ctable[ii][jj];
             atable[ii][size_t(typesGG[1][jj])]+=ctable[ii][jj];
             all_colmarg[size_t(typesGG[0][jj])]+=ctable[ii][jj];
             all_colmarg[size_t(typesGG[1][jj])]+=ctable[ii][jj];
         }
     loctotal+=locligmarg[ii];
     }
     expected.resize(nb_lig); //effectifs alleliques
     for (size_t ii=0;ii<nb_lig;ii++) {
         expected[ii].resize(0);
         for (size_t jj=0;jj<allmax+1;jj++) {
             expected[ii].push_back(double(locligmarg[ii])*all_colmarg[jj]/loctotal); //ligmarg/total est genotype count/genotype count
         }
     }
return (this->*statFnPtr)();
}

double Cctable::calc_GG() {
using namespace NS_GG;
      int obs;
     double rG=0;
       for (unsigned int ii=0;ii<nb_lig;ii++)
           for (size_t jj=0;jj<allmax+1;jj++)
               if ((obs=atable[ii][jj])>0) rG+=log(double(obs)/expected[ii][jj])*obs;
return rG;
}

double Cctable::calc_alleleNbr_trend() { //private option for EI
using namespace NS_GG;
      int obs,preobs=0;
      double rG=0;
       for (unsigned int ii=0;ii<nb_lig;ii++) {
           obs=0;
           for (size_t jj=0;jj<allmax+1;jj++) {
               if (atable[ii][jj]>0) obs+=1; // simply counts the number of alleles in population ii
           }
           if (ii>0) { // compares number of alleles in successive samples: order of samples has a meaning !
                // return value 'G' will be compared as G>=Gobs => increases P value
                if (preobs<obs) rG-=1; // decreases G => increases P value if number of alleles not decreased
                else if (preobs>obs) rG+=1; // converse
           }
           preobs=obs; // stores value for sample, then considers next sample
       }
return rG;
}

double Cctable::calc_GG_geneDiv_trend() { //private option for EI
    // the return value must be low if div decreases with increasing rank
using namespace NS_GG;
      double obs,rG=0,tmp;
      vector<double>obsQ(0);
//      if(ost!=cerr) ost<<"\nGene diversities per population:\n";
       for (unsigned int ii=0;ii<nb_lig;ii++) {
           obs=0;
           for (size_t jj=0;jj<allmax+1;jj++) {
               tmp=atable[ii][jj];
               obs+=tmp*(tmp-1.);
           }
           // actually computes Q, not 1-Q !!
           // ligmarg (margin of ctable, not atable) is # of genotypes hence 2 ligmarg # of genes
           obs/=2.*ligmarg[ii]*(2*ligmarg[ii]-1); // observed gen div in pop i
//           if(ost!=cerr) ost<<"Pop "<<ii<<": "<<1.-obs<<endl;
//           if (ii>0) { // compares number of alleles in successive samples: order of samples has a meaning !
//                // return value 'G' will be compared as G>=Gobs => increases P value
//               rG+=(obs-preobs);
//           }
//           preobs=obs; // stores value for sample, then considers next sample
           obsQ.push_back(obs);
       }
       if (NS_pairs_for_GeneDivTest::pairwiseBool)
              rG=(obsQ[2]-obsQ[1])*
                        (sequenceGeneDivRanks[NS_pairs_for_GeneDivTest::pop2]-sequenceGeneDivRanks[NS_pairs_for_GeneDivTest::pop1]);
       else for (unsigned int ii=0;ii<nb_lig;ii++) {
          for (unsigned int jj=ii+1;jj<nb_lig;jj++) {
              rG+=(obsQ[jj]-obsQ[ii])*(sequenceGeneDivRanks[jj]-sequenceGeneDivRanks[ii]);
              /*In order to test that sample 2 has lower div than sample 1, one should have declared
                sequenceGeneDivRanks[2]>sequenceGeneDivRanks[1] (second of 'higher' rank in diversity = lower diversity)
              */
          }
       }
//cout<<rG;getchar();
return rG;
}


vector<double> Cctable::GG_test() { //Genic-genotypic G (Likelihood ratio) test on the table
using namespace NS_GG;
unsigned long int prob1,prob2,cnmoy=0;
double prob,rho=0.0,hM = 0.0,hECA = 0.0;
if (alleleNbrTestBool) {  //private option for EI
    statFnPtr=&Cctable::calc_alleleNbr_trend;
} else if (geneDivTestBool) {  //private option for EI
    if (sequenceGeneDivRanks.size()!=nb_lig) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"(!) GeneDivRanks length "<<sequenceGeneDivRanks.size()<<" differs from number of subsamples "<<nb_lig<<". Exiting..."<<struc_file;
        #else
            cerr<<"(!) GeneDivRanks length "<<sequenceGeneDivRanks.size()<<" differs from number of subsamples "<<nb_lig<<". Exiting..."<<struc_file;
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "(!) GeneDivRanks length differs from number of subsamples");
    }
    statFnPtr=&Cctable::calc_GG_geneDiv_trend;
} else {
    statFnPtr=&Cctable::calc_GG;
}
double Gobs=calc_GGobs(); // resizes and fills atable, then calls (*statFnPtr)();
double G=Gobs;
time_t _StartTime,_EndTime;
float _secs;
bool affichIntermBool=false;
long int switches;


  time(&_StartTime);
  for (unsigned long int k=0;k<dem;k++) {
      choix();//choix de 2 lignes et 2 colonnes
      prob1 = ctable[lig1][col1] * ctable[lig2][col2];
      if (prob1 != 0) {
        prob2 = (ctable[lig2][col1] + 1) * (ctable[lig1][col2] + 1);
        prob = double(prob1) / prob2;
        if (!((prob < 1) && (alea()>prob))) { // not (cas ou on ne bouge pas)
           (this->*switchFnPtr)();
           rho += log(prob);
           if (fabs(rho) < DRIFT_PREC) rho = 0;
//           G=calc_GG();
           G=(this->*statFnPtr)();
           if (fabs(G-Gobs) < DRIFT_PREC) G = Gobs;   //G
        }
      }
  }

  switches=0;
  //Realisation des batches
  for (unsigned long int k=0;k<batchnbr;k++) {
       for (unsigned long int cl = 0;cl<batchlgth;cl++) {
           choix();//choix de 2 pops et 2 alleles pour swicher
           prob1 = ctable[lig1][col1] * ctable[lig2][col2];
           if (prob1 != 0) {
            prob2 = (ctable[lig2][col1] + 1) * (ctable[lig1][col2] + 1);
            prob = double(prob1) / prob2;
            if (!((prob < 1) && (alea()>prob))) { // not (cas ou on ne bouge pas)
               (this->*switchFnPtr)();
               rho += log(prob);
               if (fabs(rho) < DRIFT_PREC) rho = 0;
//               G=calc_GG();
               G=(this->*statFnPtr)();
               if (fabs(G-Gobs) < DRIFT_PREC) G = Gobs;   //G
               switches++;
            }
           }
           if (G>=Gobs) cnmoy += 1;
       }
       hM += double(cnmoy) / batchlgth;
       hECA +=pow(double(cnmoy) / batchlgth,2);
       cnmoy = 0;
       time(&_EndTime);
       _secs=float(_EndTime- _StartTime);
       if (_secs>1) {
          affichIntermBool=true;
          _gotoxy(0,18);
         noR_cout<<"Already "<<k<<" batches done out of "<<batchnbr<<"  "<<endl;
          _StartTime=_EndTime;
       }
  }
   if (affichIntermBool) {
      _gotoxy(0,18);
      noR_cout<<"Already "<<batchnbr<<" batches done out of "<<batchnbr<<"  "<<endl;
   }

  hECA = sqrt((hECA - (hM * hM) / batchnbr) / batchnbr / (batchnbr - 1.0));
  hM = hM / batchnbr;
  vector<double>zut(0);
  zut.push_back(hM);
  zut.push_back(hECA);
  zut.push_back(switches);
return zut;
} // fin GG_test


vector<double> Cctable::G_test() { //G (Likelihood ratio) or some other test on the table
  unsigned long int  prob1,prob2,cnmoy=0;
  unsigned long int obs;
  double prob,rho=0.0,hM = 0.0,hECA = 0.0;
  double Gobs=calc_Gobs();
  double G=Gobs;
  time_t _StartTime,_EndTime;
  float _secs;
  bool affichIntermBool=false;
  int switches;
  double stat,obsgeneDiv=0; // used in case test is not G-test
  time(&_StartTime);
  if(alleleNbrTestBool) {
      #ifdef COMPATIBILITYRCPP
          // Rcpp::Rcerr<<"(!) Allele number test on contingency table file not yet implemented. Exiting...";
      #else
        cerr<<"(!) Allele number test on contingency table file not yet implemented. Exiting...";
        if (cinGetOnError) cin.get();
      #endif
     genepop_exit(-1, "(!) Allele number test on contingency table file not yet implemented.");
  } else if (geneDivTestBool) {
     obsgeneDiv=calc_G_geneDiv_trend();
  } //else nothing to do here.
  for (unsigned long int k=0;k<dem;k++) {
      choix();//choix de 2 lignes et 2 colonnes
      prob1 = ctable[lig1][col1] * ctable[lig2][col2];
      if (prob1 != 0) {
        prob2 = (ctable[lig2][col1] + 1) * (ctable[lig1][col2] + 1);
        prob = double(prob1) / prob2;
        if (!((prob < 1) && (alea()>prob))) { // not (cas ou on ne bouge pas)
           if ((obs=ctable[lig1][col1])>0) G-=log(double(obs)/expected[lig1][col1])*obs; //G
           if ((obs=ctable[lig1][col2])>0) G-=log(double(obs)/expected[lig1][col2])*obs; //G
           if ((obs=ctable[lig2][col1])>0) G-=log(double(obs)/expected[lig2][col1])*obs; //G
           if ((obs=ctable[lig2][col2])>0) G-=log(double(obs)/expected[lig2][col2])*obs; //G
           (this->*switchFnPtr)();
           rho += log(prob);
           if (fabs(rho) < DRIFT_PREC) rho = 0;
           if ((obs=ctable[lig1][col1])>0) G+=log(double(obs)/expected[lig1][col1])*obs; //G
           if ((obs=ctable[lig1][col2])>0) G+=log(double(obs)/expected[lig1][col2])*obs; //G
           if ((obs=ctable[lig2][col1])>0) G+=log(double(obs)/expected[lig2][col1])*obs; //G
           if ((obs=ctable[lig2][col2])>0) G+=log(double(obs)/expected[lig2][col2])*obs; //G
           if (fabs(G-Gobs) < DRIFT_PREC) G = Gobs;   //G
        }
      }
  }
  switches=0;
  //Realisation des batches
  for (unsigned long int k=0;k<batchnbr;k++) {
       for (unsigned long int cl = 0;cl<batchlgth;cl++) {
           choix();//choix de 2 pops et 2 alleles pour swicher
           prob1 = ctable[lig1][col1] * ctable[lig2][col2];
           if (prob1 != 0) {
            prob2 = (ctable[lig2][col1] + 1) * (ctable[lig1][col2] + 1);
            prob = double(prob1) / prob2;
            if (!((prob < 1) && (alea()>prob))) { // not (cas ou on ne bouge pas)
               if ((obs=ctable[lig1][col1])>0) G-=log(double(obs)/expected[lig1][col1])*obs; //G
               if ((obs=ctable[lig1][col2])>0) G-=log(double(obs)/expected[lig1][col2])*obs; //G
               if ((obs=ctable[lig2][col1])>0) G-=log(double(obs)/expected[lig2][col1])*obs; //G
               if ((obs=ctable[lig2][col2])>0) G-=log(double(obs)/expected[lig2][col2])*obs; //G
               (this->*switchFnPtr)();
               rho += log(prob);
               if (fabs(rho) < DRIFT_PREC) rho = 0;
               if ((obs=ctable[lig1][col1])>0) G+=log(double(obs)/expected[lig1][col1])*obs; //G
               if ((obs=ctable[lig1][col2])>0) G+=log(double(obs)/expected[lig1][col2])*obs; //G
               if ((obs=ctable[lig2][col1])>0) G+=log(double(obs)/expected[lig2][col1])*obs; //G
               if ((obs=ctable[lig2][col2])>0) G+=log(double(obs)/expected[lig2][col2])*obs; //G
               if (fabs(G-Gobs) < DRIFT_PREC) G = Gobs;   //G
               switches++;
            }
           }
           if(alleleNbrTestBool) {
               #ifdef COMPATIBILITYRCPP
                   // Rcpp::Rcerr<<"(!) Allele number test on contingency table file not yet implemented. Exiting...";
               #else
                    cerr<<"(!) Allele number test on contingency table file not yet implemented. Exiting...";
                    if (cinGetOnError) cin.get();
               #endif
               genepop_exit(-1, "(!) Allele number test on contingency table file not yet implemented.");
           } else if (geneDivTestBool) {
               stat=calc_G_geneDiv_trend();
               if (stat>=obsgeneDiv) cnmoy+=1;
           } else if (G>=Gobs) cnmoy += 1;
       }
       hM += double(cnmoy) / batchlgth;
       hECA +=pow(double(cnmoy) / batchlgth,2);
       cnmoy = 0;
       time(&_EndTime);
       _secs=float(_EndTime- _StartTime);
       if (_secs>1) {
          affichIntermBool=true;
          _gotoxy(0,18);
         noR_cout<<"Already "<<k<<" batches done out of "<<batchnbr<<"  "<<endl;
          _StartTime=_EndTime;
       }
  }
   if (affichIntermBool) {
      _gotoxy(0,18);
      noR_cout<<"Already "<<batchnbr<<" batches done out of "<<batchnbr<<"  "<<endl;
   }

  hECA = sqrt((hECA - (hM * hM) / batchnbr) / batchnbr / (batchnbr - 1.0));
  hM = hM / batchnbr;
  vector<double>zut;
  zut.push_back(hM);
  zut.push_back(hECA);
  zut.push_back(switches);
return zut;
} // fin G_test


// l'autre version:
//------------------------------------------------------------------------------
int crunchLocTable(int statindic,size_t iLoc,ostream& fichier_out,vector<vector<double> > * tabF) {
  //cout<<"debut crunchLocTable()";
  using namespace NS_GG;
  //version option 1 et 3
  // tabF sert a stocker les resultats des tests dans l'option 1
  // pour le cas genotypique on commence par créer une structure cumulée des génotypes (au locus iLoc)
  // ce debut vient directement de genotip2:
  bool genobool=false,pairbool=false;
  if (statindic>2) genobool=true;  // do test on genotypic table whenever possible
  if ((statindic%2)==0) pairbool=true;
  vector<double>testresult(3);
  //int nb_sam=fichier_genepop->pops.size();
  char coding=fichier_genepop->coding[iLoc];
  bool GGbool=(coding>3 && genobool);
  int allelecoding=2+std::max(coding/2,coding%4);   // 2 3 4 6 => 4 5 4 5
  vector<unsigned long int>dummyvec;
  unsigned long int pairit,ligmarg;
  ssize_t genotype;
  int allele;
  vector<vector<unsigned long int> >table,toutable; //table pop X x genotypes counts et types alleliques
  time_t _StartTime,_EndTime;
  float _secs;
  CGenotypes *allGenotypes=NULL; // TOUS les génotypes du locus courant dans le fichier
  vector<CGenotypes *>genopops; // génotypes pour chaque population au locus courant
  vector<CGenotypes *>::iterator pG; // iterateurs sur les effectifs
  vector<CPopulation *>::iterator p,pp; // CPopulation contient les popName()
  vector<CGenotypes *>::iterator p1,p2; // CGenotypes contient toute l'info sauf les popName()
  // construire structure pour tests genotypiques versus genic
  if (GGbool) { //table genotypique
    allGenotypes = new CGenotypes;
    // instanciation des structures de stockage des génotypes
    genopops.resize(fichier_genepop->pops.size());

    // purge de la structure de génotypes globale
    allGenotypes->clear();
    // itérations sur les populations pour le locus courant
    pG = genopops.begin(); // initialisation de l'itérateur sur les populations
    for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
      allGenotypes->fillGenotypes(iLoc, *p,coding); // remplissage de la structure de génotypes cumulés total sample
      (*pG) = new CGenotypes; // instanciation de la structure de stockage de génotypes de la pop courante pour le locus courant
      (*pG)->clear(); // nettoyage
      (*pG)->fillGenotypes(iLoc, *p,coding);
      pG++;
    }
    switchFnPtr=&Cctable::switchSP_GG; //switch qui change les effectifs des tables genic et genotypiq simultanement
    testFnPtr=&Cctable::GG_test;       // pas de possibilité de proba test lé
  } else { // table genique
    switchFnPtr=&Cctable::switchSP; // geniq only
    if (genicProbaTestBool) testFnPtr=&Cctable::Proba_test;
    else testFnPtr=&Cctable::G_test;
  }
  //par pair vs table complete
  if (pairbool) table.resize(2); //par paires
  if (GGbool) fichier_out << "Pop       Genotypes:" << endl; else fichier_out << "Pop       Alleles:" << endl;
  fichier_out << "          ";
  for(int tiret=0; tiret<24; tiret ++){fichier_out << "----------";} // 240 -
  fichier_out << endl<<"          ";
  // ecriture labels colonnes tables: genotypic else genic
  if (GGbool) { // donnees diploides, tables genotypiques
    vector<vector<int> >typesloc(2);
    typesloc[0].resize(0);typesloc[1].resize(0);   //vide le contenu
    allGenotypes->resetIterator();
    while ((genotype = allGenotypes->getNext()) >= 0) { //sur les genos complets seulement
      typesloc[0].push_back(int(minAllele(genotype,coding)));
      typesloc[1].push_back(int(maxAllele(genotype,coding)));
    }
    enligne(typesloc[0],fichier_out,allelecoding);    // def dans genepop.cpp. Ecriture
    fichier_out<<endl;
    fichier_out<<"          ";
    enligne(typesloc[1],fichier_out,allelecoding);
  } else {
    std::vector<int>dummyintvec(0);
    if (coding>3) {// donnes diploides mais table genique
      fichier_genepop->loci[iLoc]->resetgIterator();
      while ((allele = fichier_genepop->loci[iLoc]->getgNext()) >= 0)   dummyintvec.push_back(allele);
    } else { // donnees haploides
      fichier_genepop->loci[iLoc]->resetIterator();
      while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0)   dummyintvec.push_back(allele);
    }
    enligne(dummyintvec,fichier_out,allelecoding);
  }
  fichier_out<<"  Total"<<endl<<endl;
  // valuation et ecriture des lignes de la toutable  (toute, pairbool ou non)
  toutable.resize(0);
  if (GGbool) { // table genotypique
    pG = genopops.begin();
    for(vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
      allGenotypes->resetIterator();
      dummyvec.resize(0);
      ligmarg=0;
      while ((genotype = allGenotypes->getNext()) >= 0) {
        dummyvec.push_back((*pG)->getEffective(genotype));
        ligmarg+=dummyvec.back();
      }
      fichier_out << setw(10) << (*ii)->popName().substr(0,9);   //FER
      enligne(dummyvec,fichier_out,allelecoding);
      fichier_out<<"   "<<ligmarg<<endl;
      toutable.push_back(dummyvec);
      pG++;
    }
    fichier_out<<endl;
    fichier_out<<"Total:    ";
    allGenotypes->resetIterator();
    dummyvec.resize(0);
    ligmarg=0;
    while ((genotype = allGenotypes->getNext()) >= 0) {
      dummyvec.push_back(allGenotypes->getEffective(genotype));
      ligmarg+=dummyvec.back();
    }
    enligne(dummyvec,fichier_out,allelecoding);
    fichier_out<<"   "<<ligmarg<<endl;
  } else { //table genique
    for(vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
      dummyvec.resize(0);
      ligmarg=0;
      if (coding>3) {//table genique sur donnees diploides
        fichier_genepop->loci[iLoc]->resetgIterator();
        while ((allele = fichier_genepop->loci[iLoc]->getgNext()) >= 0) {
          dummyvec.push_back((*ii)->loci[iLoc]->getgEffective(allele));
          ligmarg+=dummyvec.back();
        }
      } else { // donnees haploides
        fichier_genepop->loci[iLoc]->resetIterator();
        while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0) {
          dummyvec.push_back((*ii)->loci[iLoc]->getEffective(allele));
          ligmarg+=dummyvec.back();
        }
      }
      fichier_out << setw(10) << (*ii)->popName().substr(0,9);   //FER
      enligne(dummyvec,fichier_out,allelecoding);
      fichier_out<<"   "<<ligmarg<<endl;
      toutable.push_back(dummyvec);
    }
    fichier_out<<endl;
    fichier_out<<"Total:    ";
    dummyvec.resize(0);
    ligmarg=0;
    if (coding>3) {//table genique sur donnees diploides
      fichier_genepop->loci[iLoc]->resetgIterator();
      while ((allele = fichier_genepop->loci[iLoc]->getgNext()) >= 0) {
        dummyvec.push_back(fichier_genepop->loci[iLoc]->getgEffective(allele));
        ligmarg+=dummyvec.back();
      }
    } else { // donnees haploides
      fichier_genepop->loci[iLoc]->resetIterator();
      while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0) {
        dummyvec.push_back(fichier_genepop->loci[iLoc]->getEffective(allele));
        ligmarg+=dummyvec.back();
      }
    }
    enligne(dummyvec,fichier_out,allelecoding);
    fichier_out<<"   "<<ligmarg<<endl;
  } // fin genobool or not


  if (pairbool) {
    NS_pairs_for_GeneDivTest::pairwiseBool=true;
    fichier_out<<endl;
    fichier_out<<"Locus        Population pair        P-Value  S.E.     Switches"<<endl;
    fichier_out<<"-----------  ---------------------  -------  -------  --------"<<endl;
    //    	p = fichier_genepop->pops.begin();
    pairit=0;
    time(&_StartTime);
    long int nbPairs=toutable.size()*(toutable.size()-1)/2;
    // attention ordre des paires rétabli identiq a v3.4 le 25/10/06
    //        for(vector<vector<int> >::iterator p1=toutable.begin();p1<toutable.end();p1++) {
    //            pp = p+1;
    NS_pairs_for_GeneDivTest::pop2=1;
    pp = fichier_genepop->pops.begin()+1; //pointeur sur la pop d'index superieur
    for(vector<vector<unsigned long int> >::iterator pp2=toutable.begin()+1;pp2<toutable.end();pp2++) {
      NS_pairs_for_GeneDivTest::pop1=0;
      p = fichier_genepop->pops.begin(); //pointeur sur la pop d'index inferieur
      for(vector<vector<unsigned long int> >::iterator pp1=toutable.begin();pp1!=pp2;pp1++) {
        table.resize(0);
        fichier_out << setw(13)<<fichier_genepop->loci[iLoc]->locusName.substr(0,12)<<setw(9);
        fichier_out << (*pp)->popName().substr(0,9)<<" & "<<setw(11)<< (*p)->popName().substr(0,9);
        table.push_back(*pp1);
        table.push_back(*pp2);
        Cctable ctable(table); // structure pour test CT
        if (GGbool) ctable.filltypesGG(allGenotypes,coding);
        if (!ctable.purgeZeros(GGbool)) {
          fichier_out<<"No table"<<endl;
          (*tabF)[pairit][iLoc]=-1;
        } else if (!ctable.verifInfo()) {
          fichier_out<<"No information"<<endl;
          (*tabF)[pairit][iLoc]=-1;
        } else {
          testresult=(ctable.*testFnPtr)();    /// TEST HERE
          fichier_out<<fixed<<setw(7)<<testresult[0]<<"  "<<setw(7)<<testresult[1];
          fichier_out<<"  "<<setw(8)<<right<<int(testresult[2]+0.5)<<left<<endl;
          (*tabF)[pairit][iLoc]=testresult[0];
        } //verifInfo
        pairit++;
        time(&_EndTime);
        _secs=float(_EndTime- _StartTime);
        if (_secs>1) {
          _gotoxy(0,18);
          noR_cout<<"Already "<<pairit<<" pairs analyzed out of "<<nbPairs<<"  ";
          _StartTime=_EndTime;
        }
        //      		    pp++;
        p++;
        NS_pairs_for_GeneDivTest::pop1++;
      }
      //        	p++;
      pp++;
      NS_pairs_for_GeneDivTest::pop2++;
    } // fin boucle pp1
  } else { // NOT pairbool
    NS_pairs_for_GeneDivTest::pairwiseBool=false;
    Cctable ctable(toutable); // structure pour test CT
    if (GGbool) ctable.filltypesGG(allGenotypes,coding);
    fichier_out<<endl;
    if (!ctable.purgeZeros(GGbool)) {
      fichier_out<<"No table"<<endl;
      (*tabF)[0][iLoc]=-1;
    } else if (!ctable.verifInfo()) {
      fichier_out<<"No information"<<endl;
      (*tabF)[0][iLoc]=-1;
    } else {
      if (statindic<5) {//test CT
        testresult=(ctable.*testFnPtr)();
        fichier_out<<"P-value = "<<testresult[0]<<"    S.E. = "<<testresult[1];
        fichier_out<<" ("<<int(testresult[2]+0.5)<<" switches)"<<left<<endl;
        (*tabF)[0][iLoc]=testresult[0];
      } else {//
      }
    } //verifInfo
  }
  if (GGbool) { //genotypes

    for(pG = genopops.begin(); pG != genopops.end(); pG++) delete (*pG);
    delete allGenotypes;
  }
  return 0;
}


//------------------------- replacement for STRUC.BAT --------------------
//------------------------------------------------------------------------
int struc() {
string ligne;
int nLig,nCol;
unsigned long int count;
vector<unsigned long int>dummyvec;
vector<double>testresult(3);
vector<vector<unsigned long int> >table;
time_t begTime,endTime;
clock_t begClock,endClock;
float _secs;

    effacer_ecran();
    afficher_version();
    if(alleleNbrTestBool) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"(!) From struc(): Allele number test on contingency table file not yet implemented. Exiting...";
        #else
            cerr<<"(!) From struc(): Allele number test on contingency table file not yet implemented. Exiting...";
            if (cinGetOnError) cin.get();
        #endif
       genepop_exit(-1, "(!) From struc(): Allele number test on contingency table file not yet implemented.");
    } else if (geneDivTestBool) {
        noR_cout<<"\nTesting decrease in gene diversity\nInput file: "<<struc_file<<endl;
    } else {
        noR_cout<<"\nTesting independence in contingency tables\nInput file: "<<struc_file<<endl;
    }
    ifstream inputf(struc_file.c_str());
    if  (!inputf.is_open()){
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"could not open file "<<struc_file;
        #else
            cerr<<"could not open file "<<struc_file;
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "could not open file ");
    }
    getline(inputf,ligne); // passe sur la ligne de commentaires
    if (inputf.eof()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"Premature end of "<<struc_file<<" file.\nCheck first line of input";
        #else
            cerr<<"Premature end of "<<struc_file<<" file.\nCheck first line of input";
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "Premature end of file. Check first line of input");
    }
    else if (inputf.fail()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"Cannot read "<<struc_file<<" correctly.\nCheck first line of input";
        #else
            cerr<<"Cannot read "<<struc_file<<" correctly.\nCheck first line of input";
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "Cannot read correctly. Check first line of input");
    }
    inputf>>nLig;
    while (!inputf.eof() && !(isdigit(inputf.peek()))) inputf.get();
    if (inputf.eof()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"Premature end of "<<struc_file<<" file.\nCheck second line of input";
        #else
            cerr<<"Premature end of "<<struc_file<<" file.\nCheck second line of input";
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "Premature end of file. Check second line of input"); }
    else if (inputf.fail()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"Cannot read "<<struc_file<<" correctly.\nCheck second line of input";
        #else
            cerr<<"Cannot read "<<struc_file<<" correctly.\nCheck second line of input";
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "Cannot read correctly. Check second line of input"); }
    inputf>>nCol;
    for(int i=0;i<nLig;i++){
       dummyvec.resize(0);
       for(int j=0;j<nCol;j++){
           while (!inputf.eof() && !(isdigit(inputf.peek()))) inputf.get();
           if (inputf.eof()) {
               #ifdef COMPATIBILITYRCPP
                   // Rcpp::Rcerr<<"Premature end of "<<struc_file<<" file.\nCheck line "<<3+i<<" of input";
               #else
                    cerr<<"Premature end of "<<struc_file<<" file.\nCheck line "<<3+i<<" of input";
                    if (cinGetOnError) cin.get();
               #endif
               genepop_exit(-1, "Premature end of file. Check line of input"); }
           else if (inputf.fail()) {
               #ifdef COMPATIBILITYRCPP
                   // Rcpp::Rcerr<<"Cannot read "<<struc_file<<" correctly.\nCheck line "<<3+i<<" of input";
               #else
                    cerr<<"Cannot read "<<struc_file<<" correctly.\nCheck line "<<3+i<<" of input";
                    if (cinGetOnError) cin.get();
               #endif
               genepop_exit(-1, "Cannot read correctly. Check line of input"); }
           inputf>>count;
           dummyvec.push_back(count);
       }
       table.push_back(dummyvec);
    }
    inputf.close();
    ofstream outf(struc_file.c_str(),ios::app);
    if  (!outf.is_open()){
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"could not open file "<<struc_file;
        #else
            cerr<<"could not open file "<<struc_file;
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(-1, "could not open file ");
    }
    Cctable ctable(table); // structure pour test CT
    if (!ctable.purgeZeros(false)) {
       outf<<"\nNo data!"<<endl;

        noR_cout<<"\nNo data!"<<endl;

    } else if (!ctable.verifInfo()) {
       outf<<"\nNo statistical information (too low marginal counts)"<<endl;
        noR_cout<<"\nNo statistical information (too low marginal counts)"<<endl;
    } else {
       outf<<"\n\nGenepop"<< getSetting("version")<<"\n\n";
       outf<<"\nTest of association in contingency table:"<<endl;
       ctable.print(outf);
       switchFnPtr=&Cctable::switchSP; // geniq only
       if(alleleNbrTestBool) {
           #ifdef COMPATIBILITYRCPP
               // Rcpp::Rcerr<<"(!) From struc(): Allele number test on contingency table file not yet implemented. Exiting...";
           #else
                cerr<<"(!) From struc(): Allele number test on contingency table file not yet implemented. Exiting...";
                if (cinGetOnError) cin.get();
           #endif
          genepop_exit(-1, "(!) From struc(): Allele number test on contingency table file not yet implemented.");
       } else if (geneDivTestBool) {
          testFnPtr=&Cctable::G_test;
          if (sequenceGeneDivRanks.size()!=ctable.get_nb_lig()) {
              #ifdef COMPATIBILITYRCPP
                  // Rcpp::Rcerr<<"(!) GeneDivRanks length "<<sequenceGeneDivRanks.size()<<" differs from number of subsamples ";
                  // Rcpp::Rcerr<<ctable.get_nb_lig()<<". Exiting..."<<struc_file;
              #else
                cerr<<"(!) GeneDivRanks length "<<sequenceGeneDivRanks.size()<<" differs from number of subsamples ";
                cerr<<ctable.get_nb_lig()<<". Exiting..."<<struc_file;
                if (cinGetOnError) cin.get();
              #endif
             genepop_exit(-1, "(!) GeneDivRanks length differs from number of subsamples ");
          }
          ctable.calc_G_geneDiv_trend(outf);
          outf<<"\nExact test using gene diversities;"<<endl;
          noR_cout<<"\nExact test using gene diversities;"<<endl;
       } else if (genicProbaTestBool) {
          testFnPtr=&Cctable::Proba_test;
          outf<<"\nFisher's exact probability test;"<<endl;

          noR_cout<<"\nFisher's exact probability test;"<<endl;

       } else {
          testFnPtr=&Cctable::G_test;
          outf<<"\nExact G (\"Likelihood ratio\") test;"<<endl;
          noR_cout<<"\nExact G (\"Likelihood ratio\") test;"<<endl;
       }
       outf<<"\nEstimation of Pvalue by Markov chain method;\n";
       set_MC_parameters(true);
       outf<<"Markov chain parameters"<<endl<<"\tDememorisation       : "<<dem<<endl<<"\tBatches              : "<<batchnbr<<endl<<"\tIterations per batch : "<<batchlgth<<endl;

        noR_cout<<"\nEstimation of Pvalue by Markov chain method;\n";
        noR_cout<<"Markov chain parameters"<<endl<<"\tDememorisation       : "<<dem<<endl<<"\tBatches              : "<<batchnbr<<endl<<"\tIterations per batch : "<<batchlgth<<endl;
	   time(&begTime);
   	   begClock=clock();
       testresult=(ctable.*testFnPtr)();
   	   endClock=clock();
	   time(&endTime);
       if ((_secs=float(endTime- begTime))< 2147) _secs=float(endClock- begClock)/CLOCKS_PER_SEC;
       outf<<"\nP value="<<setw(9)<<testresult[0]<<"  S.E="<<setw(8)<<fixed<<testresult[1];
       outf<<" ("<<int(testresult[2]+0.5)<<" switches)"<<endl;
       outf.precision(3); //with <<fixed, ensures 3 digits after decimal point
       outf<<"Computation time= "<<fixed<<_secs<<" s.\nNormal ending.\n";

        noR_cout<<"\nP value="<<setw(9)<<testresult[0]<<"  S.E="<<setw(8)<<fixed<<testresult[1];
        noR_cout<<" ("<<int(testresult[2]+0.5)<<" switches)"<<endl;
        noR_cout.precision(3); //with <<fixed, ensures 3 digits after decimal point
        noR_cout<<"Computation time= "<<fixed<<_secs<<" s.\nNormal ending.\n";


    } //verifInfo
    outf.close();
    if (pauseGP) {
        noR_cout<<"(Return) to continue"<<endl; getchar();

    }
return 0;
}

void Cctable::cleanVar(){
    this->typesGG.clear();
    this->atable.clear();
}

void initializeCTtests() {
  //&Cctable::cleanVar;
  NS_GG::allmax = 0;
}
/*
void cleanCTtests() {
    &Cctable::cleanVar;
}*/
