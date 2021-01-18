/***************************************************************************
@ F. Rousset 2005-2006
distantly derived from earlier code by F.R.
distributed with previous versions of Genepop 1997-2003

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
#include "MersenneTwister.h"
#include <sstream>
#include <map>
#include <algorithm>
#include <cerrno>
#include "GenepopS.h" // Zegenepopsound...
#include "genepop.h"
#include "settings.h"
#include "F_est.h"
#include "myutils.h"
#include "bootstrap.h"
#include "tools.h"
#include <errno.h>
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif


#define dsqr(x) ((x)*(x));

using namespace std;


string  _logdist=""; //indic loga transfo
bool _a_stat=true,_e_stat=false; // il en faut pour a et e au cas ou une 3e viendrait...
double mindist=-1;
double maxdist=numeric_limits<double>::max();
/// for ad hoc geographic distance file
bool geoDistFromGeoFile=false;
string geoDistFile;

// vector<vector<double> >& get_geoHalfMat();
// vector<vector<double> >& get_geoHalfMat() {
//   static vector<vector<double> > geoHalfMat;
//   return(geoHalfMat);
// }

string statname=""; // statname pour affichages statInfo pour algo
//vector<double>diversities; //values computed in option 5.2 used in 6.5/6.6

namespace datamatrix {
size_t nb_sam_migf; // number of samples to be read from the multi matrix file
   vector<vector<long double> >data;
}



namespace NS_F_est {
static unsigned long int deuxsamp;
static size_t	global_geno_nbr;
    static std::size_t F_est_locIt;
    std::size_t global_pop_it,
                global_pop2_it; /* for pairs*/
    static size_t human_loc_it, /*from 1.*/
		              nonvides; /*nb of pops with data for a given locus */
		static int indic12all;  // analyse par pop ou par paire de pop ou totale
    static bool identity; //indic allele size/identity
    static bool _perf,_indiv,_pmat;
    bool _first_of_repl;
    size_t nb_sam;
    size_t nb_pair_sam;
    size_t nb_locus;
    static size_t nb_sam_sel;
    static size_t nb_pair_sam_sel;

#ifdef BIGTABLES
    static struct MStype<MSreal> MSbin; // structure definie dans le header

//    MStype *houlaMS; //serait toutes les loc*paires dans un vecteur
    static MStype<MSreal> ***MStable,*MStableit; //serait toutes les loc*paires dans un vecteur

//l'alternative est de recalcluler lecture_Xi_Xj_pmoy for each bootstrap iteration
//void relit_rho();
#endif

    static double *loc_MSG;
    //#ifdef (e stat)
    static double *loc_MSP;
    static double *Qpp;
    static double *QriQrj;
    static long int *nnn;
    static double ***houla;
    static double sumQbij,sumQriQrj,sumQpp;
    //#endif
    double MSp2P,MSg2Pw;
    static double denom_pot;//
    double MSi2P,MSi2Pw; // for analyses by pops


    static size_t **tabM; //genotype counts [pop,genotype]
    static size_t	**tabCode; //alleles de chque genotype: size_t bc mostly used as index.
    static double SSp,SSi,SSg,MSp,MSi,MSg;
    static size_t indices[2];
    //static ofstream f_mig;
    static size_t popi,popj;
//    vector<vector<long double> >data;
    //static vector<vector<long double> >selecteddata;
    //static vector<vector<long double> >indx;
    //static vector<vector<double> >sfreqs;
    //static unsigned long int Ntot;
    static double Nc,mindistorlogdist,maxdistorlogdist;
    static size_t maxallname; //  locus
    //static vector<unsigned long int>scgsp;//  sum_counts_geno_sur_pops /*effectif par pop*/

}  //namespace

namespace NS_FFF_slmt {
 double SSiTot,SSgTot,MSg2P;
}


namespace NS_tailles {
  //static vector<double>tailleMoy;
  static double tailleMoyTot;
}

int tailleOfType(size_t type);
int tailleOfType(size_t type) {
using namespace NS_F_est;
    map<int,int>::iterator ptr;
    // taille is a map<int,int> hence .find(key_type=int)
    if (taille.size()<=F_est_locIt || (ptr=taille[F_est_locIt].find(int(type)))==taille[F_est_locIt].end())
       return int(type); // taille=type si element pas trouve dans la map type->taille
                    // (avec l'avantage que 0 renvoit 0: cf option 5.3 en haploide)
    else
       return ptr->second;
}



int set_options(bool locperf,bool indiv,bool identBool) {
using namespace NS_F_est;
	_perf=locperf;
	_indiv=indiv;
	identity=identBool;
return 0;
}

int set_phylipmatrix(bool pmat) {
using namespace NS_F_est;
	_pmat=pmat;
return 0;
}


int set_ptrs() {
using namespace NS_F_est;
		loc_MSG=new double[nb_locus];
		if (_e_stat) {
//cout<<"alloc loc_MSP["<<nb_locus;getchar();
			loc_MSP=new double[nb_locus];
			Qpp=new double[nb_locus];
			QriQrj=new double[nb_locus];
			nnn=new long int[nb_locus];
		}
//else cout<<"_e_stat=false !!";getchar();

		houla=new double**[nb_locus];
	//	houla.resize(nb_locus+1);
		for(size_t locit=0;locit<nb_locus; locit++) {
	// popi va de 1 e nb_sam-1 maxi	=> nbsam-1 values
	//		houla[locit].resize(nb_sam);
			houla[locit]=new double*[nb_sam-1];
	// popj va de 2 e nb_sam=>nb_sam-1 values
			for(size_t iit=0;iit<nb_sam-1; iit++) {//houla[locit][iit].resize(nb_sam+1);
				houla[locit][iit]=new double[nb_sam-1];
			}
		}
//cout<<"ici alloc houlaMS"<<nb_locus*nb_pair_sam+1;getchar();
//		try{houlaMS=new MStype[nb_locus*nb_pair_sam+1];}
//		catch(std::bad_alloc) {
//           cerr<<"No memory available for allocation in set_first_of_repl_ptrs(.)"<<endl;
//  		   if (cinGetOnError) cin.get();
//		   genepop_exit(-1);
//		}
		MStable=new MStype<MSreal>**[nb_locus];
		for(size_t locit=0;locit<nb_locus; locit++) {
           MStable[locit]=new MStype<MSreal>*[nb_sam-1];
		   for(size_t iit=0;iit<nb_sam-1; iit++) {//houla[locit][iit].resize(nb_sam+1);
				MStable[locit][iit]=new MStype<MSreal>[iit+1];
		   }
        }
return 0;
}

int delete_ptrs() {// those allocated by set_first_of_repl_ptrs
using namespace NS_F_est;
	for(size_t locit=0;locit<nb_locus; locit++) {
		for(size_t iit=0;iit<nb_sam-1; iit++) delete[] houla[locit][iit];
		delete[] houla[locit];
	}
	delete[] houla;
//cout<<"ici delete houlaMS"<<nb_locus*nb_pair_sam+1;getchar();
//    delete[] houlaMS;
  for(size_t locit=0;locit<nb_locus; locit++) {
		for(size_t iit=0;iit<nb_sam-1; iit++) delete[] MStable[locit][iit];
		delete[] MStable[locit];
	}
    delete[] MStable;
    delete[] loc_MSG;
    if (_e_stat) {
//cout<<"delete loc_MSP";getchar();
    	delete[] loc_MSP;
    	delete[] Qpp;
    	delete[] QriQrj;
    	delete[] nnn;
    }
return 0;
}

int genotip2() { //as of 28/03/2007 this is only called by isolde_etc()
// cree des fichiers pour tous les locis, haploides et diploides
//voir genotip2new() dans CT_tests.cpp   !!! -- n'existe pas. struc(...), peut etre ?
//ecriture tous genotypes pour differentes pop et un locus dans fichier LOCUSx
using namespace NS_F_est;
	vector<vector<int> >typestip2(2);
	vector<unsigned long int>dummyvec;
    vector<CPopulation *>::iterator p; // iterateur sur les populations
	CGenotypes *allGenotypes; // TOUS les genotypes du locus courant dans le fichier
	vector<CGenotypes *>genopops; // genotypes pour chaque population au locus courant
	vector<CGenotypes *>::iterator pG; // iterateur sur les genotypes
	string fileName; // LOCUS%d
	string line;
	char coding;
    int allelecoding;   // 2 3 4 6 => 4 5 4 5
	int allele;
	ssize_t genotype;
  unsigned long int ligmarg=0;

	// instanciation des structures de stockage des genotypes
	allGenotypes = new CGenotypes;
	genopops.resize(fichier_genepop->pops.size());

	// on traite locus par locus (boucle externe) ; 1 fichier par locus
	for (F_est_locIt = 0; F_est_locIt < fichier_genepop->loci.size(); F_est_locIt ++) {
        coding=fichier_genepop->coding[F_est_locIt];
        allelecoding=2+std::max(coding/2,coding%4);   // 2 3 4 6 => 4 5 4 5
		// purge de la structure de genotypes globale
		allGenotypes->clear();
		// iterations sur les populations pour le locus courant
		pG = genopops.begin(); // initialisation de l'iterateur sur les populations
		for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
			allGenotypes->fillGenotypes(F_est_locIt, *p,coding); // remplissage de la structure de genotypes globale
			(*pG) = new CGenotypes; // instanciation de la structure de stockage de genotypes de la pop courante pour le locus courant
			(*pG)->clear(); // nettoyage
			(*pG)->fillGenotypes(F_est_locIt, *p,coding);
			pG++;
		}

		// creation et ouverture du fichier de sortie pour le locus courant
		{
		  std::stringstream stst;
		  stst<<"LOCUS"<<F_est_locIt+1;
		  fileName = stst.str();
		}
		ofstream fichier_out(fileName.c_str(), ios::out);
		// ecritures, l'horreur de reproduire les sorties basic
		fichier_out << "From file " << fichier_genepop->fileName <<"  (Locus: " << fichier_genepop->loci[F_est_locIt]->locusName << ")" << endl;
		if (coding>3) { // ie diploid data
      /// first the table that is to be read by Genepop
      // number of pops and number of genotypes
		  fichier_out << setw(4) << fichier_genepop->pops.size() << " " << setw(3) << allGenotypes->getNumber() << " " << endl;
		  fichier_out.setf(ios_base::left,ios_base::adjustfield); // calage e gauche
		  // a contingency table: loop over pops, row = vector of genotype counts
		  for(pG = genopops.begin(); pG != genopops.end(); pG++) {
		    dummyvec.resize(0);
		    line = "";
		    (*pG)->printValues(allGenotypes, &line, 2, 0); //formatage rustique mais OK pour la lecture machine
		    fichier_out << line << " " << endl;
		  }
		} else { // haploid data
		  fichier_out << setw(4) << fichier_genepop->pops.size() << " " << setw(3) << fichier_genepop->loci[F_est_locIt]->getNumber() << " " << endl;
		  fichier_out.setf(ios_base::left,ios_base::adjustfield); // calage e gauche
		  for(vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
		    dummyvec.resize(0);
		    fichier_genepop->loci[F_est_locIt]->resetIterator();
		    while ((allele = fichier_genepop->loci[F_est_locIt]->getNext()) >= 0) {
		      dummyvec.push_back((*ii)->loci[F_est_locIt]->getEffective(allele));
		      ligmarg+=dummyvec.back();
		    }
		    enligne(dummyvec,fichier_out,allelecoding);
		    fichier_out<<endl;
		  }
		}
		// la suite est un formatage propre pour lecture humaine. Simplifie e partir de crunchloclable
		fichier_out << endl;
		fichier_out << "          Genotypes:" << endl;
		fichier_out << "          ";
		for(int tiret=0; tiret<24; tiret ++){fichier_out << "----------";} // 240 -
		fichier_out << endl;
// noms des genotypes
    	if (coding>3) {
    	  typestip2[0].resize(0);typestip2[1].resize(0);   //vide le contenu
    	    allGenotypes->resetIterator();
            while ((genotype = allGenotypes->getNext()) >= 0) { //sur les genos complets seulement
              typestip2[0].push_back(int(minAllele(genotype,fichier_genepop->coding[F_est_locIt])));
              typestip2[1].push_back(int(maxAllele(genotype,fichier_genepop->coding[F_est_locIt])));
            }
    	    fichier_out<<"          ";
            enligne(typestip2[0],fichier_out,allelecoding);    // template dans genepop.h. Ecriture destypes alleliques
            fichier_out<<endl;
    	    fichier_out<<"Pops:     ";
            enligne(typestip2[1],fichier_out,allelecoding);
        } else {
            fichier_genepop->loci[F_est_locIt]->resetIterator();
            dummyvec.resize(0);
            while ((allele = fichier_genepop->loci[F_est_locIt]->getNext()) >= 0) dummyvec.push_back(size_t(allele));
            enligne(dummyvec,fichier_out,allelecoding);
        }
		fichier_out << "   All\n----" << endl;
    	if (coding>3) { // donnees diploides
          	pG = genopops.begin();
            for(vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
                fichier_out << setw(10) << (*ii)->popName().substr(0,9);   //FER
                allGenotypes->resetIterator();
                dummyvec.resize(0);
                ligmarg=0;
                while ((genotype = allGenotypes->getNext()) >= 0) {
    //cout<<genotype;getchar();
                      dummyvec.push_back((*pG)->getEffective(genotype));
                      ligmarg+=dummyvec.back();
                }
                enligne(dummyvec,fichier_out,allelecoding);
                fichier_out<<"   "<<ligmarg<<endl;
//                toutable.push_back(dummyvec);
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
        } else { //donnees haploides
            for(vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
                fichier_out << setw(10) << (*ii)->popName().substr(0,9);   //FER
                dummyvec.resize(0);
                ligmarg=0;
                fichier_genepop->loci[F_est_locIt]->resetIterator();
                while ((allele = fichier_genepop->loci[F_est_locIt]->getNext()) >= 0) {
                      dummyvec.push_back((*ii)->loci[F_est_locIt]->getEffective(allele));
                      ligmarg+=dummyvec.back();
                }
                enligne(dummyvec,fichier_out,allelecoding);
                fichier_out<<"   "<<ligmarg<<endl;
//                toutable.push_back(dummyvec);
            }
            fichier_out<<endl;
            fichier_out<<"Total:    ";
            dummyvec.resize(0);
            ligmarg=0;
            fichier_genepop->loci[F_est_locIt]->resetIterator();
            while ((allele = fichier_genepop->loci[F_est_locIt]->getNext()) >= 0) {
                  dummyvec.push_back(fichier_genepop->loci[F_est_locIt]->getEffective(allele));
                  ligmarg+=dummyvec.back();
            }
            enligne(dummyvec,fichier_out,allelecoding);
            fichier_out<<"   "<<ligmarg<<endl;
        } // fin coding
/*	    pG = genepops.begin();
		for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
			fichier_out << setw(10) << (*p)->popName().substr(0,9);   //FER
			line = "";
			(*pG)->printValues(allGenotypes, &line, 3, -1);
			fichier_out << line << "   " << (*pG)->getSum() << endl;
			pG++;
		}
		fichier_out << endl;
		fichier_out << setw(10) << "All:";
		line = "";
		allGenotypes->printValues(allGenotypes, &line, 3, -1);
		fichier_out << line << "  " << allGenotypes->getSum() << endl;*/
		fichier_out.close();
	    for(pG = genopops.begin(); pG != genopops.end(); pG++) delete (*pG);
	} // bcle sur loc
	delete allGenotypes;
//cout<< " fin genotip2 "<<flush;getchar();
	return 0;
}


void lecture_floc(){
//lecture des fichiers LOCUSn crees par genotip2 (PAS l'option 5: voir tabFtoTabM)
// loci diploides ou haploides
using namespace NS_F_est;
   string d,dummy_s;
   stringstream stst(stringstream::in | stringstream::out);
   stst<<"LOCUS";
//cout<<stst.str()<<" "<<bidon;getchar();
   stst<<F_est_locIt+1; // (evite le 0)
//cout<<stst.str();getchar();
//noter que nb_sam_lu va prendre la valeur fichier_genepop->pops.size()
//ecrite par genotip2 dans LOCUSn
size_t nb_sam_lu;
   ifstream flocusIn(stst.str().c_str(), ios::in);
   if(!flocusIn.is_open()){
       #ifdef COMPATIBILITYRCPP
           // Rcpp::R
       #else
            cerr<<"Error while reading file "<<stst.str()<<endl;
            if (cinGetOnError) cin.get();
       #endif

    	genepop_exit(-1, "Error while reading file ");
   }
//   flocusIn.getline(d, 255);
   getline(flocusIn,d);  // passe sur ligne de commentaires
   flocusIn>>nb_sam_lu>>global_geno_nbr;
//cout<<nb_sam_lu<<" "<<global_geno_nbr;getchar();
   if(nb_sam_lu == 0 || global_geno_nbr == 0){
   }else{
      // dimensionnement des tab
      tabM=new size_t*[nb_sam_lu];
      for(size_t pop = 0; pop < nb_sam_lu; pop++) {
   		tabM[pop]=new size_t[global_geno_nbr];
      }

      tabCode=new size_t*[global_geno_nbr];
      for(size_t g = 0; g < global_geno_nbr; g++) tabCode[g]=new size_t[2];
      // lecture e partir ligne 3 du fichier : effectif genotypiques (ou alleliques...) par pop
      for(size_t pop = 0; pop < nb_sam_lu; pop++){
         for(size_t g = 0; g < global_geno_nbr; g++){
            flocusIn>>tabM[pop][g];
//cout<<pop<<" "<<g<<tabM[pop][g];getchar();
         } // next g
      } // next pop
      flocusIn>>dummy_s; // passe sur "genotypes:"
      flocusIn>>dummy_s; // passe sur "----------"
      for(size_t g = 0; g < global_geno_nbr; g++){ //sert dans calc_sfreqs;Nc(), meme en haplo
         flocusIn>>tabCode[g][0]; // lit numero 1e allele du genotype
//cout<<tabCode[g][0];getchar();
      } // next j
      if (fichier_genepop->coding[F_est_locIt]>3) { // si le locus est diploide
          flocusIn>>dummy_s; // passe sur "pop:"
          for(size_t g = 0; g < global_geno_nbr; g++){
             flocusIn>>tabCode[g][1]; // lit numero 2e allele du genotype
          } // next g
      } else for(size_t g = 0; g < global_geno_nbr; g++){
             // pas vrai par defaut!
             tabCode[g][1]=0; // sert de test diploidie dans calculSSetMS()
      } // next g
   } // si n ou m
   flocusIn.close();
} // fin lecture


void tabFtotabM(vector<vector<unsigned long int> > * tabF){
// on passe des tables construites par (CT_tests)crunchLocTable() aux tabM et tabCode du code preexistant
using namespace NS_F_est;
   nb_sam=fichier_genepop->pops.size();
   global_geno_nbr=(*tabF)[0].size();
   if(nb_sam == 0 || global_geno_nbr == 0){
   }else{
      // dimensionnement des tab
//cout<<"dim a "<<nb_sam<<" "<<global_geno_nbr;getchar();
      tabM=new size_t*[nb_sam];
      for(size_t pop = 0; pop < nb_sam; pop++) {
   		tabM[pop]=new size_t[global_geno_nbr];
      }


      tabCode=new size_t*[global_geno_nbr];
      for(size_t g = 0; g < global_geno_nbr; g++) tabCode[g]=new size_t[2];
      // lecture e partir ligne 3 du fichier : effectif genotypiques par pop
      for(size_t pop = 0; pop < nb_sam; pop++){
         for(size_t g = 0; g < global_geno_nbr; g++){
            tabM[pop][g]=(*tabF)[pop][g];
         } // next g
      } // next pop
      for(size_t g = 0; g < global_geno_nbr; g++){
         tabCode[g][0]=(*tabF)[nb_sam][g]; // lit numero 1e allele du genotype
      } // next j
      for(size_t g = 0; g < global_geno_nbr; g++){
// pour option5.2, 5.3 meme haploide (*tabF)[nb_sam+1] existe (cf crunchLocTable)
         tabCode[g][1]=(*tabF)[nb_sam+1][g]; // lit numero 2e allele du genotype
      } // next g
   } // si n ou m
} // fin lecture


void FisParPop(bool identitybool,size_t iLoc,ofstream &fic_out,
               std::vector<double>& MSgTotHez,
               vector<int>& NLocHez,
               vector<int>& NLoc,
               vector<int>& NtotLoc,
               vector<double>& SSiTotLoc,
               vector<double>& SSgTotLoc) {
// also diversities, even haploid
using namespace NS_F_est;
using namespace NS_tailles;
//int jj,i;
vector<unsigned long int> scgsp;
vector<vector<double> >sfreqs;
std::vector<double> tailleMoy;
// point d'entree pour les deux variables suivantes:
  identity=identitybool;
F_est_locIt=iLoc;
nb_locus=fichier_genepop->loci.size();
double SSiLocTtesPops=0,SSgLocTtesPops=0;
double NcLoc=0; //04/2007 for weighted MSi over pops per locus
size_t Ntot=0;
if (global_geno_nbr>0) {// valued by tabFtotabM previously called; sinon plantage dans calc_sfreqs_Nc
   for (global_pop_it=0;global_pop_it<nb_sam;global_pop_it++) {
		 if (identity) {
		 	maxallname=size_t(fichier_genepop->loci[F_est_locIt]->alleleMax);
			sfreqs.resize(3);
			for(size_t ii=0;ii<3;ii++) {sfreqs[ii].resize(maxallname+1);}
//		 sfreqs=dmatrix(0,deuxsamp+1,0,maxall);
			 for (size_t i=0;i<3;i++)
					for (size_t ii=0;ii<=maxallname;ii++) sfreqs[i][ii]=0.0;
		  }
		 else {
            tailleMoy.resize(1);
            fill(tailleMoy.begin(),tailleMoy.end(),0);
            tailleMoyTot=0;
         }
//		 else taille_ipp=dvector(0,deuxsamp-1);
         calc_sfreqs_Nc(1,tailleMoy,scgsp,sfreqs,Ntot); // pour chaq pop (+tailleMoyTot)
		 if (Ntot> 0)  {
			calculSSetMS(tailleMoy,scgsp,sfreqs,Ntot);  // par locus (et par pop car deuxsamp=1)
            if (fichier_genepop->coding[iLoc]>3) {
               fic_out<<left<<setw(11)<<fichier_genepop->pops[global_pop_it]->popName().substr(0,10);
               fic_out<<setw(11)<<MSg<<setw(11)<<(MSi+MSg)/2;
               ios_base::fmtflags old=fic_out.setf(ios_base::fixed,ios_base::floatfield);
               fic_out.precision(4);
               if (MSg+MSi>0) fic_out<<internal<<setw(7)<<(MSi-MSg)/(MSg+MSi)<<endl;
               else fic_out<<"    -  \n";
               fic_out.setf(old,ios_base::floatfield);
               if (identitybool) fic_out.precision(4); else fic_out.precision(5);
             // {unweighted average diversities}
               NLocHez[global_pop_it]++;  //diploid specific
               MSgTotHez[global_pop_it]+=MSg; //diploid specific
               if (Ntot>1) {
// per locus sums over pops
                    SSiLocTtesPops+=Nc*MSi; //locus specific
                    SSgLocTtesPops+=Nc*MSg; //diploid specific
                    NcLoc+=Nc; //locus specific
// sums over loci
                    SSgTotLoc[global_pop_it]+=Ntot*MSg; //diploid specific
                    NLoc[global_pop_it]++; //diploid specific (only serves for unweighted 1-Qintra)
                    if (estimDiploidBool) {
                        SSiTotLoc[global_pop_it]+= Ntot*MSi; //ploidy specific
                        NtotLoc[global_pop_it]+=Ntot; //ploidy specific
                    }
               }
            } else { // haploid
               fic_out<<left<<setw(11)<<fichier_genepop->pops[global_pop_it]->popName().substr(0,10);
               fic_out<<setw(11)<<MSi<<endl;
               if (Ntot>1) {
// per locus sums over pops
                    SSiLocTtesPops+=Nc*MSi; //locus specific
                    NcLoc+=Nc; //locus specific
// sums over loci
//                    NLoc[global_pop_it]++; //diploid specific
                    if (!estimDiploidBool) {
                        SSiTotLoc[global_pop_it]+= Ntot*MSi; //ploidy specific
                        NtotLoc[global_pop_it]+=Ntot; //ploidy specific
                    }
               }
            }
		} else { //Ntot=0
            if (fichier_genepop->coding[iLoc]<4) fic_out<<left<<setw(11)<<fichier_genepop->pops[global_pop_it]->popName().substr(0,10)<<"    -\n";
            else fic_out<<left<<setw(11)<<fichier_genepop->pops[global_pop_it]->popName().substr(0,10)<<"    -          -          -"<<endl;
		}
   }  /* pop*/
   if (fichier_genepop->coding[iLoc]<4) { //haploid
       fic_out<<"-------------------\n";
       fic_out<<right<<setw(11)<<fichier_genepop->loci[F_est_locIt]->locusName.substr(0,10)<<", All samples:\n";
       if (NcLoc>0) fic_out<<setw(11)<<" "<<left<<setw(11)<<SSiLocTtesPops/NcLoc<<endl; //04/2007 weighted MSi's
   } else { //diploid
       fic_out<<"-------------------------------------------\n";
       fic_out<<right<<setw(11)<<fichier_genepop->loci[F_est_locIt]->locusName.substr(0,10)<<", All samples:\n";
       fic_out<<right<<setw(22)<<" ";
       if (NcLoc>0) fic_out<<left<<setw(11)<<(SSiLocTtesPops+SSgLocTtesPops)/(2*NcLoc);
       else fic_out<<setw(11)<<" "; //04/2007 weighted MSi's
       if (SSiLocTtesPops+SSgLocTtesPops>0) {
          ios_base::fmtflags old=fic_out.setf(ios_base::fixed,ios_base::floatfield);
          fic_out.precision(4);
          fic_out<<internal<<setw(7)<<(SSiLocTtesPops-SSgLocTtesPops)/(SSiLocTtesPops+SSgLocTtesPops)<<endl;
          fic_out.setf(old,ios_base::floatfield);
          if (identitybool) fic_out.precision(4); else fic_out.precision(5);
       }
       else fic_out<<"    -\n";
   }
} else {
   if (fichier_genepop->coding[iLoc]<4) fic_out<<"No genotypes.\n";
   else fic_out<<"No complete diploid genotypes.\n";
}
//blob("sortie FisParPop()");
} //end FIsParPop

void hierFstat(bool identitybool,int Indic,size_t iLoc,ofstream &fic_out, vector<vector<double> > * FFF) {
//cerr<<"debut hierFstat";
using namespace NS_F_est;
using namespace NS_FFF_slmt;
using namespace NS_tailles;
//int jj,i;
vector<unsigned long int> scgsp;
std::vector<double> tailleMoy;
vector<vector<double> >sfreqs;
// point d'entree pour les deux variables suivantes:
identity=identitybool;
F_est_locIt=iLoc;
double deno;
size_t Ntot=0;
//		 else taille_ipp=dvector(0,deuxsamp-1);
    if (Indic==2) { //sur paires ////////////////////////////////////////////////////////
        fic_out<<setprecision(4)<<internal;
        if (identity) {
            maxallname=size_t(fichier_genepop->loci[F_est_locIt]->alleleMax);
            sfreqs.resize(4);
            for(size_t ii=0;ii<4;ii++) {sfreqs[ii].resize(maxallname+1);}
            //		 sfreqs=dmatrix(0,deuxsamp+1,0,maxall);
        } else {
            tailleMoy.resize(2);
            tailleMoyTot=0;
        }
        for (global_pop_it=1;global_pop_it<nb_sam;global_pop_it++) {
            fic_out<<setw(6)<<left<<global_pop_it+1<<internal;
        	for (global_pop2_it=0;global_pop2_it<global_pop_it;global_pop2_it++) {
                if (identity) {
                     for (size_t i=0;i<4;i++)
                    		for (size_t ii=0;ii<=maxallname;ii++) sfreqs[i][ii]=0.0;
                } else {
                    fill(tailleMoy.begin(),tailleMoy.end(),0);
                    tailleMoyTot=0;
                }
//computation of locus-specific terms
                calc_sfreqs_Nc(Indic,tailleMoy,scgsp,sfreqs,Ntot); // Indic=2 => pour la pair (global_pop_it,global_pop2_it)
                calculSSetMS(tailleMoy,scgsp,sfreqs,Ntot);
                deno=MSp+(Nc-1)*MSi+Nc*MSg; // note MSg=0 en haploid; formula still correct in this case
                if (deno>0 && nonvides > 1.5 ) {
                   fic_out<<setw(7)<<(MSp-MSi)/deno<<" ";
                } else fic_out<< "   -    ";
//multilocus terms
                if (nonvides>1.5) {
                   if ((estimDiploidBool && fichier_genepop->coding[iLoc]>3) || (!estimDiploidBool && fichier_genepop->coding[iLoc]<4)) {
// sum only on loci with the appropriate ploidy
                       (*FFF)[global_pop_it][global_pop2_it]+=(MSp-MSi);
                       (*FFF)[global_pop2_it][global_pop_it]+=deno;
                   }
                }
            }
            fic_out<<endl;
        }
    } else { //  pas sur paires //////////////////////////////////////////////////////////////////////////
        if (identity) {
            maxallname=size_t(fichier_genepop->loci[F_est_locIt]->alleleMax);
            sfreqs.resize(nb_sam+2);
            for(size_t ii=0;ii<nb_sam+2;ii++) {sfreqs[ii].resize(maxallname+1);}
            //		 sfreqs=dmatrix(0,deuxsamp+1,0,maxall);
             for (size_t i=0;i<nb_sam+2;i++)
            		for (size_t ii=0;ii<=maxallname;ii++) sfreqs[i][ii]=0.0;
        } else {
            tailleMoy.resize(nb_sam);
            fill(tailleMoy.begin(),tailleMoy.end(),0);
            tailleMoyTot=0;
        }
        calc_sfreqs_Nc(Indic,tailleMoy,scgsp,sfreqs,Ntot); // Indic=0 => sur ttes les pops
        if (Ntot> 0)  {
//computation of locus-specific terms
            calculSSetMS(tailleMoy,scgsp,sfreqs,Ntot);  // par locus sur ttes les pops ///// \u2259 = "estimates" not supported by many text editors
            if (fichier_genepop->coding[iLoc]>3) { //diploid
                if (MSi+MSg>0) {
                    if (identitybool) fic_out<<"\nFis^= "; else fic_out<<"\nRhois^= ";
                    fic_out<<(MSi-MSg)/(MSi+MSg);
                    (*FFF)[iLoc].push_back((MSi-MSg)/(MSi+MSg));
                } else (*FFF)[iLoc].push_back(-666);
                deno=MSp+(Nc-1)*MSi+Nc*MSg;
                if (deno>0 && nonvides > 1.5 ) {
                   if (identitybool) fic_out<<"\nFst^= "; else fic_out<<"\nRhost^= ";
                   fic_out<<(MSp-MSi)/deno;
                   (*FFF)[iLoc].push_back((MSp-MSi)/deno);
                   if (identitybool) fic_out<<"\nFit^= "; else fic_out<<"\nRhoit^= ";
                   fic_out<<1.-(2.*Nc*MSg)/deno;
                   (*FFF)[iLoc].push_back(1.-(2.*Nc*MSg)/deno);
                } else {(*FFF)[iLoc].push_back(-666);(*FFF)[iLoc].push_back(-666);}
                if (identitybool) fic_out<<"\n1-Qintra^= "; else fic_out<<"\nSSDintra^= ";
                fic_out<<MSg<<" (gene diversity intra-individuals)";
                if (identitybool) fic_out<<"\n1-Qinter^= "; else fic_out<<"\nSSDinter^= ";
                fic_out<<(MSi+MSg)/2<<" (gene diversity inter-individuals intra-pop)";
        //cout << Nc<<" "<<MSi<<" "<<MSg<<" "<<SSiTot<<" "<<SSgTot;getchar();
            } else { //haploid
                (*FFF)[iLoc].push_back(-666); // Fis
                deno=MSp+(Nc-1)*MSi;
                if (deno>0 && nonvides > 1.5 ) {
                   if (identitybool) fic_out<<"\nFst^= "; else fic_out<<"\nRhost^= ";
                   fic_out<<(MSp-MSi)/deno;
                   (*FFF)[iLoc].push_back((MSp-MSi)/deno); //Fst
                } else (*FFF)[iLoc].push_back(-666);  //Fst
                (*FFF)[iLoc].push_back(-666); //Fit
                if (identitybool) fic_out<<"\n1-Qinter^= "; else fic_out<<"\nSSDinter^= ";
                fic_out<<MSi<<" (gene diversity inter-individuals intra-pop)";
            }
//multilocus terms
// sum only on loci with the appropriate ploidy
            if ((estimDiploidBool && fichier_genepop->coding[iLoc]>3) || (!estimDiploidBool && fichier_genepop->coding[iLoc]<4)) {
                SSiTot+=Nc*MSi;
                SSgTot+=Nc*MSg;
                if (nonvides>0 && Ntot-nonvides>0) {
                    MSi2P+=MSi;
                    MSp2P+=MSp;
                    MSg2P+=MSg;
                    if (nonvides>1) { // ne sert que pour Fst
                       MSi2Pw+=(Nc-1)*MSi;
                       MSg2Pw+=Nc*MSg;
                    }
                }
            }
        } else {
               (*FFF)[iLoc].push_back(-666);(*FFF)[iLoc].push_back(-666);(*FFF)[iLoc].push_back(-666);
        }
        fic_out<<endl<<endl;
    } // end test sur Indic ////////////////////////////////////////////////////
//cerr<<"fin hierFstat";
}


void calc_sfreqs_Nc(int Indic,vector<double>& tailleMoy,
                    vector<unsigned long int>& scgsp,
                    vector<vector<double> >& sfreqs,
                    size_t& Ntot) { //sfreqs ou tailles selon identity
// ne pas appeller quand global_geno_nbr=0 !!!!!
// indic=1 => Fis par pop
// indic=2 => F par paires
// indic=0 => F total table
using namespace NS_F_est;
using namespace NS_tailles;
    indic12all=Indic; //copie dans variable membre de NS_F_est
    double moy;
    size_t *tabMind=NULL;
    size_t tabMindj;
    nonvides=0;
    if (indic12all==2) {
       indices[0]=global_pop_it;
       indices[1]=global_pop2_it;
       deuxsamp=2;
    } else if (indic12all==1) {
       indices[0]=global_pop_it;
       deuxsamp=1;
    } else {deuxsamp=nb_sam;}
//cerr<<global_pop_it<<" "<<global_pop2_it<<" "<<deuxsamp<<endl;

// ATTENTION les tabcodes sont les noms d'alleles et donc l'iterateur est borne par le nom max,
//pas par le nombre max d'alleles (et sferqs dimensionne en consequence e maxallname+1)

    Ntot=0;	
    double Ntot2=0;
    scgsp.resize(deuxsamp);
	for (size_t i=0;i<deuxsamp;i++) {
		 scgsp[i]=0; //effectif sample [i]...
//evid le code suivant plante si tabM n'a pas ete alloue quand global_geno_nbr=0
		 if (global_geno_nbr>0) {
		 	if (indic12all==1 ||indic12all==2) tabMind=tabM[indices[i]]; // pointe sur le vecteur des effectifs genos pour la pop donnee
		 	else tabMind=tabM[i];
		 } // else the following loop does nothing
		 for (size_t j=0;j<global_geno_nbr;j++) {
//			  fscanf(f_pairw,"%d",&geno[i][j]);
              		  tabMindj=(tabMind)[j];
			  scgsp[i]+=tabMindj;
			  /*Fst:*/
			  if (identity) {
//sfreqs a ete reinit e 0 juste avant chaque appel de calc_sfreqs_Nc() dans lecturePaires()
//cout<<indices[i]<<" "<<j<<" "<<tabM[indices[i]][j];getchar();
//cout<<"ident=true! !!";getchar();
            	  		sfreqs[i][tabCode[j][0]]+=tabMindj; //2 lignes de rajouts d'effectifs alleliq pour chaque geno dans chaque pop
            	  		sfreqs[i][tabCode[j][1]]+=tabMindj; // en haploide tabCode[j][1]=0 et seulement dans ce cas
            	  		sfreqs[deuxsamp][tabCode[j][0]]+=tabMindj; // effectifs alleliq marginal sur toutes les pops
            	  		sfreqs[deuxsamp][tabCode[j][1]]+=tabMindj; //(bis) en haploide tabCode[j][1]=0 et seulement dans ce cas
            	  		if (tabCode[j][0]!=tabCode[j][1]) { // junk mais potentiellement vrai en haploid !
            		  		sfreqs[deuxsamp+1][tabCode[j][0]]+=tabMindj; // nb hez par allele marginal sur ttes les pops
            		  		sfreqs[deuxsamp+1][tabCode[j][1]]+=tabMindj;
				}
			  } else {
//cout<<"ident=false OK";getchar();
        		  	if (tabCode[0][1]>0) //diploidy
                     		moy=(tailleOfType(tabCode[j][0])+tailleOfType(tabCode[j][1]))/2.;
                  		else moy=tailleOfType(tabCode[j][0]); //tailleOfType(tabCode[j][0]=0)=0 en principe
                  		tailleMoy[i]+=moy*tabMindj;
              		  }/* identity */
//if (deuxsamp!=2) cerr<<"("<<j<<" "<<tabCode[j][0]<<" "<<tabCode[j][1]<<" "<<(tabMind)[j]<<")";
		 } /*j*/
		 if (scgsp[i]>0) nonvides++;
		 Ntot+=scgsp[i]; //n1+n2 //FER vire (+)1.0*(scgsp[i]) bizarre int-> double -> int. Idem ci dessous
		 Ntot2+=scgsp[i]*scgsp[i]; //
	} /*i*/
//cerr<<"c";
	if (Ntot==0) Nc=Ntot;  //variable nonvides remise 09/006
	else {
// probleme pour les paires d'individus
        if (nonvides==1) {
            if(Indic==1) Nc=Ntot; //Nc sert pour les Fis multilocus; il adonc un sens si nonvides=1
            else Nc=0;// mais il ne faut surtout pas que Nc=1 pour les paires d'indiv (indic=2) incompletes car:
// =1 pour les paires d'individus completes:
        } else Nc= (Ntot-Ntot2/Ntot)/(nonvides-1);
		  /*Fst:*/
        if (identity) {
            for (size_t i=0;i<deuxsamp;i++) {
            	if (scgsp[i]>0) {
            		if (tabCode[0][1]>0) {//diploidy
                       for (size_t ii=1;ii<=maxallname;ii++) sfreqs[i][ii]/=(2.*scgsp[i]);
                    } else for (size_t ii=1;ii<=maxallname;ii++) sfreqs[i][ii]/=scgsp[i];
            	}
            }
            for (size_t ii=1;ii<=maxallname;ii++)  {
//cerr<<"d";
        		  if (tabCode[0][1]>0) {//diploidy
            	     sfreqs[deuxsamp][ii]/=(2.*Ntot);
                  } else sfreqs[deuxsamp][ii]/=Ntot;
            	  sfreqs[deuxsamp+1][ii]/=Ntot;
            /*				  printf("%4.2f %d/%d\n",sfreqs[deuxsamp][ii],locus,nblocus);*/
     	    } /*for ii*/
        } else {
            for (unsigned int i=0;i<deuxsamp;i++) {
           		tailleMoyTot+=tailleMoy[i];
            	if (scgsp[i]>0) {
            		tailleMoy[i]/=(scgsp[i]);
                }
          	}
       		tailleMoyTot/=Ntot;
        }
	} /*if Ntot else*/
}


void calculSSetMS(std::vector<double>& tailleMoy,
                  vector<unsigned long int>& scgsp,
                  vector<vector<double> >& sfreqs,
                  size_t& Ntot) {    /* calcul SSp et MSp; appele une fois par locus */
// here from hierFstat (options 6.1--6.4) or from isolde_etc()
/* Noter les sommes explicites sur les alleles [ii=1;ii<=maxallname;ii++]
si on ne prend qu'un terme de la somme on arrive finalement e un Fst par l'allele donne.
Expressions semblables e celles de Weir, GDA2, p.185.
*/
// les indices[] ont ete valu'es dans calc_sfreqs_Nc
using namespace NS_F_est;
using namespace NS_tailles;
//int i,j,ii;
    double moy;
    size_t *tabMind;
	/* calcul SSp et MSp */
    if (global_geno_nbr==0) {
	Ntot=0;MSg=0.;MSp=0.;MSi=0;SSg=0.;SSi=0.;SSp=0.;
	return;
    } //05/2011 because otherwise the following code access tabcode which has not been correctly initialised
    SSp=0.;
    if (identity) {/*Fst:*/
         for (size_t i=0;i<deuxsamp;i++) {
        	  for (size_t ii=1;ii<=maxallname;ii++) {
        			  SSp=SSp + pow(sfreqs[i][ii]-sfreqs[deuxsamp][ii],2)*scgsp[i];
//cout<<sfreqs[i][ii]<<" "<<sfreqs[deuxsamp][ii];getchar();
                    }
        	 }
	      if (tabCode[0][1]>0) // not haploid  04/2007
             SSp=2*SSp;
    } else {/*Rho:*/
        for (size_t i=0;i<deuxsamp;i++)
        		SSp=SSp + pow(tailleMoy[i]-tailleMoyTot,2)*scgsp[i];
	      if (tabCode[0][1]>0) // not haploid 04/2007
             SSp=4*SSp; //for the 4, see below for SSi and SSg
          else SSp=2*SSp;
    }
    if (nonvides<2) MSp=0.; else MSp=SSp/(nonvides-1);  /*une global_pop_nbr vide =>MSp=0*/

	 /* calcul SSi et MSi */
	 SSi=0.;
	 if (!identity) {	 /*rho:*/
		  for (size_t i=0;i<deuxsamp;i++) {
                if (indic12all==1 ||indic12all==2) tabMind=tabM[indices[i]]; // pointe sur le vecteur des effectifs genos pour la pop donnee
                else tabMind=tabM[i];
				for (size_t j=0;j<global_geno_nbr;j++) {
        		    if (tabCode[0][1]>0) {//diploidy
                       moy=(tailleOfType(tabCode[j][0])+tailleOfType(tabCode[j][1]))/2.;
					   moy= pow(moy-tailleMoy[i],2);
				       SSi+=4.0*(tabMind)[j]*moy; // to fit with SSg below
                    } else { //haploidy
                       moy=tailleOfType(tabCode[j][0]); //tailleOfType(tabCode[j][1]=0)=0 en principe
					   moy= pow(moy-tailleMoy[i],2);
				       SSi+=2.0*(tabMind)[j]*moy;
                    }
				} /* j */
		  } /* i */
	 }
	 else { /*Fst:*/
		  for (size_t i=0;i<deuxsamp;i++) {
				for (size_t ii=1;ii<=maxallname;ii++) {
//sommes sur les individus. la sommation est deje operee dan les formules de WeirC84,
// d'oe leur terme supplementaire
        		  if (tabCode[0][1]>0) //diploidy
					 SSi=SSi+2.0*scgsp[i]*(sfreqs[i][ii]-pow(sfreqs[i][ii],2)
											  -sfreqs[deuxsamp+1][ii]/4);
                  else SSi=SSi+scgsp[i]*(sfreqs[i][ii]-pow(sfreqs[i][ii],2));  // for haploid indiv data sfreq[][]= 0 or 1 hence SSi=0
             }
		  } /* i */
	 } /* condition */

	 if ((Ntot-nonvides)<0.5) MSi=0.; else MSi=SSi/(Ntot-nonvides);
/* test Vrai si un ind maxi /pop =>MSi=0 mais aussi si une des deux pop est vide
  mais on peut avoir MSi>0 et Nc<=1; ce MSi n'intervient pas dans le calcul
du ST global dans mat_ST grace au test Nc>1.0001 */

	 /* calcul SSg et MSg */
	 SSg=0.;
	 if (tabCode[0][1]>0) {// i.e. not haploid; si haploid SSg reste nul
    	 if (!identity) {	 /*rho:*/
    		  for (size_t i=0;i<deuxsamp;i++) {
                 if (indic12all==1 ||indic12all==2) tabMind=tabM[indices[i]]; // pointe sur le vecteur des effectifs genos pour la pop donnee
                 else tabMind=tabM[i];
    			 for (size_t j=0;j<global_geno_nbr;j++) {
                    moy= (tailleOfType(tabCode[j][0])+tailleOfType(tabCode[j][1]))/2.;
                    SSg= SSg+ (tabMind)[j]*pow(tailleOfType(tabCode[j][0])-moy,2);
                    SSg= SSg+ (tabMind)[j]*pow(tailleOfType(tabCode[j][1])-moy,2);
    		     } /* j */
    		  } /*i */
    	 SSg*=2; //for sizes k and k+1: adds one for each hez => identical MSs whether (identity) or not
    	 }
    	 else 	 /*Fst:*/
    	 { for (size_t i=0;i<deuxsamp;i++) {
    				for (size_t ii=1;ii<=maxallname;ii++) {
//if (deuxsamp!=2) cerr<<"("<<scgsp[i]<<" "<<sfreqs[deuxsamp+1][ii]<<")";
                        SSg=SSg+scgsp[i]*sfreqs[deuxsamp+1][ii]/2;
    					  /*pas economique, mais clair*/
                     }
    		 } /* i */
    	 } /* identity */
     }
//if (deuxsamp!=2) {cerr<<"SSg "<<SSg;getchar();}
	 if (Ntot==0) MSg=0.; else MSg=SSg/Ntot; // test change de Nc e Ntot 21/10/2006
/*if (MSg*(1-MSg)*(0.5-MSg)>0.001) {
cout<<endl<< global_pop_it<<" "<<global_pop2_it<<endl;
for (size_t j=0;j<global_geno_nbr;j++) {cout<<" "<<tabCode[j][0]<<" "<<tabCode[j][1]<<endl;}
cout<<endl;
for (int i=0;i<2;i++) {
for (size_t j=0;j<global_geno_nbr;j++) {cout<<" "<<tabM[indices[i]][j];}
cout<<endl;
for (size_t j=0;j<maxallname;j++) {cout<<" "<<sfreqs[i][j];}
cout<<endl;
}
for (size_t j=0;j<maxallname;j++) {cout<<" "<<sfreqs[2][j];}
cout<<endl;
getchar();
}*/
  /*qui est aussi=0 si
			on a aucun vrai heteroz*/
//cout<<endl<<MSi<<": "<<MSp;getchar();
}  /*calculSSetMS*/


void lecturePaires(void) {
//cout<<"lecturePaires";
vector<unsigned long int> scgsp;
vector<double>tailleMoy;
vector<vector<double> >sfreqs;
using namespace NS_F_est;
using namespace NS_tailles;
//int jj,i;
size_t Ntot=0;

		 if (identity) {
//cerr<<"identity=true au debut lecturePaires()"<<identity<<endl;getchar();
		 	maxallname=size_t(fichier_genepop->loci[F_est_locIt]->alleleMax);
		 	{size_t li=nb_sam*nb_sam*maxallname;
		 	if(li>100000) {
                noR_cout<<"\nMay be slow...";
            }
            else {

                noR_cout<<"\n                         ";

            }
			}
			sfreqs.resize(4);
			for(size_t ii=0;ii<4;ii++) {sfreqs[ii].resize(maxallname+1);}
         } else {
//cerr<<"identity=false au debut lecturePaires()"<<identity<<endl;getchar();
            tailleMoy.resize(2);
         }

for (global_pop_it=1;global_pop_it<nb_sam;global_pop_it++) {
	for (global_pop2_it=0;global_pop2_it<global_pop_it;global_pop2_it++) {
		 if (identity) {
			 for (size_t i=0;i<4;i++)
					for (size_t ii=0;ii<=maxallname;ii++) sfreqs[i][ii]=0.0;
         } else {
            tailleMoyTot=0;
         }
		 calc_sfreqs_Nc(2,tailleMoy,scgsp,sfreqs,Ntot); // pour chaque paire
         MStableit=&MStable[human_loc_it-1][global_pop_it-1][global_pop2_it];
		 if (Ntot> 0)  {
			 calculSSetMS(tailleMoy,scgsp,sfreqs,Ntot);  //tjrs par locus
			MStableit->bbordel=human_loc_it;
			MStableit->mmsp=MSp;
			MStableit->mmsi=MSi;
			MStableit->mmsg=MSg;
			MStableit->nnc=Nc;
		} else {
			MStableit->bbordel=human_loc_it;
			MStableit->mmsp=0;
			MStableit->mmsi=0;
			MStableit->mmsg=0;
			MStableit->nnc=0;
		}
//cout<<global_pop_it<<" "<<global_pop2_it;
	}/* de jj*/
//cout<<global_pop_it<<" ";
 }  /* de iii*/
/*fclose(f_res); */
}

void delete_tabM_tabCode() {
// pour la fin de 5.2     // sinon aussi dans main2x2
using namespace NS_F_est;
//cout<<"delete a"<<nb_sam<<" "<<global_geno_nbr;getchar();
   if(nb_sam == 0 || global_geno_nbr == 0){
   }else{
        for(size_t pop = 0; pop < nb_sam; pop++) delete[] tabM[pop];
        delete[] tabM;
        for(size_t g = 0; g < global_geno_nbr; g++) delete[] tabCode[g];
        delete[] tabCode;
    }
}





int main2x2(vector<bool>& ploidBool) { //called only in main2x2 hence in isolde_etc -> ... -> pairwMS -> main2x2
//cout<<"main2x2";getchar();
using namespace NS_F_est;
using namespace NS_GP;
   if (!_perf) effacer_ecran();
#ifdef COMPATIBILITYRCPP
   Rcpp::Rcerr<<"Computing pairwise Fst's or analogous correlations:"<<endl;
   Progress progbar(nb_locus, true);
#else   
   _gotoxy(0,5);
   noR_cout<<" Computing pairwise Fst's or analogous correlations";
   /* RHO */
   //if (!identity) {cout<<"!identity pas implemente";getchar();} //lectureTaille();
   _gotoxy(19,6);
   
   noR_cout<<"Already            loci considered, out of "<<nb_locus;
#endif


for (F_est_locIt=0;F_est_locIt<nb_locus;F_est_locIt++) {
	 human_loc_it=F_est_locIt+1; // :-)
		  /*FST:*/
	global_geno_nbr=0; // will read a new value in LOCUSn:
	if(ploidBool[F_est_locIt]) {
         lecture_floc(); //reads LOCUSn files
    }
	if (global_geno_nbr>1) { // s'il a trouve des genos dans LOCUSn
	  deuxsamp=2; //global variable: info for lecturePaires-> (calc_sfreqs_Nc and calculSSetMS)
		lecturePaires();// calcule les MS de l'analyse de variance
	} else {
		for (size_t k=1;k<nb_sam;k++) //FR corrected 10/2010
			for (size_t kpop=0;kpop<k;kpop++)
#ifdef RHO_FILE
			  fprintf(f_rho,"%d    0        0         0         0\n",human_loc_it);
#else
			{
                MStableit=&MStable[human_loc_it-1][k-1][kpop];
			    MStableit->bbordel=human_loc_it;
			    MStableit->mmsp=0;
			    MStableit->mmsi=0;
			    MStableit->mmsg=0;
			    MStableit->nnc=0;
			}
#endif

	 } /*if global_geno_nbr>1*/
#ifdef COMPATIBILITYRCPP
	progbar.increment();
#else   
	_gotoxy(31,6);
	noR_cout<<human_loc_it<<"  ";
#endif

// libere les allocations effectuees dans lecture_floc
      if(nb_sam == 0 || global_geno_nbr == 0){
// *meme test que dans lecture_floc* sauf si difference entre nb_sam et nb_sam_lu....
      }else{
          for(size_t pop = 0; pop < nb_sam; pop++) delete[] tabM[pop];
          delete[] tabM;
          for(size_t g = 0; g < global_geno_nbr; g++) delete[] tabCode[g];
          delete[] tabCode;
      }
  } /*de l'iteration sur les locus*/


#ifdef RHO_FILE
  fclose(f_rho);
#endif
return 0;
} /* main2x2 */

int pairwMS(vector<bool>& ploidBool){ // called only in isolde_etc
using namespace NS_F_est;
    if(nb_sam > 1){
      if(nb_locus != 0){
      main2x2(ploidBool);  //contient boucle sur loci
      } // si locus
    } // si pop
return 0;
} // fin main

void readGGFile(const char nom_fich_mig[]) {
//cout<<"begin readGGFile ";
//cout<<nom_fich_mig<<endl;
// reecrit 08/006  			fscanf(ipt, "%lg", &data[s2][s1]); ne lit pas bien sous dev-c++
using namespace datamatrix;
using namespace NS_F_est;
//char atuer[50];

/*   \  gene
		 \
	 geog \
*/
//cout<<nb_sam_migf<<" "<<nb_pair_sam_migf;getchar();
char ch;
long int Nmiss=0;
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
            // Rcpp::R
        #else
            cerr<<"\n This file exists but is empty. I must exit.";
            ipt.clear();
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(1, "This file exists but is empty.");
    } else while (ipt.get()!= '\n') continue; //skipln
    ipt>>nb_sam_migf;
//cout<<"\anb_sam_migf"<<nb_sam_migf;getchar();
    while (ipt.get()!= '\n') continue; //skipln
    while (ipt.get()!= '\n') continue; //skipln
    data.resize(nb_sam_migf);
    for (vector<vector<long double> >::iterator ii=data.begin();ii!=data.end();ii++) (*ii).resize(nb_sam_migf);
//cout<<nb_sam_migf;getchar();
    if (_first_of_repl) {  /// then we print a *warning* for NaN's
        /// context is: isolde_etc creates a migFile for each bootstrap replicate then calls this function
    	for (size_t s1=1; s1<nb_sam_migf; s1++) {
    		for (size_t s2=0; s2<s1; s2++) {
//cout<<"\r"<<s1<<" "<<s2;
// C++ does not seem to read a quiet_NaN even if it has written a quiet_NaN. So we write the string "NaN" rather than a quiet_NaN and read a string
                ipt>>prefiltre;
                if(prefiltre=="NaN") {
     				Nmiss++;

                        if (Nmiss<4) noR_cout<<"\a\n No genetic information for pair "<<s1+1<<" and "<<s2+1;
     				    else if (Nmiss==4) noR_cout<<"\a\n (more pairs without genetic information)";


                    data[s2][s1]=numeric_limits<long double>::quiet_NaN();
                } else {
                    filtre.str(prefiltre);
                    filtre>>data[s2][s1];
                    filtre.str("");
                    filtre.clear();
                }
//cout<<"\n"<<s1<<" "<<s2<<" "<<data[s2][s1];
     		}
//getchar();
    	}
//cout<<"\r Finished reading genetic distance matrix";
    	if (Nmiss>0) {

                if (Nmiss>3) noR_cout<<"\a\n  "<<Nmiss<<" pairs without genetic information"<<endl;
        		noR_cout<<"\n\n For pairs of individuals, this will typically occur when such pairs"<<endl;
        		noR_cout<<"have no genotyped locus in common that are polymorphic in the population."<<endl;
        		noR_cout<<"\n The analysis can nevertheless proceed."<<endl;
                if (pauseGP) { noR_cout<<"\n\n(Return) to continue"<<endl; getchar();}

    	}
    } else { //bootstrap replicates : pas besoin d'annoncer les NaN's
    	for (size_t s1=1; s1<nb_sam_migf; s1++)
    		for (size_t s2=0; s2<s1; s2++) {
                ipt>>prefiltre;
                if(prefiltre=="NaN") {
                    data[s2][s1]=numeric_limits<long double>::quiet_NaN();
                } else {
                    filtre.str(prefiltre);
                    filtre>>data[s2][s1];
                    filtre.str("");
                    filtre.clear();
                }
     		}
    }
// puis lecture demi matrice geo. Only needed if first replicate, otherwise the transformed geo matrix are kept in datamatrix::data
    if (_first_of_repl) {
        while (ipt.get()!= '\n') continue; //skipln
        while (ipt.get()!= '\n') continue; //skipln
        for (size_t s1=1; s1<nb_sam_migf; s1++)
            for (size_t s2=0; s2<s1; s2++) {
                ipt>>data[s1][s2];
                if (ipt.fail()) {
                    #ifdef COMPATIBILITYRCPP
                        // Rcpp::R
                    #else
                        cerr<<"\a\n Incomplete distance matrix! Check input file.\n...the program must terminate...\n";
                        if (cinGetOnError) cin.get();
                        ipt.clear();
                    #endif
                    genepop_exit(1, "Incomplete distance matrix! Check input file.");
                }
            }
    }
    ipt.close();
//cout<<"\aend readGGFile";getchar();
} /*end, readGGFile */


void writedat(vector<vector<long double> >m,const char nom_fich_mig[]);
void writedat(vector<vector<long double> >m,const char nom_fich_mig[])	{
using namespace datamatrix;
using namespace NS_F_est;
char char_gra[]=".GRA";
string graphname;
//  char dumname[70],graphname[70];
double truncated;
//  string graphname;
//  strcpy(dumname,nom_fich_mig);
//  outname=strcat(strtok(dumname,"."),char_iso);
    outname=gp_file+char_iso;
    //  strcpy(graphname,strcat(strtok(dumname,"."),char_gra));
    graphname=gp_file+char_gra;
    ofstream opt(outname.c_str()); // iso file
    ofstream gra(graphname.c_str()); // gra file
    opt<<nb_sam_migf<<" populations ("<<nom_fich_mig<<")\n";
    opt<<"\ngenetic estimates ("<<statname<<"):\n";
    for (size_t s1=1; s1<nb_sam_migf; s1++) {
        for (size_t s2=0; s2<s1; s2++)
            if (isnan(m[s2][s1])) opt<<"     -    ";
            else {
//"C++ never truncates data"... donc il va systematiquement sortir 6 chiffres signif
//meme s'il commence e 0.00... et ea va depasser la width()
                truncated=double(floor(m[s2][s1]*1000000.)/1000000.);
// rounding by int(.+0.5) does ugly things with negative numbers
    	        opt<<" ";opt.width(9);opt<<truncated;
            }
        opt<<"\n";
    }
    if (_logdist.compare("identity")==0)    opt<<"\ndistance:\n"; else  opt<<"\nLn(distance):\n";
    for (size_t s1=1; s1<nb_sam_migf; s1++) {
        for (size_t s2=0; s2<s1; s2++) {
            if (isnan(m[s1][s2])) opt<<"     -    ";
//              if (isinf(m[s1][s2])) opt<<"     -    ";
            else {opt<<" ";opt.width(9);opt<<double(m[s1][s2]);}
            if(!isnan(m[s2][s1])&&!isnan(m[s1][s2])) {
                gra<<double(m[s1][s2])<<" "<<double(m[s2][s1])<<endl;
//                cout<<double(m[s1][s2])<<" "<<double(m[s2][s1])<<endl;
            }
        }
        opt<<"\n";
    }
    opt.close();
    gra.close();
}  /*  end, writedat  */

void writeGraOnly(const char nom_fich_mig[])	{
using namespace datamatrix;
using namespace NS_F_est;
  string graphname=nom_fich_mig;
  graphname+=".GRA";
  ofstream gra(graphname.c_str()); // gra file
  for (size_t s1=1; s1<nb_sam_migf; s1++) {
	 for (size_t s2=0; s2<s1; s2++) {
  	    if(!isnan(data[s2][s1])&&!isnan(data[s1][s2])) gra<<double(data[s1][s2])<<" "<<double(data[s2][s1])<<endl;
//	 	if(!isnan(data[s2][s1])&&!isinf(data[s1][s2])) gra<<data[s1][s2]<<" "<<data[s2][s1]<<endl;
	 }
  }
  gra.close();
}  /*  end, writeGraOnly  */



void writepma()	{
//140206: write distance matrix for Phylip
using namespace datamatrix;
using namespace NS_F_est;
	long double tabmin=0;
	long double truncated;
	string pmaname;
    char char_pma[]=".PMA";
    pmaname=gp_file+char_pma;
	ofstream pma(pmaname.c_str()); // pma file
	pma<<nb_sam_migf<<endl;
	for (size_t s1=0; s1<nb_sam_migf; s1++) {
		tabmin=std::min(tabmin,*min_element(data[s1].begin(),data[s1].end()));
	} // tabmin finit negatif sinon 0
	for (size_t s1=0; s1<nb_sam_migf; s1++) {
	 pma.setf(ios_base::left,ios_base::adjustfield);
	 pma.width(11); //"The species name is ten characters long, and must be padded out with blanks if shorter."
	 pma<<(fichier_genepop->pops[s1]->popName()).substr(0,10);
	 pma.setf(ios_base::right,ios_base::adjustfield);
	 for (size_t s2=0; s2<nb_sam_migf; s2++) {
	//"C++ never truncates data"... donc il va systematiquement sortir 6 chiffres signif
	//meme s'il commence e 0.00... et ea va depasser la width()
		if (s1==s2) truncated=0.0;
		else {
			if (s1>s2) truncated=data[s2][s1]; else truncated=data[s1][s2];
			truncated=floor((truncated-tabmin)*1000000.)/1000000.;
		}
	// rounding by int(.+0.5) does ugly things with negative numbers
	 	pma<<" ";pma.width(9);pma<<double(truncated);
	 }
	 pma<<"\n";
	}
	pma.close();

        noR_cout<<"Genetic distance matrix written in Phylip format in file "<<pmaname<<"."<<endl;


}  /*  end, writepma  */

void conversionFst() { // only called when isoldeFileBool or MultiMigFileBool is true
// conversion F->F/(1-F)
using namespace datamatrix;
using namespace NS_F_est;
bool indic_isnan=false;
char ch;
//int ignoredvalue;
        if (statname.compare("")==0) {

                noR_cout<<"\a\n Convert genetic distance from F to F/(1-F)?\n";
                noR_cout<<"\n Enter 'y' or  'n':\n";
                cin>>ch;
                cin.ignore();


		   if (ch=='y' || ch=='Y') statname="F/(1-F)"; else statname="identity";
        }
        if (statname.compare("F/(1-F)")==0) {
              for (size_t s1=0; s1<nb_sam_migf; s1++)
	             for (size_t s2=0; s2<s1; s2++)  {
                    if (!isnan(data[s2][s1])) { // sinon le laisse tel quel
                       if (data[s2][s1]!=1.0) data[s2][s1]/=(1.0-data[s2][s1]);
                       else {
                          data[s2][s1]=numeric_limits<double>::quiet_NaN();
                          indic_isnan=true;
                       }
                    }
                 }
        }
        if (indic_isnan) {
                noR_cout<<"(!) Some genetic distances=1 converted to missing information.\n";
                if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
        }
} /* end, conversionFst */

bool includedfn(size_t s1,size_t s2);
bool includedfn(size_t s1,size_t s2) {// function to exclude pairs from regression and Mantel tests
    if(cmp_nocase(typeSelection,"all")==0) return(true);
    else if(cmp_nocase(typeSelection,"only")==0) {
        return((poptypes[s1]==poptypes[s2]) && (poptypes[s1]==typeindex1));
    } else if(cmp_nocase(typeSelection,"inter")==0) {
        if ((poptypes[s1]==typeindex1 && poptypes[s2]==typeindex2)) return(true);
        if ((poptypes[s2]==typeindex1 && poptypes[s1]==typeindex2)) return(true);
    } else if(cmp_nocase(typeSelection,"inter_all_types")==0) {
      return(poptypes[s1]!=poptypes[s2]);
    } else if(cmp_nocase(typeSelection,"intra_all_types")==0) {
      return(poptypes[s1]==poptypes[s2]);
    }
    // poids nul si on arrive ici:
    return(false);
}

void conversionGeo() {
  // conversion of geographic distances
  using namespace datamatrix;
  using namespace NS_F_est;
  //cout<<"begin conversionGeo "<<nb_sam_migf;getchar();
  bool indic_lognan=false;
  char ch;
  if(_logdist.compare("")==0) {
    //if ((pauseGP) && (!_perf) && (!GeometryInSettingsBool) && (_first_of_repl)) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    printf("\a\n Fit to (l) Ln(distance) or to (d) distance?\n");
    printf("\n Enter 'l' or  'd':\n");
#endif
    
    if (scanf("%s",&ch)==0) {} // avoiding Wunused-result; but no particular scanf error to check for?
    cin.ignore();
    if (ch=='D' || ch=='d') _logdist="identity"; else _logdist="log";
  }
  for (size_t s1=0; s1<nb_sam_migf; s1++) {
    for (size_t s2=0; s2<s1; s2++)  {
      if ( ! includedfn(s1,s2)) // tests for inclusion of data
        data[s1][s2]= numeric_limits<double>::quiet_NaN();
      else if (_logdist.compare("log")==0) {
        if (data[s1][s2]>0.) {
          data[s1][s2] = log(data[s1][s2]);
          //cout<<"|"<<data[s1][s2];
        } else {
          //cout<<"$"<<data[s1][s2];
          data[s1][s2] = numeric_limits<double>::quiet_NaN();
          //            data[s1][s2] = -numeric_limits<double>::infinity(); //MINUS...
          if (_first_of_repl) {
            if(!indic_lognan) {
              noR_cout<<" Pair(s) \n";
            }
            noR_cout<<"("<<s1+1<<","<<s2+1<<") ";
          }
          indic_lognan=true;
        }
      }
    }
  }
  if ((_first_of_repl) &&(indic_lognan)) {
    noR_cout<<"\n had geographic distance <= 0.";
    noR_cout<<"\n The ln(distance) will appears as '-' in these (row,column) positions";
    noR_cout<<"\n in the output matrix.";
    if (pauseGP) { noR_cout<<"\n(Return) to continue"<<endl; getchar();}
  }
  //cout<<"end conversionGeo";getchar();
} /* end, conversionGeo */



vector<double> isoldeproc(const char nom_fich_mig[])  {
using namespace datamatrix;
using namespace NS_F_est;
vector<double> output(3);
    readGGFile(nom_fich_mig);  // note that  isolde_etc writes a mig file for each bootstrap 'replicate'... then read here
    // but the geo distances are read only the first time in such mig files
    if (_first_of_repl) {
        conversionGeo(); // note that this operates typeSelection on population (habitats) types by setting some geo dist to NaN
        if (isoldeFileBool || multiMigFileBool) { // no need to write matrix files
            writeGraOnly(isolde_file.c_str()); // call this after conversion of data !
        } else {
            writedat(data,nom_fich_mig);
        }
    }
    output=calcwritecorw(); // remplit output (regression values)  // sets _first_of_repl to false
//cout<<"fin isoldeproc: "<<output[1]<<endl;getchar();
return output;
}  /*  end, isoldeproc */



/***********************************************************************/

/********************************************************************/
void MS_for_isolde() {
using namespace NS_F_est;
long int controle;

//vector<MStype>:: iterator MSit,MSini;

//MStype<MSreal> *MSit;
//MStype<MSreal> **MStableit1,*MStableit2;


//cout<<"nb_locus "<<nb_locus<<endl;getchar();
    for(size_t loc=0;loc<nb_locus;loc++)  {
        loc_MSG[loc]=0.0;
        if (_e_stat) { loc_MSP[loc]=0.0;} // memory not allocated otherwise
        controle=0;
        for(size_t ii1=0;ii1<nb_sam-1;ii1++) {//boucle sur toute les paires de pop for one locus
            for(size_t ii2=0;ii2<ii1+1;ii2++) {
                // FR->FR we need typeSelection here...
//cout<<"MSit->nnc "<<MSit->nnc<<endl;getchar();
                MStableit=&MStable[loc][ii1][ii2];
                if(MStableit->nnc>0.5) {
                    if ( ! estimDiploidBool) {
                        if (_indiv)
                            loc_MSG[loc]+=MStableit->mmsp; //20/01/11 for IBD; recall loc_MSP is not allocated in this case
                        // this is 0 or 1, ie diversity 1-q for a pair. The loop sums over n(n-1)/2 pairs
                        else
                            loc_MSG[loc]+=MStableit->mmsi; //21/01/11 for IBD; recall loc_MSP is not allocated in this case
                    } else loc_MSG[loc]+=MStableit->mmsg;  // any other case, original code
                    if (_e_stat) loc_MSP[loc]+=MStableit->mmsp;
                    controle++;
                }
            }
        }//fin boucle sur pair

        if(controle>0.0) {
            loc_MSG[loc]/=controle;
            if (_e_stat) {
                loc_MSP[loc]/=controle;
                Qpp[loc]=double(nb_pair_sam)*(2.-loc_MSG[loc]-loc_MSP[loc])/2.+nb_sam*(1.-loc_MSG[loc]/2.);
                Qpp[loc]/=(nb_pair_sam+nb_sam);
            }
        } else { // si controle=0
            if (_e_stat) {
                Qpp[loc]=0.;
            }
        }
//if (loc==1) {cout<<"loc_MSG[1]= "<<loc_MSG[1]<<endl;getchar();}
    }//fin boucle sur loc
} //end MS_for_isolde()
/********************************************************************/

void readGeoFile(const char nom_fich_geo[], vector<vector<double> >& geoHalfMat);
void readGeoFile(const char nom_fich_geo[],vector<vector<double> >& geoHalfMat) {
//cout<<"begin readGeoFile ";getchar();
//cout<<nom_fich_mig<<endl;
using namespace NS_F_est;

/*   \
		 \
	 geog \
*/
//cout<<nb_sam_migf<<" "<<nb_pair_sam_migf;getchar();
size_t s1, s2;
char ch;
double tmp;
//////vector<vector<double> > geoHalfMat=get_geoHalfMat();
//long int Nmiss=0;
    ifstream ipt(nom_fich_geo);
    while(!ipt.is_open()) {
        noR_cout<<"\n Cannot open file "<<nom_fich_geo<<". Give another input file again: ";
	 string nfs;
	 cin>>nfs;
	 cin.ignore();
	 ipt.clear();
	 ipt.open(nfs.c_str());
    }   	/*  Data file  */
    ipt.get(ch);
    if (ipt.eof()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
            cerr<<"\n This file exists but is empty. I must exit.";
            ipt.clear();
            if (cinGetOnError) cin.get();
        #endif
        genepop_exit(1, "This file exists but is empty.");
    } else while (ipt.get()!= '\n') continue; //skip a line of comments
    geoHalfMat.resize(nb_sam);
    geoHalfMat[0].resize(0); // but this first line will not be used
    if (_first_of_repl) {
        for (s1=1; s1<nb_sam; s1++) {
            geoHalfMat[s1].resize(0);
            for (s2=0; s2<s1; s2++) {
                ipt>>tmp;
                if (ipt.fail()) {
                    #ifdef COMPATIBILITYRCPP
                        // Rcpp::R
                    #else
                        cerr<<"\a\n Incomplete geographic distance matrix! Check input file.\n...the program must terminate...\n";
                        if (cinGetOnError) cin.get();
                        ipt.clear();
                    #endif
                    genepop_exit(1, "Incomplete geographic distance matrix! Check input file.");
                } else geoHalfMat[s1].push_back(tmp);
            }
        }
    }
    ipt.close();
//cout<<"\aend readGeoFile";getchar();
} /*end, readGeoFile */

int CheckWriteDistMat(const char outfile[], vector<vector<double> >& geoHalfMat);
int CheckWriteDistMat(const char outfile[], vector<vector<double> >& geoHalfMat) {
    /// sur le modele de CFichier_genepop::computeCheckWriteDistMat
#ifdef VERBOSE
	// cout<<"debut CheckWriteDistMat"<<endl<<flush;
#endif
double maxpot=0;
///////vector<vector<double> >geoHalfMat=get_geoHalfMat();
  double pot;
  ofstream stre;
  stre.open(outfile,ios::app);
  if(! stre.is_open()) {
    RnoR_cerr<<"CheckWriteDistMat cannot open file "<<outfile;
    if (cinGetOnError) cin.get();
    genepop_exit(1, "CheckWriteDistMat cannot open file ");
  } 
  if (geoDistFromGeoFile) {
     for (unsigned int s1=1;s1<geoHalfMat.size();s1++) {
       for (unsigned int s2=0;s2<s1;s2++) {
           pot=geoHalfMat[s1][s2];
           if (pot>maxpot) maxpot=pot;
           stre<<fixed <<setprecision(15)<<pot<<" "; 
       }
       stre<<endl;
     }
  }
  stre.close();
#ifdef VERBOSE
	//cout<<"fin CheckWriteDistMat"<<endl<<flush;
#endif
if (maxpot==0) return -1; else return 0;
}


/*******************************************************************/
int create_matrices(const char nom_fich_mig[]) {
//cout<<"debut create_matrices";getchar();
//equivalent de malecot.exe
/* Pour e_r (Loiselle) les expressions sont justifies par un calcul manuscrit que j'ai
photocopie dans un maximum d'endroits...*/
  using namespace NS_F_est;
  ofstream f_mig;
  f_mig.open(nom_fich_mig,ios::out);
    if( ! f_mig.is_open()) {
      RnoR_cerr<<"\n create_matrices() cannot open file '"<<nom_fich_mig<<"' for writing: ";
      RnoR_cerr<< strerror(errno) << endl;
      // if( remove(nom_fich_mig) != 0 ) {
      //   RnoR_cerr<<"\n create_matrices() failed to delete '" << nom_fich_mig << "': " << strerror(errno) << endl;
      // }
      if (cinGetOnError) cin.get();
      genepop_exit(1, "create_matrices() cannot open file for writing.");
    } //else RnoR_cerr<<"+"; //(FIXME)
    f_mig<<"From File: "<<gp_file.c_str()<<endl;
    f_mig<< fichier_genepop->pops.size()<<" populations"<<endl;
    f_mig<<"Genetic statistic ("<<statname.c_str()<<"):"<<endl; //FR->FR if no locus of appropriate ploidy, no message here. A bit more clear in .ST2

#ifdef BIGTABLES //v4.0 defines it in GenepopS.h
#ifdef RHO_FILE // not def in v4.0
//cout<<"b";getchar();
//relit_rho();//lit result.rho
#endif
MS_for_isolde();//calcul MSG de la pop pour chaque locus
#endif


#ifdef BINFILE //not def in v4.0
    if((f_binary=fopen("binary.sss","wb"))==0) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        printf("\n program cannot open file : binary.ss");
        printf("\n end of the program");
        getchar();
        #endif

        genepop_exit(1, "program cannot open file : binary.ss");
    }
    ecriture_binfile();//relit result.rho et convertit en binaire
    rewind(f_binary);

    if(fclose(f_binary)!=0)
        {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::R
            #else
            printf("\n program cannot close file binary.sss");
            printf("\n end of the program");
            getchar();
            #endif

        genepop_exit(1, "program cannot close file binary.sss");
    }

    if((f_binary=fopen("binary.sss","rb"))==0) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        printf("\n program cannot open file : binary.sss for reading");
        printf("\n end of the program");
        getchar();
        #endif

        genepop_exit(1, "program cannot open file : binary.sss for reading");
	}
#endif

    if (nb_sam<2) {
        noR_cout<<"\nOnly "<< nb_sam<<" population. No pairwise estimation.\n";
        if (pauseGP) getchar();
    } else for(popj=2;popj<=nb_sam;popj++) {
//	cout<<"\r Creating matrices for population "<<popj;
        for(popi=1;popi<popj;popi++) {
            if ((_e_stat) && (_first_of_repl)) lecture_Xi_Xj_pmoy();//calcul (Xi+Xj)pmoy$
            lecture_popi_popj();//calcul les MS (ponderees pour ABCboot) pour la paire i_j
            ecriture_pop_tot(f_mig); // writes to f_mig
        }//fin popi
        f_mig<<endl;
    }//fin else for popj

#ifdef BINFILE
    if(fclose(f_binary)!=0) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        printf("\n program cannot close file binary.sss");
        printf("\n end of the program");
        getchar();
        #endif

        genepop_exit(1, "program cannot close file binary.sss");
	}
#endif

    f_mig<<"distances:"<<endl; //FER \n rajoute le 11/05/06
    if( f_mig.is_open()) {
      f_mig.close();
    } //else RnoR_cerr<<"!"; // FIXME
    int indicnodist=0;
    vector<vector<double> >geoHalfMat; // create_matrices -> code below writes a .mig file in all cases 
    if (_first_of_repl) { // rajout 11/2012
      if (geoDistFromGeoFile) {
        readGeoFile(geoDistFile.c_str(),geoHalfMat);
        indicnodist=CheckWriteDistMat(nom_fich_mig,geoHalfMat);
      } else indicnodist=fichier_genepop->computeCheckWriteDistMat(nom_fich_mig); //-1 (means max dist is 0...) ou 0
    }
    return indicnodist;
} //end create_matrices

vector<double> calcwritecorw() {
using namespace datamatrix;
using namespace NS_F_est; //_perf
vector<double> output(3);
unsigned int s1, s2;
ofstream iso;
//   	if (!_perf) {
  		iso.open(outname.c_str(),ios::app);   /*  Output to iso file  */
  		if (!iso.is_open()) {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::R
            #else
                cerr<<"calcwritecorw() could not reopen iso file"<<endl;
                if (cinGetOnError) cin.get();
            #endif
            genepop_exit(1, "calcwritecorw() could not reopen iso file");
        }
//  	}
  long double xsum, ysum, x2sum, xysum, x, y, xcount;
  long double xbar, ybar, varx, covxy;
  const float minifloat=numeric_limits<float>::epsilon();
  xsum = ysum = x2sum = xysum = xcount = 0.0;
  if (!_perf && _first_of_repl) effacer_ecran();
  gasp:
    /** provide a mindist value **/
    if (_perf || !pauseGP) {
        if (_logdist.compare("log")==0) {
            if (mindist<=minifloat) mindist=minifloat; //this is the NOT-log transf minimal positive value
        } else {
            if (mindist<=minifloat) mindist=0; // 0 is a valid minimal default value
        }
    } // default pour perf ou batch mode si rien dans settings
    else if(_first_of_repl && (mindist<=minifloat || alwaysAskBool)) { // -1 a la declaration, a pu etre modif par settings
             	_gotoxy(0,1);

                noR_cout<<" Minimum distance between samples\n to be taken in account for regression? ";
             	_gotoxy(41,2);
             	cin>>mindist;
             	cin.ignore();
    } //sinon il a deja une valeur>0
    /** Only now we are sure to have a mindist value**/
    if(_first_of_repl) {
       if (maxdist < mindist) {
               effacer_ecran();
               #ifdef COMPATIBILITYRCPP
                   // Rcpp::R
               #else
                 printf("maxdist < mindist. Check options.\n");
               #endif

    		   goto gasp;
        }
       if (_logdist.compare("log")==0) {
          if (mindist>0) mindistorlogdist = log(mindist);
    	  else {
               effacer_ecran();
               #ifdef COMPATIBILITYRCPP
                   // Rcpp::R
               #else
                printf("Only positive minimal distance if log transformation is used, please.\n");
               #endif

    		   goto gasp;
          }
          maxdistorlogdist = log(maxdist);
       } else {
           mindistorlogdist = mindist; // case not log transfo
           maxdistorlogdist = maxdist;
       }
    }
    //else mindistorlogdist should already be set
//cout<<mindistorlogdist;getchar();
   for (s1=0; s1<nb_sam_migf; s1++)
	 for (s2=0; s2<s1; s2++) //{cout<<"("<<s1<<","<<s2<<","<<mindist<<") ";
///////////////cout<<"("<<data[s1][s2]<<","<<data[s2][s1]<<") ";getchar();
//	  if (!isnan(data[s2][s1])&& !isinf(data[s1][s2])&& data[s1][s2]>mindistorlogdist) {
	  if (!isnan(data[s2][s1])&& !isnan(data[s1][s2])&& data[s1][s2]>mindistorlogdist && data[s1][s2]<maxdistorlogdist) {
       //cout<<s1<<" "<<s2<< " "<<data[s1][s2]<<" "<<data[s2][s1];
        xcount+=1;
		x = data[s1][s2];
		y = data[s2][s1];
//cout<<" "<<x;
		xsum += x;
		ysum += y;
		x2sum += dsqr(x);
		xysum += x * y;
      } //else cout<<"("<<s1<<","<<s2<<","<<mindist<<") ";
//}
//cout<<"(!) getchar():";getchar();
//cout<<xsum<<" "<<ysum<<" "<<x2sum<<" "<<xysum<<" "<<xcount<<endl;getchar(); getchar();
  xbar = xsum / xcount;
  ybar = ysum / xcount;
  varx = x2sum / xcount - dsqr(xbar);
  covxy = xysum / xcount - xbar * ybar;
  output[2] = double(ybar);  //mean y
  output[1] = double(covxy / varx);  //b
  output[0] = double(ybar - output[1] * xbar); //a
  // if (chh=='n') fprintf(opt,"\nFitting Fst "); else
	if ((!_perf) && (_first_of_repl)) {
        iso<<"\n________________________________________________________________\n";
	    if(cmp_nocase(typeSelection,"only")==0) {
            iso<<"Only [type "<<typeindex1<<"] pairs retained ("<<double(xcount)<<" pairs)"<<endl;
        } else if(cmp_nocase(typeSelection,"inter")==0) {
            iso<<"Only [inter types "<<typeindex1<<","<<typeindex2<<"] pairs retained ("<<double(xcount)<<" pairs)"<<endl;
        }
           iso<<"\nMinimum distance between pairs in regression: "<<mindist<<"\n";
           if (maxdist< numeric_limits<double>::max()) iso<<"Maximum distance between pairs in regression: "<<maxdist<<"\n";
        if (varx<numeric_limits<double>::epsilon()) { // comparing long double to double for safety
            iso<<"\nNo variation in geographic distance. No further analysis possible\n";
            output.resize(0);
            output.resize(2,numeric_limits<double>::quiet_NaN());
        } else {
            iso<<"\nFitting "<<statname;
            if (_logdist.compare("log")==0)
            iso<<" to a + b ln(distance)\na = "<<output[0]<<", b = "<<output[1]<<"\n";
            else
            iso<<" to a + b (distance)\na = "<<output[0]<<", b = "<<output[1]<<"\n";
        }
//		iso.close();
	}
    _first_of_repl=false;
    iso.close();
//cout<<"fin calcwritecorw: "<<output[1]<<endl;getchar();
return output;
}  /*  end, calcwritecorw  */

void mantelTest(bool clearscreen, bool rankBool) { // no declaration defaults for the arguments
using namespace datamatrix;
using namespace NS_F_est;
//Mantel test
  size_t s1, s2, t, k, z, lo, up;
  MTRand aleam;
  vector<vector<long double> >selecteddata;
  long np;
  long double ssd,ssdobs,x,y;
  long double Pvalueneg,Pvaluepos;
  long double tmp; // local var for immediate use
  vector<vector<long double> > *matfortestnan;

  vector<size_t>selectedpops1(0);
  vector<size_t>selectedpops(0);
   if(cmp_nocase(typeSelection,"all")!=0 && poptypes.size()!=nb_sam_migf) {
       #ifdef COMPATIBILITYRCPP
           // Rcpp::R
       #else
            cerr<<"(!) Mantel test: Type Selection requested but length of 'poptypes' vector does not match number of samples "<<endl;
            cerr<<"(!) 'poptypes' length is "<<poptypes.size()<<" and number of samples is "<<nb_sam_migf<<endl;
            if (cinGetOnError) cin.get();
       #endif
       genepop_exit(-1, "(!) Mantel test: Type Selection requested but length of 'poptypes' vector does not match number of samples ");
  }
  if(cmp_nocase(typeSelection,"only")==0) {
        for(size_t it = 0; it < poptypes.size(); it++) {if(poptypes[it] == typeindex1) {selectedpops1.push_back(it);}}
        nb_sam_sel=selectedpops1.size();
        selecteddata.resize(nb_sam_sel);
        for (vector<vector<long double> >::iterator ii=selecteddata.begin();ii<selecteddata.end();ii++) (*ii).resize(nb_sam_sel);
        for (s1=0; s1<nb_sam_sel; s1++) for (s2=0; s2<nb_sam_sel; s2++) selecteddata[s1][s2]=data[selectedpops1[s1]][selectedpops1[s2]];
  } else  if(cmp_nocase(typeSelection,"inter")==0) {
        for(size_t it = 0; it < poptypes.size(); it++) {
            if(poptypes[it] == typeindex1) {
                selectedpops1.push_back(it);
                selectedpops.push_back(it);
            } else if(poptypes[it] == typeindex2) {selectedpops.push_back(it);}
        }
        nb_sam_sel=selectedpops.size();
        selecteddata.resize(nb_sam_sel);
        for (vector<vector<long double> >::iterator ii=selecteddata.begin();ii<selecteddata.end();ii++) (*ii).resize(nb_sam_sel);
        for (s1=0; s1<nb_sam_sel; s1++) for (s2=0; s2<nb_sam_sel; s2++) selecteddata[s1][s2]=data[selectedpops[s1]][selectedpops[s2]];
  } else {
      nb_sam_sel=nb_sam_migf;
 }
  nb_pair_sam_sel=(nb_sam_sel)*(nb_sam_sel-1)/2;
  bool checkrangeselection=false;
  vector<size_t>p(nb_sam_sel);
  vector<vector<long double> >idistmat;
  idistmat.resize(nb_sam_sel);
  for (vector<vector<long double> >::iterator ii=idistmat.begin();ii<idistmat.end();ii++)
      (*ii).resize(nb_sam_sel);
  vector<vector<long double> >indx;
  indx.resize(nb_sam_sel);
  for (vector<vector<long double> >::iterator ii=indx.begin();ii<indx.end();ii++)
      (*ii).resize(nb_sam_sel);
  aleam.seed(static_cast<unsigned long int>(mantelSeed));
  // ici pb est que boostrap -> (*estimatingFnPtr)~slope() -> datamatrix::data has been modified...

//  cout <<double(datamatrix::data[0][1]);getchar();
  if(cmp_nocase(typeSelection,"all")==0) {matfortestnan=&data;} else {matfortestnan=&selecteddata;}
  if (rankBool) {
    idxsup((*matfortestnan),indx); //creates a semi matrix of ranks (in indx) for 'x' var (equal values => noninteger ranks) which has no NaN's; NaN are treated as lowest ranks
    // no effect of mindist or maxdist here
    for (s1=0; s1<nb_sam_sel; s1++)
        for (s2=0; s2<s1; s2++){
            /// test the distance before replacing it by its rank (new post-4.4.1): only to warn the user that Mantel test and regression are on different points
            if ((*matfortestnan)[s1][s2]<mindistorlogdist || (*matfortestnan)[s1][s2]>maxdistorlogdist) checkrangeselection=true;
            /// then replaces by its rank
            idistmat[s1][s2]=indx[s2][s1]; /* INVERSION VOLONTAIRE */  /// x
        }
    idxinf((*matfortestnan),indx); //same for 'y'; again no NaN's but matfortestnan used to detect original y NaN's below
    for (s1=0; s1<nb_sam_sel; s1++)
      for (s2=s1+1; s2<nb_sam_sel; s2++)
            idistmat[s1][s2]=indx[s1][s2]; /// y
  } else { // not ranks
    // not rank: data selection based on x range !
    // devel code '4.4 - 4.4.1' deleted here
/**
A wrong approach would be that with e.g. a low maxdist a lot of x's were set to NaN.
Then x indices woudl be permuted, creating (x,y) pairs involving a non-NaN x and an y not present in the original regression.
(with, moreover, a clear "bias" in the test as the ECDF of resampled x would include larger distances than maxDist, and the test statistic assumes that the ECDFs of both x and y are fixed.

The more general problems is that we cannot select permutable units by a one-to-one correspondance with properties of pairs of units.
Most importantly, it's not clear which distinct processes can be inferred by comparing Mantel tests on different geographic spans.
So there is no motive for defining more complex selection procedures
(except in a hierarchical setting, not relevant here).

 Hence => no selection based on x range
**/
    for (s1=0; s1<nb_sam_sel; s1++)
          for (s2=0; s2<s1; s2++) {
             idistmat[s2][s1]=(*matfortestnan)[s2][s1]; /// x (upper triang)
             idistmat[s1][s2]=(*matfortestnan)[s1][s2]; /// y (lower triang)
             if (idistmat[s1][s2]<mindistorlogdist || idistmat[s1][s2]>maxdistorlogdist) checkrangeselection=true; /// check on x
          } // x
  }
/* test on regression = test on cov = test on x*y*/
  Pvalueneg = Pvaluepos = 0.0;
  ssdobs=0.0;


  if (clearscreen) effacer_ecran();
  bool needmarginals=false; /// eviter cette complication en utilisant tjrs l'algo general ?
    for (s1=0; s1<nb_sam_sel-1; s1++)
        for (s2=s1+1; s2<nb_sam_sel; s2++) {
            if (!isnan((*matfortestnan)[s1][s2])) { //(full y) there is genetic info (!isnan) and the pair is not excluded by poptypeselection => is in matfortestnan
                /// Rank Test:
                ///    idxinf/sup sorts NaN's as -inf => lowest rank and ties NaN's with NaN's => all pairs are retained according to distance, only y NaN's are removed identically in all permutations
                ///    => marginal x sum changes over permutations if there are NaN y's
                /// Not Rank: additional problem of NaN's in x if log scale. => neither marginal sum may be held constant over permutations
                x = idistmat[s1][s2];
                y = idistmat[s2][s1];
                tmp = x*y;
    //cout<<double(x)<<" "<<double(y)<<" "<<double(tmp)<<endl;
                if (! isnan(tmp)) {
                    ssdobs -= tmp;
//                npairs++;
                } else needmarginals=true;
            }
        }

        int npairs=0;
        long double sumx=0,sumy=0,sumx2=0;
        if (needmarginals) {
            for (s1=0; s1<nb_sam_sel-1; s1++)
                for (s2=s1+1; s2<nb_sam_sel; s2++) {
                    if (!isnan((*matfortestnan)[s1][s2])) { //(full y) there is genetic info (!isnan) and the pair is not excluded by poptypeselection => is in matfortestnan
                        x = idistmat[s1][s2];
                        y = idistmat[s2][s1];
                        if (! isnan(x*y)) {
                            sumx +=x;
                            sumy +=y;
                            sumx2 +=x*x;
                            npairs +=1;
                        }
                    }
                }
            ssdobs += (sumx*sumy)/npairs;
            ssdobs /= sumx2;
        }
//cout<<npairs;getchar();
//cout<<double(ssdobs);getchar();
  if (pauseGP && mantelPerms<0) { // il est >=0 seulement si on lui a dit dans settings
            noR_cout<<"\a\n Mantel test: give the number of permutations (0 to skip test): ";
            cin>>mantelPerms;
            cin.ignore();
  }
  if (mantelPerms>0) {
          noR_cout<<"Computing the Mantel test (number of permutations: "<<mantelPerms<<")...\n";
          if (checkrangeselection) {
              noR_cout<<"NOTE: Mantel procedure not meaningful on minDist/maxDist selection.\n";
              noR_cout<<"      Using all distances, not directly comparable to regression estimates.\n";
          }

		for(t=0;t<nb_sam_sel;t++){p[t]=t;}
//ostream_vector(selectedpops1,cout);getchar();
		for(np=0;np<mantelPerms;np++){
		  if(cmp_nocase(typeSelection,"inter")==0) {// even for 'inter', permutation of type 1 samples is enough
		    size_t locint=selectedpops1.size();
		    for(t=0;t<locint-1;t++){
		      k=aleam.randInt(locint-t-1)+t; //dans t+[0,locint-t-1] = [t,locint[
		      //cout<<t<<" "<<k<<" "<<selectedpops1[t]<<" "<<selectedpops1[k]<<endl;
		      z=p[selectedpops1[k]];p[selectedpops1[k]]=p[selectedpops1[t]];p[selectedpops1[t]]=z;
		      //note that p may refer to two selected types of pops but still only type 1 positions are permuted
		      //ostream_vector(p,cout);getchar();
		    }
		  } else { // "all" or "only": permutation of all elements of p
		    for(t=0;t<nb_sam_sel-1;t++){
		      k=aleam.randInt(nb_sam_sel-t-1)+t; //dans t+[0,nb_sam_sel-t-1] = [t,nb_sam_sel[
		      z=p[k];p[k]=p[t];p[t]=z;
		    }
		  }
		  ssd=0.0;
		  
		  if (needmarginals) {
		    sumx=0.;
		    sumy=0.;
		    sumx2=0.0;
		    npairs=0;
		    
		    for (s1=0; s1<nb_sam_sel-1; s1++)
		      for (s2=s1+1; s2<nb_sam_sel; s2++) {
		        if (p[s2]>p[s1]) {
		          lo=p[s1]; up=p[s2];
		        } else {lo=p[s2]; up=p[s1];}
		        y = idistmat[lo][up]; // genetic !
		        if (!isnan(y)) { //(full y) there is genetic info (!isnan) and the pair is not excluded by poptypeselection => is in matfortestnan
		          x = idistmat[s2][s1]; // geographic !
		          tmp=x*y;
		          if (! isnan(tmp)) {
		            ssd -=tmp;
		            sumx +=x;
		            sumy +=y;
		            sumx2 +=x*x;
		            npairs +=1;
		          }
		        }
		      }
		    ssd += (sumx*sumy)/npairs;
		    ssd /= sumx2;
		  } else {
		    for (s1=0; s1<nb_sam_sel; s1++)
		      for (s2=s1+1; s2<nb_sam_sel; s2++) {
		        if (p[s2]>p[s1]) {
		          lo=p[s1]; up=p[s2];
		        } else {lo=p[s2]; up=p[s1];}
		        y = idistmat[lo][up]; // genetic !
		        if (!isnan(y)) { // there is genetic info (!isnan) and the pair is not excluded by poptypeselection (=> is in matfortestnan)
		          x = idistmat[s2][s1]; // geographic !
		          ssd -= x*y;
		          //cout <<double(x)<<" "<<double(y)<<" "<<double(tmp)<<" ";
		          //                npairs++;
		        }
		      }
		  }
		  //cout<<npairs;getchar();
		  if (ssd >= ssdobs-numeric_limits<float>::epsilon()) {Pvalueneg += 1.0;}
		  if (ssd <= ssdobs+numeric_limits<float>::epsilon()) {Pvaluepos += 1.0;}
		}   /* end, FOR nperm */
        Pvalueneg=Pvalueneg/mantelPerms;
        Pvaluepos=Pvaluepos/mantelPerms;
//cout<<outname;getchar();
        ofstream opt(outname.c_str(), ios::app);    /*  Output matrix of Mhat values  */
        opt<<"\nMantel test, "<<mantelPerms<<" permutations\n";
	    if (checkrangeselection) {
	        opt<<"NOTE: Mantel procedure not meaningful on minDist/maxDist selection.\n";
	        opt<<"      Using all distances, not directly comparable to regression estimates.\n";
	    }
        if (rankBool){
            opt<<" (statistic: Spearman Rank correlation coefficient)\n";
            opt<<" Rank-based tests may be discrepant with CIs, particularly in small samples.\n";
        } else {
            // no messgae, although discrepancies between bootstrap and permutation are possible in all cases
        }
        opt<<"\nTest of isolation by distance (One tailed Pvalue):\n";
        opt<<" Pr(correlation > observed correlation) ="<<double(Pvaluepos)<<" under null hypothesis\n";
        opt<<"\nOther one tailed Pvalue:\n Pr(correlation < observed correlation) ="<<double(Pvalueneg)<<" under null hypothesis\n";
        opt.close();
	} /*end, IF nperm */
  }  /*  end, mantelTest  */


void idxsup(vector<vector<long double> >& locdata,
            vector<vector<long double> >& indx)
//after Press et al indexx: create a semimatrix of ranks indx[][]
//for a (j>i) semi matrix (i,j)
// comme dans idxinf les NaN sont traites comme des -Inf.
// ici de toute faeon les rangs correspondants ne sont pas utilises par le test de
// Mantel; ea evite juste de donner des rangs aleatoires sinon
{
using namespace NS_F_est;
  size_t i,j,jj, l,ir,indxt,jl;
long double q;
size_t locdatasize=locdata.size();
vector<size_t>ivect(nb_pair_sam_sel+2);    //bornes ?
vector<size_t>jvect(nb_pair_sam_sel+2);
vector<size_t>jjvect(nb_pair_sam_sel+2);
vector<long double>dvect(nb_pair_sam_sel+2);
i=1;
for (j=0;j<locdatasize-1;j++) {
   for(jj=j+1;jj<locdatasize;jj++)	{
		ivect[i]=i; dvect[i]=locdata[j][jj];
		jvect[i]=j;jjvect[i]=jj;
		i++;
	 }
 } //for j
l=(nb_pair_sam_sel>>1)+1;
ir=nb_pair_sam_sel;
    for (;;) {
	    if (l>1) q=dvect[(indxt=ivect[--l])];
	    else {
		    q=dvect[(indxt=ivect[ir])];
		    ivect[ir]=ivect[1];
            if (--ir==1) {
	           ivect[1]=indxt;
// dvect[ivect[i]] contains ith value from lowest
//  TIES: (juste avant return donc ne peut perturber le reste)
	           l=1;i=1;
	           for (;i<=nb_pair_sam_sel;) {
	               while ((i+l<=nb_pair_sam_sel) &&
                            (fabs(dvect[ivect[i+l]]-dvect[ivect[i]])<numeric_limits<float>::epsilon()
                              || (isnan(dvect[ivect[i+l]])&& isnan(dvect[ivect[i]])) // autre forme de TIE
                            )
                         ) l++;
                   for(jl=0;jl<l;jl++) {
	                  indx[jvect[ivect[i+jl]]][jjvect[ivect[i+jl]]]= i+(l-1.)/2.;
                   }
                   i=i+l;l=1;
               } //for i
               return;
            } //if --ir==1
        } // (l<1) else
        i=l;
        jl=l<<1;
        while (jl<=ir){
               if (jl<ir &&
                  (dvect[ivect[jl]] < dvect[ivect[jl+1]] ||
                   (isnan(dvect[ivect[jl]]) && !(isnan(dvect[ivect[jl+1]]))) // isnan signif inf petit
                  )
               ) jl++;
//cout<<jl;getchar();
               if (q<dvect[ivect[jl]] || (isnan(q) && !isnan(dvect[ivect[jl]]))) {
                  ivect[i]=ivect[jl];
                  jl+= (i=jl);
               }
               else jl=ir+1;
         }
	     ivect[i]=indxt;
    } //for(;;)
} // idxsup

void idxinf(vector<vector<long double> >& locdata,
            vector<vector<long double> >& indx)
//after Press et al indexx: create a semimatrix of ranks indx[][]
//for a (j<i) semi matrix (i,j)
{
using namespace NS_F_est;
size_t l,ir,indxt,i,jl, j,jj;
long double q;
size_t locdatasize=locdata.size();
vector<size_t>ivect(nb_pair_sam_sel+2);
vector<size_t>jvect(nb_pair_sam_sel+2);
vector<size_t>jjvect(nb_pair_sam_sel+2);
vector<long double>dvect(nb_pair_sam_sel+2); // case O pas utilisee
i=1; // noter la valeur 0 dans la case zero pas utilisee
for (j=0;j<locdatasize-1;j++) { // remplit les cases 1 e...
	for(jj=j+1;jj<locdatasize;jj++)	{
		ivect[i]=i;
        dvect[i]=locdata[jj][j];
		jvect[i]=j;jjvect[i]=jj;
      i++;
	 }
 }
l=(nb_pair_sam_sel>>1)+1;
ir=nb_pair_sam_sel;
   for (;;) {
       if (l>1)
	      q=dvect[(indxt=ivect[--l])];
	   else {
		    q=dvect[(indxt=ivect[ir])];
		    ivect[ir]=ivect[1];
            if (--ir==1) {
  	           ivect[1]=indxt;
  // TIES:
	           l=1;i=1;
	           for (;i<=nb_pair_sam_sel;) {
	               while ((i+l<=nb_pair_sam_sel) &&
                            (fabs(dvect[ivect[i+l]]-dvect[ivect[i]])<numeric_limits<float>::epsilon()
                              || (isnan(dvect[ivect[i+l]])&& isnan(dvect[ivect[i]])) // autre forme de TIE
                            )
                         ) {
                             l++;
//cout<<double(dvect[ivect[i+l]])<<" "<<double(dvect[ivect[i]]);getchar();
                         }
//cout<<nb_pair_sam_sel<<" "<<l<<" "<<i;getchar();
                   for(jl=0;jl<l;jl++)
                       indx[jvect[ivect[i+jl]]][jjvect[ivect[i+jl]]]= i+(l-1.)/2.;
                   i=i+l;l=1;
               } //for i
	           return;
            }
        }
	    i=l;
	    jl=l<<1;
	    while (jl<=ir){
		           if (jl<ir &&
                      (dvect[ivect[jl]] < dvect[ivect[jl+1]] ||
                       (isnan(dvect[ivect[jl]]) && !(isnan(dvect[ivect[jl+1]]))) // isnan signif inf petit
                      )
                   ) jl++;
		           if (q<dvect[ivect[jl]] || (isnan(q) && !isnan(dvect[ivect[jl]]))) {
			          ivect[i]=ivect[jl];
			          jl+= (i=jl);
		           }
		           else jl=ir+1;
        }
	    ivect[i]=indxt;
   } //end for(;;)
} //end idxinf




void lecture_popi_popj() {
    using namespace NS_F_est;
#ifdef BIGTABLES
    //MStype* toto=new MStype;
    //vector<MStype>::iterator it;
#endif
// cout<<"debut lecture_popi_popj() "<<popi<<" "<<popj<<endl;getchar();
    MSp2P=0.0;
    MSg2Pw=0.0;
    denom_pot=0.;
    if (_e_stat) {
        sumQbij=0;
        sumQpp=0;
        sumQriQrj=0;
    }
    if (!_indiv) {
        MSi2P=0;
        MSi2Pw=0;
    }
//    lebon=long(popj-2)*(popj-1)/2+popi;//OK
//printf("lebon=%ld; nb_pair_tot=%ld",lebon,nb_pair_pop);
//getchar();


//va maintenant chercher les MS pour chaque pop (ou ind) aux differents loci dans binary.sss
//cout<<"ICI"<<nb_locus;getchar();
    for(size_t loc1=0;loc1<nb_locus;loc1++) {
#ifdef BINFILE
        if(fseek(f_binary,sizeof(struct MStype)*(nb_pair_sam*(loc1)+lebon-1),SEEK_SET)) {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::R
            #else
                printf("\n \n Binary file seek error in lecture_popi_popj.");
            #endif

        }
        if(fread(&MSbin,sizeof(struct MStype),1,f_binary)!=1) {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::R
            #else
                printf("\n \n Binary file read error in lecture_popi_popj.");
            #endif

        }
#endif
#ifdef BIGTABLES
        MSbin=MStable[loc1][popj-2][popi-1];
#endif
        if (_indiv) {//cout<<MSbin.mmsp<<" "<<MSbin.mmsi<<" "<<MSbin.mmsg<<" "<<MSbin.nnc<<" "<<loc_MSG[loc1]<<endl;getchar();
            if(MSbin.nnc>0.5) {
                MSg2Pw+=MSbin.mmsg*ABCweight[loc1];
                MSp2P+=MSbin.mmsp*ABCweight[loc1];
        //if (isnan(MSp2P)) {cout<<MSp2P<<" "<<MSbin.mmsp<<" "<<ABCweight[loc1]<<" "<<loc1<<endl;getchar();}
                if (_e_stat) {//estimates Qb for pair j!=r. MSg=1-Qw; MSp=1+Qw-2Qb => Q3=(1+Qw-MSp)/2=1.0-(MSg+MSp)/2.
                    sumQbij+=(1.0-(MSbin.mmsg+MSbin.mmsp)/2.)*ABCweight[loc1];
                    sumQpp+=Qpp[loc1]*ABCweight[loc1];
        //		sumQriQrj+=(QriQrj[loc1]/nnn[loc1])*ABCweight[loc1];
                    sumQriQrj+=houla[loc1][popi-1][popj-2]*ABCweight[loc1];
                    // (: l'aternative est de recalculer QriQrj et nnn e chaque iteration du bootstrap par lecture_Xi_Xj_pmoy;
                }
                //Note hack on loc_MSG in MS_for_isolde: it contains MSp or MSi if haploid singleGeneDiv
                denom_pot+=ABCweight[loc1]*loc_MSG[loc1];
            }
        } else { //groups
            if( (singleGeneDiv && MSbin.nnc>0.00000001)                             // this is tested independly for each locus
                 || ( ! singleGeneDiv && MSbin.nnc>1.00000001) ) {
/* //Weir's weighting
	if (MSbin.mmsi>-8) {
        MSi2P+=MSbin.mmsi*ABCweight[loc1]/MSbin.nnc;
        MSi2Pw+=MSbin.mmsi*ABCweight[loc1]*(1.-1./MSbin.nnc);
     }
     if (MSbin.mmsg>-8) MSg2Pw+=MSbin.mmsg*ABCweight[loc1];
     if (MSbin.mmsp>-8) MSp2P+=MSbin.mmsp*ABCweight[loc1]/MSbin.nnc;*/
 //Genepop's weighting
                if (MSbin.mmsi>-8) {
                    MSi2P+=MSbin.mmsi*ABCweight[loc1];
                    MSi2Pw+=MSbin.mmsi*ABCweight[loc1]*(MSbin.nnc-1);
                }
                if (MSbin.mmsg>-8) MSg2Pw+=MSbin.mmsg*ABCweight[loc1]*MSbin.nnc;
                if (MSbin.mmsp>-8) MSp2P+=MSbin.mmsp*ABCweight[loc1];
                if (singleGeneDiv) { // whether indiv or not
                     denom_pot+=ABCweight[loc1]*loc_MSG[loc1]*MSbin.nnc;  // nc weighting cf Genepop weighting above...
                    //the Rieux et al paper may have used a 'Weir's weighting' version
                } // else denom_pot is computed below
            } //ncc
        } //indiv else
    } //fin for if (sum sur loci)
    if (singleGeneDiv) { // whether indiv or not
        /** three cases here :
            Haploid Indiv (always singleGeneDiv) =>
            Haploid Pop // test dans testHaplPop IF singleGeneDiv in settings
            Diploid Pop // test dans testDiplPop IF singleGeneDiv in settings
        **/

        if (estimDiploidBool) denom_pot*=2;
    } else if (_indiv) {// but not the singleGeneDiv option; a_stat, e_stat... nc=1
        // if this were haploid indiv that would have been singleGeneDiv so this is always diploid here
        denom_pot*=2; /** test a ou e sur w2 dans testDiplInd **/
        if (_e_stat) sumQriQrj*=2.;
    } else { //traditional F/(1-F) without common denom, whether haploid or diploid
        //cout<<MSp2P<<" "<<MSi2Pw<<" "<<MSg2Pw;getchar();
        denom_pot=MSp2P+MSi2Pw+MSg2Pw;  /** test dans testDiplPop testHaplPop if NOT singleGeneDiv in settings **/
        // note that denom_pot is not zero if traditional F(1-F) but only one indiv per sample
    }
#ifdef BIGTABLES
//delete toto;
#endif
} //end lecture_popi_popj
/*********************************************************************/

/*********************************************************************/
int ecriture_pop_tot(ofstream& f_mig) {
using namespace NS_F_est;
double pot;
// denom pot may be negative in bootstrap when the only informative loci has <0 weight
    /**in all cases we have some meaningful denominator and use it. We first catch any zero value **/
    if(fabs(denom_pot)<0.000001) {
        f_mig<<"NaN                  "; // FR10/2010
    } else {//pot=(MSp2P-MSg2P)/denom_pot;//ancien estimateur
        if (_e_stat) pot=-2.*(sumQbij-sumQriQrj+sumQpp)/denom_pot; //noter sumQriQrj*=2.; a la fin de lecture_popi_popj();
        else if (_a_stat || singleGeneDiv) {
            pot=MSp2P/denom_pot-0.5;
        } else { //cout<<denom_pot<<" "<<MSi2P<<" "<<MSp2P;getchar(); //du meme ordre => donne scaling MSI ici
            pot=(MSp2P-MSi2P)/denom_pot;	//Fst WeirC84: (MSp-MSi)/(MSp+MSi)
            pot/=(1.-pot); //(MSp-MSi)/(2*MSi)
        }
        f_mig<<fixed<<setprecision(15)<<pot<<" ";
    }
return 0;
} // end ecriture_pop_tot()
/*********************************************************************/

void lecture_Xi_Xj_pmoy() { //pair ij compared to other individuals
using namespace NS_F_est;
//	vector<MStype>::iterator it;
	for(size_t loc=0;loc<nb_locus;loc++) {
		QriQrj[loc]=0.;
		nnn[loc]=0;
	}
	for(size_t r=1;r<=nb_sam;r++)
	if (r!=popi) {
			for(size_t loc=0;loc<nb_locus;loc++) {
    			if (popi>r) MSbin=MStable[loc][popi-2][r-1];
				else MSbin=MStable[loc][r-2][popi-1];
				if(MSbin.nnc>0.5) {
					QriQrj[loc]+=(1.-(MSbin.mmsg+MSbin.mmsp)/2);
					nnn[loc]+=1;
     			} //if nnc
			} //for loc
	} //if r vs popi ET for r
	for(size_t r=1;r<=nb_sam;r++)
	if (r!=popj) {
			for(size_t loc=0;loc<nb_locus;loc++) {
    			if (popj>r) MSbin=MStable[loc][popj-2][r-1];
				else MSbin=MStable[loc][r-2][popj-1];;
				if(MSbin.nnc>0.5) {
					QriQrj[loc]+=(1.-(MSbin.mmsg+MSbin.mmsp)/2);
					nnn[loc]+=1;
     			} //if nnc
			} //for loc
	} //if r vs popj ET for r
	// r=i xor r=j
	for(size_t loc=0;loc<nb_locus;loc++) {
        MSbin=MStable[loc][popj-2][popi-1];;
		if(MSbin.nnc>0.5) {
			QriQrj[loc]+=(2.-(MSbin.mmsg)); //X_i^2+X_j^2
			nnn[loc]+=2;
   		} //if nnc
   	houla[loc][popi-1][popj-2]=QriQrj[loc]/nnn[loc];	// keep for next use
//!!!
		if (houla[loc][popi-1][popj-2]<0.) {

                noR_cout<<"Value <0... in lecture_Xi_Xj_pmoy()\n";
			    noR_cout<<houla[loc][popi-1][popj-2]<<endl;getchar();
			    noR_cout<<loc<<" "<<popi-1<<" " <<popj-2<<" "<<nnn[loc]<<endl;getchar();


		}
    } //for loc
} // end lecture_Xi_Xj_pmoy

void initializeFest() {
  _logdist.clear();
  _a_stat=true;
  _e_stat=false;
  mindist=-1;
  maxdist=numeric_limits<double>::max();
  geoDistFromGeoFile=false;
  geoDistFile.clear();
  //geoHalfMat.clear();
  statname.clear();

  datamatrix::nb_sam_migf = 0;
  datamatrix::data.clear();

  NS_F_est::global_geno_nbr = 0L;
  NS_F_est::F_est_locIt = 0u;
  NS_F_est::global_pop_it = 0;
  NS_F_est::global_pop2_it = 0;
  NS_F_est::human_loc_it = 0;
  NS_F_est::nonvides = 0;
  NS_F_est::indic12all = 0;
  NS_F_est::identity = false;
  NS_F_est::_perf= false;
  NS_F_est::_indiv = false;
  NS_F_est::_pmat = false;
  NS_F_est::_first_of_repl = false;
  NS_F_est::nb_sam = 0;
  NS_F_est::nb_pair_sam = 0L;
  NS_F_est::nb_locus = 0u;
  NS_F_est::nb_sam_sel = 0;
  NS_F_est::nb_pair_sam_sel = 0L;
  NS_F_est::sumQbij=0.0;
  NS_F_est::sumQriQrj=0.0;
  NS_F_est::sumQpp=0.0;
  NS_F_est::MSp2P =0.0;
  NS_F_est::MSg2Pw =0.0;
  NS_F_est::denom_pot =0.0;
  NS_F_est::MSi2P =0.0;
  NS_F_est::MSi2Pw =0.0;
  NS_F_est::SSp=0.0;
  NS_F_est::SSi=0.0;
  NS_F_est::SSg=0.0;
  NS_F_est::MSp=0.0;
  NS_F_est::MSi=0.0;
  NS_F_est::MSg=0.0;
  NS_F_est::popi=0;
  NS_F_est::popj=0;
  //NS_F_est::selecteddata.clear();
  //NS_F_est::indx.clear();
  //NS_F_est::sfreqs.clear();
  //NS_F_est::Ntot = 0L;
  //NS_F_est::Ntot2 = 0.;
  NS_F_est::Nc=0.0;
  NS_F_est::mindistorlogdist=0.0;
  NS_F_est::maxdistorlogdist=0.0;
  NS_F_est::maxallname=0;
  NS_F_est::deuxsamp=0u;
  //NS_F_est::scgsp.clear();

  //NS_Fis_slmt::NLocHez.clear();
  //NS_Fis_slmt::MSgTotHez.clear();
  //NS_Fis_slmt::NLoc.clear();
  //NS_Fis_slmt::NtotLoc.clear();
  //NS_Fis_slmt::SSiTotLoc.clear();
  //NS_Fis_slmt::SSgTotLoc.clear();

  NS_FFF_slmt::SSiTot=0.0;
  NS_FFF_slmt::SSgTot=0.0;
  NS_FFF_slmt::MSg2P=0.0;

  //NS_tailles::tailleMoy.clear();
  NS_tailles::tailleMoyTot=0.0;

}

void cleanFest() {
  //geoHalfMat.clear();

  datamatrix::data.clear();

  //NS_F_est::selecteddata.clear();
  //NS_F_est::indx.clear();
  //NS_F_est::sfreqs.clear();
  //NS_F_est::scgsp.clear();

  //NS_Fis_slmt::NLocHez.clear();
  //NS_Fis_slmt::MSgTotHez.clear();
  //NS_Fis_slmt::NLoc.clear();
  //NS_Fis_slmt::NtotLoc.clear();
  //NS_Fis_slmt::SSiTotLoc.clear();
  //NS_Fis_slmt::SSgTotLoc.clear();

  //NS_tailles::tailleMoy.clear();

}
