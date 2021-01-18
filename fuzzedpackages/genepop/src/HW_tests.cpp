/***************************************************************************
@ F. Rousset 2005-2006

francois.rousset@umontpellier.fr

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
/* F.R. pinxit
 06/006 portage rapide de l'option 1 de Genepop
 09/006 changements cosmetiques
 10/006 revision Proba test
+
plus EM algo (option 8.1) 09/006
pour am�liorer le code
* on pourrait virer la plupart des �critures en fichiers (sauf ceux des chaines)
et la lecture des fichiers P_L
* La var effallnbr inclut l'info de la var indic
*/
#include <ctime>
#include <cmath> // std::max
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <limits>
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
#endif
#include "MersenneTwister.h"
#include "proba.h"
#include "GenepopS.h" //set_MC, print_p... au moins
#include "genepop.h"
#include "settings.h"
#include "HW_tests.h"
#include "tools.h"

using namespace std;

namespace NS_HW {
static size_t nb_sam,nb_locus; // from genepop.cpp
    static bool deficitbool,probtestbool,globtestbool,hwfilebool;
    //static vector<vector<bool> >indic; // test faisable ou non
}


namespace NS_HW2 {
static bool enumBool;
static long int compt;
static double pr,pech,ptot,pU,pmult;
static long int    hom[4];                 // : effectif des homozygotes
static long int    het[6];         // : effectif des h�t�rozygotes
static long int    al[4];                  // : effectif observ� des all�les
static double f[4];
//static vector<short int>effallnbr;
}

namespace NS_HW3 {
static float UInf, ecaUinf;
static size_t ii1,ii2,jj1,jj2; //exists fn double j1(double) qq part
static unsigned long int allele; // alleleest lu dans un fichier est ne peut donc �tre char
static unsigned long int tot;
static double Uu,*p,Uobs,logLR;
static unsigned long int **geno;
static double lr,lr2,seuil,seuil2;
static long int switches;
}

namespace NS_HW4 {
static float*** proba; // new in hardy1, delete in delete_proba.
} 

int hardymin() {
//option 1.x (
// approche un peu obsolete pour manipuler les donnees par rapport � celle de 5.x ds CT_tests.cpp
	using namespace NS_HW;
	vector<CPopulation *>::iterator p; // it�rateur sur les populations
	map<int,CAllele* >::iterator mit,mjt;
	CGenotypes genos; // g�notypes pour chaque population au locus courant
	char fileName[255]; // Px_Ly
	size_t nb_all;
	int popit,longueur,codi;
	size_t invcoding,s_alli,s_allj;
	// instanciation des structures de stockage des g�notypes

	// on traite locus par locus (boucle externe) ;
	for (unsigned int iLoc = 0; iLoc < fichier_genepop->loci.size(); iLoc ++) {
	  // it�rations sur les populations pour le locus courant
	  popit=0;
	  if ((codi=fichier_genepop->coding[iLoc])>3) { // sinon les fichier PL ne seront pas lus
	    if (codi==6) {
	      invcoding=1000;
	    } else invcoding=100;
	    for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
	      genos.clear(); // nettoyage
	      genos.fillGenotypes(iLoc, *p,fichier_genepop->coding[iLoc]);
	      sprintf(fileName,"P%d_L%d",popit+1,iLoc+1);
	      ofstream PL(fileName,ios::out);
	      if (!PL.is_open()) {
#ifdef COMPATIBILITYRCPP
	        // Rcpp::R
#else
	        cerr<<"Could not open file "<<fileName;
	        if (cinGetOnError) cin.get();
#endif
	        
	        genepop_exit(-1, "Could not open file ");}
	      PL<<"File: "<<gp_file.c_str()<<"   Pop: "<<(*p)->popName();
	      PL<<"   Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
	      nb_all=(*p)->loci[iLoc]->alleles.size();
	      PL<<nb_all<<endl; //nb alleles dans des genos diploides complets � loc ... dans pop ...
	      for (mit=(*p)->loci[iLoc]->alleles.begin();
            mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
	        s_alli=size_t(mit->second->_identif);
	        for (mjt=(*p)->loci[iLoc]->alleles.begin();mjt!=mit;mjt++) {
	          s_allj=size_t(mjt->second->_identif);
	          PL<<genos.getEffective(invcoding*s_alli+s_allj)<<" ";
	        }
	        PL<<genos.getEffective(invcoding*s_alli+s_alli)<<endl;
	      }
	      PL<<endl<<"  ";
	      for (mit=(*p)->loci[iLoc]->alleles.begin();
            mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
	        PL.width(4);PL<<mit->second->_identif;
	      }
	      PL<<endl<<"    ";
	      if (nb_all>62) longueur = 62 * 4; else longueur = int(4 * nb_all);
	      for (int k=0;k<longueur;k++) PL<<"_";
	      PL<<endl;
	      for (mit=(*p)->loci[iLoc]->alleles.begin();
            mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
	        s_alli=size_t(mit->second->_identif);
	        PL.width(2);PL<<mit->second->_identif;
	        for (mjt=(*p)->loci[iLoc]->alleles.begin();mjt!=mit;mjt++) {
	          s_allj=size_t(mjt->second->_identif);
	          PL.width(4);PL<<genos.getEffective(invcoding*s_alli+s_allj);
	        }
	        PL.width(4);PL<<genos.getEffective(invcoding*s_alli+s_alli);
	        PL<<endl;
	      }
	      popit++;
	      PL.close();
	    } //iter sur pops
	  } // if coding
	} // iter sur loc
return 0;
}

int hardy1(bool defbool,bool prbool, bool globbool,bool hwbool,string& hw_outfile) {
	using namespace NS_HW;
    using namespace NS_HW4;
/// this at entry:
    deficitbool=defbool;
    probtestbool=prbool;
    globtestbool=globbool;
    hwfilebool=hwbool;
    if (hwfilebool) {
       nb_sam=nb_locus=1; // se substitue � une lecture de fichier genepop
	} else {
        nb_sam=fichier_genepop->pops.size();
        nb_locus=fichier_genepop->loci.size();

        #ifdef COMPATIBILITYRCPP
        check_restriction(nb_locus, nb_sam);
        #endif

        hardymin(); //cree fichier P_L
    	proba=new float**[nb_sam];
    	for (size_t ii=0;ii<nb_sam;ii++) {
            proba[ii]=new float*[nb_locus];
            for (size_t jj=0;jj<nb_locus;jj++) {
                proba[ii][jj]=new float[5]; // Pvalue, SEPvalue, FisWC , nb matrices/switches, FisRH
                proba[ii][jj][3]=-1; // indic que pas de test (0 ne suffit pas car on peut faire le test avec 0 switches)
            }
        }
    }
//
    hw_outfile=gp_file+".HW"; // seule valuation var du namespace NS_HW4; sera renomm� en .D / .E / .P / .DG / .EG � la fin
return 0;
}

//------------------------------------------------------------------------------
//---------------------- Ecriture des donn�es dans le *.HW ---------------------
//------------------------------------------------------------------------------
void ecriture_sample_HW(string& hw_outfile){
using namespace NS_HW;
using namespace NS_HW4;
ofstream outfile;
outfile.open(hw_outfile.c_str());
    if  (!outfile.is_open()){
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        cerr<<"ecriture_sample_HW() cannot open "<<hw_outfile;
        if (cinGetOnError) cin.get();
        #endif

        genepop_exit(-1, "ecriture_sample_HW() cannot open ");
    } else {
        outfile<<"Genepop "<<getSetting("version")<<": Hardy-Weinberg test\n\n";
        outfile<<"File: "<<gp_file.c_str()<<" ("<<fichier_genepop->fileTitle.c_str()<<")\n\n";
        outfile<<"Number of populations detected:    "<<nb_sam<<"\n";
        outfile<<"Number of loci detected:           "<<nb_locus<<"\n";
        outfile.close();
    }
}



//------------------------------------------------------------------------------
//-  Fonction qui lit les fichiers pour d�terminer ceux analysables et comment -
//------------------------ Fichiers du type P1_L1 ------------------------------
//------------------------------------------------------------------------------
int lecture_fich_PL(bool testBool,vector<int>& effallnbr){
// testBool: lecture pour test ou non
//cerr<<"debut lecture_fich_PL()";
using namespace NS_HW;
using namespace NS_HW2;

int possible_exact=0, necess_mark=0;

bool analyzeBool;
ifstream fichier;
size_t kfi;
int n, n11, n12, n22;
string nom_fic;
string ligne;
    effallnbr.resize(nb_sam*nb_locus);
    kfi=0;
//cerr<<nb_sam<<" !lecture_fich_PL! "<<nb_locus;
    for(size_t ifi=0;ifi<nb_sam;ifi++){
        for(size_t jfi=0;jfi<nb_locus;jfi++){
            analyzeBool=false;
            if (hwfilebool) {
               nom_fic=hw_file; // semi matrix file
               analyzeBool=true;
            } else if (fichier_genepop->coding[jfi]>3) {
              stringstream locstst;
              locstst<<"P"<<ifi+1<<"_L"<<jfi+1;
              nom_fic=locstst.str();
              analyzeBool=true;
            }
            if (analyzeBool) {
                fichier.open(nom_fic.c_str(),ios::in);
                if  (!fichier.is_open()){
                    #ifdef COMPATIBILITYRCPP
                        // Rcpp::R
                    #else
                    cerr<<"could not open file "<<nom_fic;
                    if (cinGetOnError) cin.get();
                    #endif

                    genepop_exit(-1, "could not open file ");
                }
                else{
                    getline(fichier,ligne); // ligne de titre
                    fichier>>n; // recupere n
                    if (fichier.eof()) {
                        #ifdef COMPATIBILITYRCPP
                            // Rcpp::R
                        #else
                        cerr<<"Premature end of "<<nom_fic<<" file.\nCheck second line of input";
                        if (cinGetOnError) cin.get();
                        #endif

                        genepop_exit(-1, "Premature end of file."); }
    //cerr<<"hum "<<n;getchar();

                    n12=0; n11=0; n22=0;
                    if ( n==2 ) {
                        fichier>>n11>>n12>>n22;
                        if ( n12==1 && ( n11==0 || n22==0 ) ){
    // une fa�on pas limpide de dire que pas d'info dans le fichier
                            effallnbr[kfi]=1;
                        } else effallnbr[kfi]=2;
                    } else {
                           effallnbr[kfi]=n; // n<>2
    //cout<<endl<<kfi<<" "<<effallnbr.size()<<" "<<effallnbr[kfi]<<" "<<n;getchar();
                    }
// do nothing if effallnbr<2 !!!!!!!!!!!!!!!
                    if (effallnbr[kfi]>1) {if (effallnbr[kfi]<5) possible_exact++; else necess_mark++;}
    //cout<<endl<<kfi<<" "<<effallnbr.size()<<" "<<effallnbr[kfi]<<" "<<n;getchar();
                    fichier.close();
               } // if coding...
            } // analyzeBool
        kfi++;
        }
    }
    if (testBool && possible_exact!=0 ){
       string var;
        if (globtestbool==true) enumBool=false;
        else {
             if (alwaysAskBool || (pauseGP && enumMCindic==0)) { //Ask || (default mode et rien dans settings)

                    noR_cout<<"The complete enumeration test can be performed for some loci\n";
                    noR_cout<<"Do you want it whenever it is possible (y/n, default=";
                    if (enumMCindic==2) noR_cout<<"no"; else noR_cout<<"yes"; // enum est le d�fault si indic==1 ou 0
                    noR_cout<<") ?\n";
                    var=char(cin.get());


                 if(cmp_nocase(var,"n")==0) enumBool=false;
                 if(cmp_nocase(var,"y")==0) enumBool=true;
             } // sinon default+setting ou (batch avec ou sans setting):
             else if (enumMCindic==1) enumBool=true;  //default+setting ou (batch+setting)
             else if (enumMCindic==2) enumBool=false; //default+setting ou (batch+setting)
             else enumBool=true; // default in Batch mode; not in default mode (1e cas de ces if)
        }
/*        if(enumBool) {
            mark=necess_mark;
            exact=possible_exact;
        } else {
            mark=necess_mark+possible_exact;
            exact=0;
        }*/
    }
return 0;
}

/*------------------------------------------------------------------------------
Traitement de resultats des fichiers du type P1_L1
then writes which P_L_ files to analyze exact and MCMC computations in two vectors
then writes blabla in <data>.HW file
then writes more blabla in <data>.HW file
/------------------------------------------------------------------------------*/
void traitement_result_fichiers(vector<string>& markName,vector<string>& exacName,vector<int>& effallnbr,string& hw_outfile){
  using namespace NS_HW;
  using namespace NS_HW2;
  using namespace NS_HW4;
  ofstream outfile;
  string nom_fic;
  size_t k;
  int allnbr;
  
  exacName.resize(0);
  markName.resize(0);
  if (hwfilebool) {
    nom_fic=hw_file;
    noR_cout<<"\nInput file: "<<nom_fic<<endl;
    if (probtestbool) noR_cout<<"\nProbability test\n";
    else if (deficitbool) noR_cout<<"\nHeterozygote deficiency test\n";
    else noR_cout<<"\nHeterozygote excess test\n";
    if (enumBool && effallnbr[0]<=4) {
      exacName.push_back(nom_fic);
    }
    else {
      markName.push_back(nom_fic);
    }
  } else {
    k=0;
    for (size_t i=0;i<nb_sam;i++){
      for (size_t j=0;j<nb_locus;j++){
        stringstream locstst;
        locstst<<"P"<<i+1<<"_L"<<j+1;
        nom_fic=locstst.str();
        allnbr=effallnbr[k]; // nb d'alleles sauf �ventuellement dans le cas de 2 alleles
        if (fichier_genepop->coding[j]>3 && allnbr>1) {
          // enumBool signifie qu'on a choisi de faire le test exact whenever possible
          if (enumBool && allnbr<=4) exacName.push_back(nom_fic);
          else markName.push_back(nom_fic);
        }
        k++;
      }
    }
    nom_fic=hw_outfile; // de namespace HW_NS4; .HW
  }
  
  outfile.open(nom_fic.c_str(),ios::app);
  if  (!outfile.is_open()){
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cerr<<"traitement_result_fichiers() cannot open "<<nom_fic;
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "traitement_result_fichiers() cannot open ");
  }
  if (hwfilebool) {
    if (probtestbool) outfile<<"\n\nGenepop"<<getSetting("version")<<", Probability test:";
    else if (deficitbool) outfile<<"\n\nGenepop"<<getSetting("version")<<", Heterozygote deficiency test:";
    else outfile<<"\n\nGenepop"<<getSetting("version")<<", Heterozygote excess test:";
  }
  if (exacName.size()>0){
    outfile<<"\n\nHardy-Weinberg exact test for up to four alleles.\n (complete enumeration)";
  }
  
  if (markName.size()>0){
    set_MC_parameters(false);
    outfile<<"\n\nEstimation of exact P-Values by the Markov chain method. \n---------------------------------------------";
    if (hwfilebool) outfile<<"\nMarkov chain parameters:";
    else outfile<<"\nMarkov chain parameters for all tests:";
    outfile<<"\nDememorization:              "<<dem;
    outfile<<"\nBatches:                     "<<batchnbr;
    outfile<<"\nIterations per batch:        "<<batchlgth;
  }
  outfile.close();
  effallnbr.resize(0);    //aussi bord�lique qu'un new + delete
}

//---------------------------------------------------------------------------------------------
//--------------- calcul de la probabilit� (d'un echantillon et de chaq config dans l'enum complete) ------------
//---------------------------------------------------------------------------------------------

int calcul_proba(int nn){
using namespace NS_HW;
using namespace NS_HW2;

long i2; 
int i3;

	pr = 1;
	for (int ii1=0; ii1<nn * (nn - 1) / 2; ii1++){
		i2=het[ii1];
		if (i2!=0) {
			if (i2==1) pr=pr*2*pmult;
			else{
				for(i3=1; i3<=i2; i3++) pr=pmult*2*pr/double(i3);
			}
		}
	}

	for (int ii1=0; ii1<nn; ii1++){
		i2=hom[ii1];
		if (i2!=0) {
			if (i2==1) pr=pr*pmult;
			else{
				for(i3=1; i3<=i2; i3++) pr=pmult*pr/double(i3);
			}
		}
	}
return 0;
}

int enumeration_test(int nn, double uobs){
// attention al[0] doit �tre >0
using namespace NS_HW;
using namespace NS_HW2;


long int i, j, maxvar, imax=0, temp1, ma, mb, mc, md, all, zero, na, nb, nc;
long int C3, C4, C5, C6, C7, C8,  C10, C11, C12, DL, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, Bl, CL;
double Uu, temp2, temp;

	j=0;
	do{

		 maxvar = 0;
		 for (i=0; i<(4-j); i++){
				if ( al[i] > maxvar ){
						maxvar=al[i];
						imax=i;
				}
		 }
		 //�change des valeurs
		 temp1=al[imax];
		 al[imax]=al[3-j];
		 al[3-j]=temp1;

		 //�change des valeurs
		 temp2=f[imax];
		 f[imax]=f[3-j];
		 f[3-j]=temp2;

		 j++;
	}while (j!=4);

	ma = al[0]; mb = al[1]; mc = al[2]; md = al[3];
	ptot=0; pU=0; zero=0; compt=0;
	C3 = (md / 2) + 1;
	for(C4=1; C4<=C3; C4++){
		hom[3]=C4-1;
		DL = md - 2 * hom[3];
        C5 = ma + mb + mc;
        if ( DL>C5 ) goto fin1;
        C6 = DL + 1;
        for(C7=1; C7<=C6; C7++){
    	   C8 = DL - C7 + 1;
    	   //C9 = C7 - 1;
    	   for(C10=1; C10<=C7; C10++){
    		  C11 = C7 - C10;
              C12 = C10 - 1;
              if ( C8>ma || C11>mb || C12>mc )goto fin2;
              het[3]=C8; het[4]=C11; het[5]=C12;
              na = ma - C8; nb = mb - C11; nc = mc - C12;
              K1 = (na / 2) + 1;
              for(K2=1; K2<=K1; K2++){
      	         hom[0] = K2 - 1;
      	         all = na - 2 * hom[0];
        	     K3 = nb + nc;
      	         if (all >K3)goto fin3;
         	     K10 = all + 1;
       	         for(K3=1; K3<=K10; K3++){
         		    K4 = all - K3 + 1;
   	                K11 = K3 - 1;
           	        if (K4 > nb)goto fin4;
           	        if (K11 > nc)goto fin3;
           	        het[0] = K4;
          	        het[1] = K11;
           	        Bl = nb - K4;
           	        CL = nc - K11;
           	        K5 = Bl / 2;
           	        K6 = CL / 2;
           	        K7 = Bl - 2 * K5;
           	        if (K5 < K6) K8 = K5 + 1;
           	        else K8 = K6 + 1;
           	        for(K9=1; K9<=K8; K9++){
           		       hom[1] = K5 - K9 + 1;
             	       het[2] = K7 + 2 * K9 - 2;
             	       hom[2] = K6 - K9 + 1;
          	           calcul_proba(nn);
             	       if (pr == 0) zero = zero + 1;
             	       ptot = ptot + pr;//cumul des probas
             	       compt++;
//             	       if (compt % 10000 == 0 ) printf("%d\n",compt);
		               Uu = 0;
		               for(i=0; i<4; i++){
		                  if (al[i] != 0) Uu = Uu + double(hom[i]) / double(f[i]);
		               }
	                   if (probtestbool) {
// si je compare des probas del'ordre de 10^-10, je dois le faire en relatif...
                           temp=pech / pr - 1;
//compar relative et compar absolue de pobs >= pr dans les 2 cas      // afaire: vraiment besoin de l'absolue ?
	                       if (temp > -HW_PR_PREC || pech >= pr) pU+= pr; // proba test malgre nom pU
	                   } else {
                          if (deficitbool) {if (uobs - Uu <= HW_U_PREC) pU+= pr;}
                          else {if (Uu - uobs <= HW_U_PREC) pU+= pr;}
                       }

           	}
fin4: ; }
fin3: ; }
fin2: ; }
    }
fin1: ; }

return 0;
}


//------------------------------------------------------------------------------
//-------------- main de Fis_D : divers affichages + enumeration  ------------
//------------------------------------------------------------------------------


int enum_test_et_affich(vector<string>& exacName){
//cerr<<"debut enum_test_et_affich()";
using namespace NS_HW;
using namespace NS_HW2;
using namespace NS_HW4;

//FILE *fichier3;
ifstream inputf;
string nom_fic;
string buf;
stringstream strstr;
string ligne;
long deuxn;
int n, nn, i, j, k, inc,valpop,vallocus;
string::size_type pos;
double uobs,ptestU;
time_t midTime,_EndTime;
float _secs;

	time(&midTime);
    if (exacName.size()>0) {
        if (hwfilebool) {
                noR_cout<<" (Hardy-Weinberg exact test)";
        }
        else {
           _gotoxy(0,16);

           #ifdef COMPATIBILITYRCPP
               // Rcpp::R
           #else
                printf("Hardy-Weinberg exact tests:\n");
                noR_cout<<"     Already      exact tests done out of "<<exacName.size();
           #endif
        }
        for (size_t test=0;test<exacName.size();test++) { //boucle sur les tests exacts
            nom_fic=exacName[test];
//cerr<<nom_fic<< " "<<exacName.size();        getchar();
            inputf.open(nom_fic.c_str());
            if (!inputf.is_open()){
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Cannot open "<<nom_fic;
                if (cinGetOnError) cin.get();
                #endif

                genepop_exit(-1, "Cannot open file");
            }

            //on saute une ligne
            getline(inputf,ligne);
            inputf>>n; // recup�re nb alleles
            if (inputf.eof()) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Premature end of "<<nom_fic<<" file.\nCheck second line of input";
                if (cinGetOnError) cin.get();
                #endif

                genepop_exit(-1, "Premature end of file. Check second line of input"); }
            else if (inputf.fail()) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Cannot read "<<nom_fic<<" correctly.\nCheck second line of input";
                if (cinGetOnError) cin.get();
                #endif

                genepop_exit(-1, "Cannot read correctly. Check second line of input"); }

            if (n<4) nn=4;
            else nn=n;

            //Initialisation des matrices
            for(i=0;i<6;i++) het[i]=0;
            for(i=0;i<4;i++) hom[i]=0;
            for(i=0;i<4;i++) al[i]=0;

            k=0;
            deuxn=0;
            inputf>>hom[0];
            if (inputf.eof()) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Premature end of "<<nom_fic<<" file.\nCheck third line of input";
                if (cinGetOnError) cin.get();
                #endif

                genepop_exit(-1, "Premature end of file. Check third line of input"); }
            else if (inputf.fail()) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Cannot read "<<nom_fic<<" correctly.\nCheck third line of input";
                if (cinGetOnError) cin.get();
                #endif

                genepop_exit(-1, "Cannot read correctly. Check third line of input"); }
            al[0]=al[0]+2*hom[0];
            deuxn=deuxn+2*hom[0];

            for(i=1;i<n;i++){
               for(j=0;j<=i-1;j++){
                   inputf>>het[k];
                   if (inputf.eof()) {
                       #ifdef COMPATIBILITYRCPP
                           // Rcpp::R
                       #else
                       cerr<<"Premature end of "<<nom_fic<<" file.\nCheck line "<<3+i<<" of input";
                       if (cinGetOnError) cin.get();
                       #endif

                       genepop_exit(-1, "Premature end of file. Check lines of input"); }
                   else if (inputf.fail()) {
                       #ifdef COMPATIBILITYRCPP
                           // Rcpp::R
                       #else
                       cerr<<"Cannot read "<<nom_fic<<" correctly.\nCheck line "<<3+i<<" of input";
                       if (cinGetOnError) cin.get();
                       #endif

                       genepop_exit(-1, "Cannot read correctly. Check lines of input"); }
                   al[i]=al[i]+het[k];
                   al[j]=al[j]+het[k];
                   deuxn=deuxn+2*het[k];
                   k++;
               }
               inputf>>hom[i];
               if (inputf.eof()) {
                   #ifdef COMPATIBILITYRCPP
                       // Rcpp::R
                   #else
                   cerr<<"Premature end of "<<nom_fic<<" file.\nCheck line "<<3+i<<" of input";
                   if (cinGetOnError) cin.get();
                   #endif

                   genepop_exit(-1, "Premature end of file. Check lines of input"); }
               else if (inputf.fail()) {
                   #ifdef COMPATIBILITYRCPP
                       // Rcpp::R
                   #else
                   cerr<<"Cannot read "<<nom_fic<<" correctly.\nCheck line "<<3+i<<" of input";
                   if (cinGetOnError) cin.get();
                   #endif

                   genepop_exit(-1, "Cannot read correctly. Check lines of input"); }
               al[i]=al[i]+2*hom[i];
               deuxn=deuxn+2*hom[i];
            }
            inputf.close();

            uobs=0;
            for(i=0;i<n;i++){
                f[i]=(double(al[i])/double(deuxn));
                uobs = uobs + (double(hom[i])/double(f[i]));
            }

            //Calcul de probabilit� de l'�chantillon et ajustement de la variable mult
            pmult=1;
            inc=2;

            do{
                calcul_proba(nn);
                if ( pr < 1E-50 )	pmult= pmult +1;
                if ( pr > 1E+50 ) pmult = pmult / double(inc);
            }while ( (pr < 1E-50) || (pr > 1E+50) );

            pech=pr;
            enumeration_test(nn, uobs);

            //printf("%d\n",test);
            ptestU = pU / ptot;
            if (hwfilebool) {
              ofstream shw_outfile(hw_file.c_str(),ios::app);
              if(!shw_outfile.is_open()){
                RnoR_cerr<<"Error while reopening file "<<hw_file<<endl;
                genepop_exit(-1, "Error while reopening file ");
              }
              shw_outfile<<endl<<"P-value="<<ptestU<<" ("<<compt<<" matrices)\nNormal ending.\n";
              shw_outfile.close();
              noR_cout << "\nNormal ending.\nEdit the file " << hw_file << " for results.";
              if (pauseGP) getchar();
            } else {
              //            fichier3=fopen("result.hw","a");
              //          	fprintf(fichier3,"\n%s , %.17f D , %d", nom_fic.c_str(), ptestU, compt); //D: lu comme hetero %s
              //          	fclose(fichier3);
              // on va sauter la lecture de ce fichier; tout n'est pas tr�s �l�gant
              buf=nom_fic;
              pos=buf.find('_');
              strstr<<buf.substr(1,pos-1);
              strstr>>valpop;
              //			    cout<<buf<<" "<<buf.substr(1,pos-1)<<" "<<valpop;getchar();
              strstr.clear();
              buf=buf.substr(pos+1);
              pos=buf.find('_');
              strstr<<buf.substr(1,pos-1);
              strstr>>vallocus;
              //			    cout<<buf<<" "<<buf.substr(1,pos-1)<<" "<<vallocus;getchar();
              strstr.clear();
              //                if (ptestU == 0) proba[valpop-1][vallocus-1][0] = 9; //pas distinguable de 0 en enum complete ??
              //                else proba[valpop-1][vallocus-1][0] = ptestU;
              proba[valpop-1][vallocus-1][0] = float(ptestU);
              proba[valpop-1][vallocus-1][1] = -1; // indic enum complete
              proba[valpop-1][vallocus-1][3] = compt;
              time(&_EndTime);
              _secs=float(_EndTime- midTime);
              if (_secs>0) {
                _gotoxy(13,17);
                noR_cout<<test+1; midTime=_EndTime;
              }
            } // if hwfilebool...
        } // fin bcle sur test
        if (!hwfilebool) {
            _gotoxy(13,17);
            noR_cout<<exacName.size();
        }
        exacName.resize(0); // cf HWtest()
    } // if il y a des tests
    return 0;
}

double matrice(ifstream& sfich_PL);
double matrice(ifstream& sfich_PL) {
  // code d'erreur de lecture necess si on acc�de ici par HWfile setting
  using namespace NS_HW3;
  string bidon;
  //    cin.flush();
  
  getline(sfich_PL,bidon);//Ligne d'info
  if (sfich_PL.eof()) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cerr<<"Premature end of P<p>_L<l> file.\nCheck first line of input";
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "Premature end of file. Check first line of input"); }
  sfich_PL>>allele;//nb d'alleles
  p=new double[allele+1];
  geno=new unsigned long int*[allele+1];
  for (size_t ii=0;ii<=allele;ii++) geno[ii]=new unsigned long int[allele+1];
  if (sfich_PL.eof()) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cerr<<"Premature end of P<p>_L<l> file.\nCheck second line of input";
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "Premature end of file. Check second line of input"); }
  else if (sfich_PL.fail()) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cerr<<"Cannot read P<p>_L<l> correctly.\nCheck second line of input";
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "Cannot read correctly. Check second line of input"); }
  
  for (size_t c1=1;c1<=allele;c1++) {
    p[c1] = 0;
    for (size_t c2=1;c2<=allele;c2++) {
      geno[c1][c2] = 0;
      geno[c2][c1] = 0;
    }
  }
  
  for (size_t c1=1;c1<=allele;c1++) {
    for (size_t c2=1;c2<=c1;c2++) {
      sfich_PL>>geno[c1][c2];
      if (sfich_PL.eof()) {
#ifdef COMPATIBILITYRCPP
        // Rcpp::R
#else
        cerr<<"Premature end of P<p>_L<l> file.\nCheck line "<<2+c1<<" of input";
        if (cinGetOnError) cin.get();
#endif
        
        genepop_exit(-1, "Premature end of file. Check lines of input"); }
      else if (sfich_PL.fail()) {
#ifdef COMPATIBILITYRCPP
        // Rcpp::R
#else
        cerr<<"Cannot read P<p>_L<l> correctly.\nCheck line "<<2+c1<<" of input";
        if (cinGetOnError) cin.get();
#endif
        
        genepop_exit(-1, "Cannot read correctly. Check lines of input"); }
    }
  }
  
  tot = 0;
  
  for (size_t c1=1;c1<=allele;c1++) {
    for (size_t c2=1;c2<=c1;c2++) {
      tot = tot + 2*geno[c1][c2];
      geno[c2][c1] = geno[c1][c2]; //afaire: pourquoi ?
      p[c1]+=geno[c1][c2];
      p[c2]+=geno[c1][c2];
    }
  }
  
  Uobs = 0;
  
  for (size_t c1=1;c1<=allele;c1++) {
    p[c1]/=tot;
    Uobs=Uobs+double(geno[c1][c1])/p[c1];
  }
  
  return Uobs;
}

void choix();
void choix() {
using namespace NS_HW3;
size_t hasard;
//cout<<"debut choix() ";

   ii2 = size_t(alea.randExc(allele))+1;
   do {
       hasard = size_t(alea.randExc(allele))+1;
   } while (hasard==ii2);

   if (ii2<hasard) {
      ii1=ii2;
      ii2=hasard;
   }
   else ii1=hasard;

   jj2 = size_t(alea.randExc(allele))+1;
   do {
        hasard = size_t(alea.randExc(allele))+1;
   } while (hasard==jj2);

   if (jj2<hasard) {
      jj1=jj2;
      jj2=hasard;
   }
   else jj1=hasard;
}

void deuxheteroD();
void deuxheteroD() { //deux heteros redonnent deux autres heteros
using namespace NS_HW;
using namespace NS_HW3;
lr=double(geno[ii1][jj1])*geno[ii2][jj2] / ( double(geno[ii1][jj2]+1.0)*(geno[ii2][jj1]+1.0) );

//    if (crr<=1.0) crr/=2; else crr=0.5;
    if (lr<=1.0) seuil=lr/2; else seuil=0.5;
    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]-1;
      geno[ii2][jj2]=geno[ii2][jj2]-1;
      geno[ii1][jj2]=geno[ii1][jj2]+1;
      geno[ii2][jj1]=geno[ii2][jj1]+1;
      if (probtestbool) logLR+=log(lr);
    }
}

void deuxheteroR();
void deuxheteroR() { //deux heteros redonnent deux autres heteros
using namespace NS_HW;
using namespace NS_HW3;
lr=  double(geno[ii1][jj2])*geno[ii2][jj1] / ( double(geno[ii1][jj1]+1.0)*(geno[ii2][jj2]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]+1;
      geno[ii2][jj2]=geno[ii2][jj2]+1;
      geno[ii1][jj2]=geno[ii1][jj2]-1;
      geno[ii2][jj1]=geno[ii2][jj1]-1;
      if (probtestbool) logLR+=log(lr);
    }
}

void deuxhetero();
void deuxhetero() { //deux heteros redonnent deux autres heteros
using namespace NS_HW;
using namespace NS_HW3;
double uni;
lr= double(geno[ii1][jj2])*geno[ii2][jj1] / ( double(geno[ii1][jj1]+1.0)*(geno[ii2][jj2]+1.0));
lr2= double(geno[ii1][jj1])*geno[ii2][jj2]  / ( double(geno[ii1][jj2]+1.0)*(geno[ii2][jj1]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;
    if (lr2<=1.0) seuil2=lr2/2; else seuil2=0.5;
    uni=alea();

    if (uni<=seuil+seuil2) {
       switches++;
       if (uni<=seuil) {
         geno[ii1][jj1]=geno[ii1][jj1]+1;
         geno[ii2][jj2]=geno[ii2][jj2]+1;
         geno[ii1][jj2]=geno[ii1][jj2]-1;
         geno[ii2][jj1]=geno[ii2][jj1]-1;
         if (probtestbool) logLR+=log(lr);
       }
       else {
         geno[ii1][jj1]=geno[ii1][jj1]-1;
         geno[ii2][jj2]=geno[ii2][jj2]-1;
         geno[ii1][jj2]=geno[ii1][jj2]+1;
         geno[ii2][jj1]=geno[ii2][jj1]+1;
         if (probtestbool) logLR+=log(lr2);
       }
    }
}

void unhomobisD();
void unhomobisD() { //deux heteros donnent un hom et un het
using namespace NS_HW;
using namespace NS_HW3;
lr=   0.5*geno[ii1][jj1]*geno[ii2][jj2] / (double(geno[ii1][jj2]+1.0)* (geno[ii2][jj1]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]-1;
      geno[ii2][jj2]=geno[ii2][jj2]-1;
      geno[ii1][jj2]=geno[ii1][jj2]+1;
      geno[ii2][jj1]=geno[ii2][jj1]+1;
      if (probtestbool) logLR+=log(lr);
      else {
           if (ii1==jj2) Uu=Uu+(1.0/p[ii1]);
           else Uu=Uu+(1.0/p[ii2]);
      }
    }
}

void unhomobisR();
void unhomobisR() { //deux heteros donnent un hom et un het
using namespace NS_HW;
using namespace NS_HW3;
lr= ( 2.0*geno[ii1][jj2]*geno[ii2][jj1] ) / ( double(geno[ii1][jj1]+1.0)* (geno[ii2][jj2]+1.0) );
    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]+1;
      geno[ii2][jj2]=geno[ii2][jj2]+1;
      geno[ii1][jj2]=geno[ii1][jj2]-1;
      geno[ii2][jj1]=geno[ii2][jj1]-1;
      if (probtestbool) logLR+=log(lr);
      else {
           if (ii1==jj2) Uu=Uu-(1.0/p[ii1]);
           else Uu=Uu-(1.0/p[ii2]);
      }
    }
}

void unhomobis();
void unhomobis() { //deux heteros donnent un hom et un het
using namespace NS_HW;
using namespace NS_HW3;
double uni;
lr=2.0*geno[ii1][jj2]*geno[ii2][jj1] / ( double(geno[ii1][jj1]+1.0)* (geno[ii2][jj2]+1.0) );
lr2=0.5*geno[ii1][jj1]*geno[ii2][jj2] / ( double(geno[ii1][jj2]+1.0)*(geno[ii2][jj1]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;
    if (lr2<=1.0) seuil2=lr2/2; else seuil2=0.5;

    uni=alea();

    if (uni<=seuil+seuil2) {
       if (uni<=seuil) {
         switches++;
         geno[ii1][jj1]=geno[ii1][jj1]+1;
         geno[ii2][jj2]=geno[ii2][jj2]+1;
         geno[ii1][jj2]=geno[ii1][jj2]-1;
         geno[ii2][jj1]=geno[ii2][jj1]-1;
         if (probtestbool) logLR+=log(lr);
         else {
              if (ii1==jj2) Uu=Uu-(1.0/p[ii1]);
              else Uu=Uu-(1.0/p[ii2]);
         }
       }
       else {
         geno[ii1][jj1]=geno[ii1][jj1]-1;
         geno[ii2][jj2]=geno[ii2][jj2]-1;
         geno[ii1][jj2]=geno[ii1][jj2]+1;
         geno[ii2][jj1]=geno[ii2][jj1]+1;
         if (probtestbool) logLR+=log(lr2);
         else {
              if (ii1==jj2) Uu=Uu+(1.0/p[ii1]);
              else Uu=Uu+(1.0/p[ii2]);
         }
       }
    }
}

void unhomoD();
void unhomoD() {
using namespace NS_HW;
using namespace NS_HW3;
lr=  2.0*geno[ii1][jj1]*geno[ii2][jj2] / ( double(geno[ii1][jj2]+1.0)*(geno[ii2][jj1]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]-1;
      geno[ii2][jj2]=geno[ii2][jj2]-1;
      geno[ii1][jj2]=geno[ii1][jj2]+1;
      geno[ii2][jj1]=geno[ii2][jj1]+1;
      if (probtestbool) logLR+=log(lr);
      else {
           if (ii1==jj1) Uu=Uu-(1.0/p[ii1]);
           else Uu=Uu-(1.0/p[ii2]);
      }
    }
}

void unhomoR();
void unhomoR() {
using namespace NS_HW;
using namespace NS_HW3;
lr= 0.5*geno[ii1][jj2]*geno[ii2][jj1] / ( double(geno[ii1][jj1]+1.0)*(geno[ii2][jj2]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]+1;
      geno[ii2][jj2]=geno[ii2][jj2]+1;
      geno[ii1][jj2]=geno[ii1][jj2]-1;
      geno[ii2][jj1]=geno[ii2][jj1]-1;
      if (probtestbool) logLR+=log(lr);
      else {
           if (ii1==jj1) Uu=Uu+(1.0/p[ii1]);
           else Uu=Uu+(1.0/p[ii2]);
      }
   }
}

void unhomo();
void unhomo() {
using namespace NS_HW;
using namespace NS_HW3;
double uni;
lr=  0.5*geno[ii1][jj2]*geno[ii2][jj1] / ( double(geno[ii1][jj1]+1.0)*(geno[ii2][jj2]+1.0) );
lr2= 2.0*geno[ii1][jj1]*geno[ii2][jj2] / ( double(geno[ii1][jj2]+1.0)*(geno[ii2][jj1]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;
    if (lr2<=1.0) seuil2=lr2/2; else seuil2=0.5;

    uni=alea();

    if (uni<=seuil+seuil2) {
       switches++;
       if (uni<=seuil) {
         geno[ii1][jj1]=geno[ii1][jj1]+1;
         geno[ii2][jj2]=geno[ii2][jj2]+1;
         geno[ii1][jj2]=geno[ii1][jj2]-1;
         geno[ii2][jj1]=geno[ii2][jj1]-1;
         if (probtestbool) logLR+=log(lr);
         else {
              if (ii1==jj1) Uu=Uu+(1.0/p[ii1]);
              else Uu=Uu+(1.0/p[ii2]);
         }
       }
       else {
         geno[ii1][jj1]=geno[ii1][jj1]-1;
         geno[ii2][jj2]=geno[ii2][jj2]-1;
         geno[ii1][jj2]=geno[ii1][jj2]+1;
         geno[ii2][jj1]=geno[ii2][jj1]+1;
         if (probtestbool) logLR+=log(lr2);
         else {
              if (ii1==jj1) Uu=Uu-(1.0/p[ii1]);
              else Uu=Uu-(1.0/p[ii2]);
         }
       }
     }
}

void deuxhomoD();
void deuxhomoD() {
using namespace NS_HW;
using namespace NS_HW3;
lr= 4.0*geno[ii1][jj1]*geno[ii2][jj2] / ( double(geno[ii1][jj2]+2.0)*(geno[ii1][jj2]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]-1;
      geno[ii2][jj2]=geno[ii2][jj2]-1;
      geno[ii1][jj2]=geno[ii1][jj2]+2;
      if (probtestbool) logLR+=log(lr);
      else Uu=Uu-(1.0/p[ii1])-(1.0/p[ii2]);
    }
}

void deuxhomoR();
void deuxhomoR(){
using namespace NS_HW;
using namespace NS_HW3;
lr= 0.25*geno[ii1][jj2]*(geno[ii1][jj2]-1.0)  / ( double(geno[ii1][jj1]+1.0)*(geno[ii2][jj2]+1.0) );

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;

    if (alea()<=seuil) {
      switches++;
      geno[ii1][jj1]=geno[ii1][jj1]+1;
      geno[ii2][jj2]=geno[ii2][jj2]+1;
      geno[ii1][jj2]=geno[ii1][jj2]-2;
      if (probtestbool) logLR+=log(lr);
      else Uu=Uu+(1.0/p[ii1])+(1.0/p[ii2]);
    }
}

void deuxhomo();
void deuxhomo() {
using namespace NS_HW;
using namespace NS_HW3;
double uni;
lr= 0.25*geno[ii1][jj2]*(geno[ii1][jj2]-1.0) / ( double(geno[ii1][jj1]+1.0)*(geno[ii2][jj2]+1.0) );
lr2= 4.0*geno[ii1][jj1]*geno[ii2][jj2] / ( double(geno[ii1][jj2]+2.0)*(geno[ii1][jj2]+1.0));

    if (lr<=1.0) seuil=lr/2; else seuil=0.5;
    if (lr2<=1.0) seuil2=lr2/2; else seuil2=0.5;

    uni=alea();

    if (uni<=seuil+seuil2) {
       switches++;
       if (uni<=seuil) {
         geno[ii1][jj1]=geno[ii1][jj1]+1;
         geno[ii2][jj2]=geno[ii2][jj2]+1;
         geno[ii1][jj2]=geno[ii1][jj2]-2;
         if (probtestbool) logLR+=log(lr);
         else Uu=Uu+(1.0/p[ii1])+(1.0/p[ii2]);
       }
       else {
         geno[ii1][jj1]=geno[ii1][jj1]-1;
         geno[ii2][jj2]=geno[ii2][jj2]-1;
         geno[ii1][jj2]=geno[ii1][jj2]+2;
         if (probtestbool) logLR+=log(lr2);
         else Uu=Uu-(1.0/p[ii1])-(1.0/p[ii2]);
       }
     }
 }


void alonzy() {
using namespace NS_HW3;
//cout<<"debut alonzy() : ";

    choix();
//cout<<"apres choix() ";
//cout<<"apres choix "<<ii1<<" "<<jj1<<" "<<ii2<<" "<<jj2<<" "<<geno[ii1][jj1]<<" "<<geno[ii2][jj2]<<" "<<geno[ii1][jj2];getchar();

  if ((ii1==jj1) && (ii2==jj2)) {
       if ((geno[ii1][jj1]*geno[ii2][jj2])==0) {
          if (geno[ii1][jj2] < 2) return;
          else deuxhomoR();
       }
       else {
          if (geno[ii1][jj2] < 2) deuxhomoD();
          else deuxhomo();
       }
    }

    else {


          if ((ii1==jj1) || (ii2==jj2)) {
               if (geno[ii1][jj1]*geno[ii2][jj2]==0) {
                   if (geno[ii1][jj2]*geno[ii2][jj1]==0) return;
                   else unhomoR();
               }
               else {
                  if (geno[ii1][jj2]*geno[ii2][jj1]==0) unhomoD();
                  else unhomo();
               }
           }

           else {

                if ((ii1==jj2) || (ii2==jj1)) {
                   if (geno[ii1][jj1]*geno[ii2][jj2]==0) {
	                   if (geno[ii1][jj2]*geno[ii2][jj1]==0) return;
                       else unhomobisR();
                    }
                    else {
	                   if (geno[ii1][jj2]*geno[ii2][jj1]==0) unhomobisD();
	                   else unhomobis();
                    }
		        }

		        else {

			          if (geno[ii1][jj1]*geno[ii2][jj2]==0) {

			             if (geno[ii1][jj2]*geno[ii2][jj1]==0) return;
			             else deuxheteroR();
			          }
		              else {
			               if (geno[ii1][jj2]*geno[ii2][jj1]==0) deuxheteroD();
			               else deuxhetero();
		              }
                }
           }
       }

						     geno[jj1][ii1]=geno[ii1][jj1];
						     geno[jj2][ii2]=geno[ii2][jj2];
						     geno[jj2][ii1]=geno[ii1][jj2];
						     geno[jj1][ii2]=geno[ii2][jj1];
}

int dememorisation();
int dememorisation() {
using namespace NS_HW;
using namespace NS_HW2;
using namespace NS_HW3;

   for (unsigned long int de=0;de<dem;de++) {
     alonzy();
     if (probtestbool) {if (fabs(logLR)<DRIFT_PREC) logLR=0;}
     else {if (fabs(Uu-Uobs)<DRIFT_PREC) Uu=Uobs;}
     }
   switches=0;
return 0;
}


//------------------------------------------------------------------------------
//-------------- main de Wein_D: MCMC ------------
//------------------------------------------------------------------------------

int HW_Pvalues_chains(vector<string>& markName) {
  using namespace NS_HW;
  using namespace NS_HW2;
  using namespace NS_HW3;
  using namespace NS_HW4;
  char inlocname[25],inpopname[25];
  // FILE* fiche4;
  ifstream inpop,inloc,intot,infile;
  ofstream outpop,outloc,outtot;
  double U_obs=0.0,Uloc,Upop,Utot; 
  float zone;
  string buf;
  stringstream strstr;
  int valpop,vallocus;
  string::size_type pos;
  time_t initesTime,inibatTime,_EndTime=0;
  clock_t initClock,endClock;
  float _secs;
  
  if (markName.size()>0) {
    if (!globtestbool) {
      _gotoxy(0,18);
      noR_cout<<"Estimation of the exact P-value (Markov chain algorithm)";
    }
    _gotoxy(0,19);
    noR_cout<<"    Already       batches done out of "<<batchnbr<<".";
    
    if (!hwfilebool) {
      _gotoxy(0,20);
      noR_cout<<"    Already       MC tests done out of "<<markName.size()<<"."; // m�me pour le global il refait tout �a
    }
    for (size_t test=0; test<markName.size();test++) {
      string nom_ficPL;
      time(&initesTime);
      endClock=initClock=clock();
      logLR=0.0;
      UInf=0.0;
      ecaUinf=0.0;
      //        ro2=0;
      //        de=0;
      nom_ficPL=markName[test];//On r�cupere le nom du fichier P_L_
      // inelegantly extract pop and locus index from the file name:
      buf=nom_ficPL;
      pos=buf.find('_');
      strstr<<buf.substr(1,pos-1);
      strstr>>valpop;
      //			    cout<<buf<<" "<<buf.substr(1,pos-1)<<" "<<valpop;getchar();
      strstr.clear();
      buf=buf.substr(pos+1);
      pos=buf.find('_');
      strstr<<buf.substr(1,pos-1);
      strstr>>vallocus;
      //			    cout<<buf<<" "<<buf.substr(1,pos-1)<<" "<<vallocus;getchar();
      strstr.clear();
      infile.open(nom_ficPL.c_str());
      if (!infile.is_open()) {
        RnoR_cerr<<"Error opening "<<nom_ficPL<<endl;
        if (cinGetOnError) cin.get();
        genepop_exit(0, "Error opening file");
      }
      else {
        // debut du test: calcul de la stat observ�e; pour le proba test on reste � la valeur relative logLR=0
        U_obs=matrice(infile); // contient p=new...[] !
        infile.close();
      }
      if (globtestbool) {
        sprintf(inpopname,"popc%d",valpop);
        inpop.open(inpopname,ios::binary);
        sprintf(inlocname,"locc%d",vallocus);
        inloc.open(inlocname,ios::binary);
        intot.open("poploc",ios::binary);
        outpop.open("outpop",ios::binary);
        outloc.open("outloc",ios::binary);
        outtot.open("outtot",ios::binary);
        //            		    outpop<<fixed;outloc<<fixed;outtot<<fixed;
        //                        outpop.precision(8); //with <<fixed, ensures 8 digits after decimal point
        //                        outloc.precision(8); //la pr�cision doit se d�grader au cours des lectures/�critures
        //                        outtot.precision(8); //=> marge conservative lors de la comparaison � Uobs
      }
      
      //cout<<"milieu HW_Pvalues_chains()";getchar();
      
      
      Uu=U_obs;
      //cout<<"dememorisation() : ";
      dememorisation();
      
      //cout<<"bcle batches : ";
      for (unsigned long int bn=0;bn<batchnbr;bn++){
        time(&inibatTime);
        
        zone=0.0;
        if (globtestbool) {
          for (unsigned long int comp=0 ;comp<batchlgth;comp++) {
            alonzy(); //one step of the MC
            if (fabs(Uu-U_obs)<DRIFT_PREC) Uu=U_obs;
            // pour les fichiers marginaux loc pop tot, il lit la valeurs cumul�e des �chant pr�c, y ajoute celle de la chaine en cours,
            //puis �crit dans 3 fichiers de sortie qui remplaceront ensuite les fichiers d'entr�e quand on changera de locus.
            // afaire: lire par blocs de batchnumber (plus petit) pour gagner du temps-- non c'est con on est dans un batch. Trouver autre chose
            inpop.read(reinterpret_cast<char *>(&Upop),sizeof Upop);
            inloc.read(reinterpret_cast<char *>(&Uloc),sizeof Uloc);
            intot.read(reinterpret_cast<char *>(&Utot),sizeof Utot);;
            Upop+=Uu;
            Uloc+=Uu;
            Utot+=Uu;
            //cout<<Uu<<" ";
            outpop.write(reinterpret_cast<char *>(&Upop),sizeof Upop);
            outloc.write(reinterpret_cast<char *>(&Uloc),sizeof Uloc);
            outtot.write(reinterpret_cast<char *>(&Utot),sizeof Utot);
          }
          //                _gotoxy(0,20);
          //                cout<<"\r "<<bn+1<<" batches done out of "<< batchnbr<<"; for pop "<<valpop<<" and locus "<<vallocus;
        } else {// pas test U global
          for (unsigned long int comp=0 ;comp<batchlgth;comp++) {
            alonzy();
            if (probtestbool) {
              if (fabs(logLR)<DRIFT_PREC) logLR=0;
              if (logLR<COMP_PREC) zone++;
            } else {
              if (fabs(Uu-U_obs)<DRIFT_PREC) Uu=U_obs;
              if (deficitbool) {if (Uu>=U_obs) zone++;}
              else if (Uu<=Uobs) zone++;
            }
          }
          UInf=UInf+ zone/batchlgth;
          ecaUinf=ecaUinf+float(pow(zone/batchlgth,2));
          //                cout<<"\r "<<bn+1<<" batches done out of "<< batchnbr<<"; ";
        } // if globtestbool else
        endClock=clock();
        time(&_EndTime);
        _secs=float(_EndTime- inibatTime);
        if (_secs>0) { // doit etre hors if globtestBool else
          _gotoxy(13,19);
          noR_cout<<"     ";
          _gotoxy(13,19);
          noR_cout<<bn+1;
          
          inibatTime=_EndTime;
        }
      } //batchnbr
      _gotoxy(13,19);
      noR_cout<<batchnbr;
      
      //        if (test<mark-1) cout<<"\v\v";getchar();
      //cout<<"apr�s bcle batches : ";
      if (globtestbool) {
        proba[valpop-1][vallocus-1][3] = switches;
        inloc.close();
        inpop.close();
        intot.close();
        outpop.close();
        outloc.close();
        outtot.close();
        remove(inpopname); // makes outfiles infiles for the next steps
        remove(inlocname);
        remove("poploc");
        rename("outloc",inlocname);
        rename("outpop",inpopname);
        rename("outtot","poploc");
      } else {
        ecaUinf=ecaUinf-float(pow(UInf,2))/float(batchnbr);
        UInf=UInf/float(batchnbr) ;
        ecaUinf=ecaUinf/ ( float(batchnbr)*float(batchnbr-1)) ;
        ecaUinf=float(std::max(float(0),sqrt(ecaUinf))); // peut �tre n�gatif par erreur num�rique
        if ((_secs=float(_EndTime- initesTime))< 2147) _secs=float(endClock- initClock)/CLOCKS_PER_SEC;
        if (hwfilebool) {
          
          noR_cout<<endl<<"P-value="<<UInf<<"; S.E="<<ecaUinf<<" ("<<switches<<" switches)\n";
          if (switches<lowSwitchNbr) noR_cout<<"BEWARE number of switches was suspiciously low. Increase MC length?\n";
          
          
          noR_cout<<"Computation time= "<<_secs<<" s.\nNormal ending.\n";
          
          
          ofstream outfile(hw_file.c_str(),ios::app);
          if(!outfile.is_open()){
            
            noR_cout<<"Error while reopening file "<<hw_file<<endl;
            
            
            genepop_exit(-1, "Error while reopening file ");
          }
          outfile<<endl<<"P-value="<<UInf<<"; S.E="<<ecaUinf<<" ("<<switches<<" switches)\n";
          if (switches<lowSwitchNbr) outfile<<"BEWARE number of switches was suspiciously low. Increase MC length?\n";
          outfile<<"Computation time= "<<_secs<<" s.\nNormal ending.\n";
          outfile.close();
        } else {
          /*                if ( (fiche4 = fopen(nom4,"a")) == 0) {
          printf("\nError opening %s ...press any key to quit\n", nom4);
          getchar();
          exitp(-1);
        }
          fprintf(fiche4,"%s , %f %f D\n",nom_ficPL.c_str(),UInf,ecaUinf);
          fclose(fiche4);*/
          //                    if (ecaUinf == 0) proba[valpop-1][vallocus-1][1] = 9; //Pour indicateur de vrai z�ro
          //                    else
          proba[valpop-1][vallocus-1][0] = UInf;
          proba[valpop-1][vallocus-1][1] = ecaUinf;
          proba[valpop-1][vallocus-1][3] = switches;
      } // if hwfilebool esle
    } //if globtestbool else
      time(&_EndTime);
      _secs=float(_EndTime- initesTime);
      if (!hwfilebool && _secs>0) { // doit etre hors if globtestBool else
        _gotoxy(13,20);
        noR_cout<<"     ";
        _gotoxy(13,20);
        noR_cout<<test+1;
        initesTime=_EndTime;
      }
      // delete ptrs allocated in matrice();
      delete[] p;
      for (unsigned int ii=0;ii<=allele;ii++) delete[] geno[ii];
      delete[] geno;
  } // test=markName.size()
    if (!hwfilebool) {
      _gotoxy(13,20);
      noR_cout<<markName.size();
    }
    markName.resize(0); //cf HWtest
  } // if markname.size...
  return 0;
}

//---------------- translation of option 5.1: basic info, frequencies, Fis for an HWfile
int HWfile_info() {
  using namespace NS_HW3; // aussi 'tot'
  using namespace NS_HW4; 
  double *allhet=NULL,homo,hetero,sumb,sumc,U=0.0,freq;
  int homobs,heterobs;
  ifstream infile;
  ofstream outfile;
  infile.open(hw_file.c_str()); 
  if (!infile.is_open()) {
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cerr<<"HWfile_info() cannot open "<<hw_file<<endl;
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(0, "HWfile_info() cannot open ");
  } else {
    U=matrice(infile); // contient p=new...[] !
    allhet=new double[allele+1];
    for (unsigned int ii=0;ii<=allele;ii++) allhet[ii]=0;
    infile.close();
  }
  outfile.open(hw_file.c_str(),ios::app);
  if  (!outfile.is_open()){
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cerr<<"HWfile_info() cannot reopen "<<hw_file;
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "HWfile_info() cannot reopen ");
  } else {
    outfile.setf(ios_base::fixed,ios_base::floatfield);
    outfile.precision(4);
    outfile<<"\n\nGenepop"<<getSetting("version")<<", basic information:\n";
    outfile<<"'Expected' numbers of homozygotes or heterozygotes\nare computed using Levene's correction\n\n";
    outfile<<"    Genotypes  Obs.      Expected\n";
    homo=0;hetero=0;heterobs=0;homobs=0;
    for (unsigned int ii=1;ii<=allele;ii++) {
      for (unsigned int jj=1;jj<ii;jj++) {
        allhet[ii]+=geno[ii][jj];
        allhet[jj]+=geno[ii][jj];
        heterobs+=geno[ii][jj];
        outfile<<right<<setw(6)<<ii<<" , "<<left<<setw(7)<<jj;
        outfile<<left<<setw(4)<<geno[ii][jj];
        freq=tot*p[ii]*tot*p[jj] / (tot - 1);
        outfile<<right<<setw(11)<<freq<<endl;
        hetero+=freq;
      }
      homobs+=geno[ii][ii];
      outfile<<right<<setw(6)<<ii<<" , "<<left<<setw(7)<<ii;
      outfile<<left<<setw(4)<<geno[ii][ii];
      freq=p[ii]*tot * (tot*p[ii] - 1) / (2*(tot - 1));
      outfile<<right<<setw(11)<<freq<<endl;
      homo+=freq;
    }
    outfile<<endl<<endl<<"    Expected number of homozygotes  : "<<homo<<endl;
    outfile<<"    Observed number of homozygotes  : "<<homobs<<endl;
    outfile<<"    Expected number of heterozygotes: "<<hetero<<endl;
    outfile<<"    Observed number of heterozygotes: "<<heterobs<<endl<<endl;
    
    outfile<<"Fis: computed as in Weir & Cockerham (1984);\nalso as in Robertson & Hill (1984).\n";
    outfile<<endl<<endl<<"    Allele frequencies and Fis:"<<endl;
    outfile<<"    -------------------------------------------------------"<<endl;
    outfile<<"                                           Fis"<<endl;
    outfile<<"                                           ----------------"<<endl;
    outfile<<"    Allele     Sample count     Frequency   W&C      R&H"<<endl;
    sumb = 0; sumc = 0;
    for (unsigned int ii=1;ii<=allele;ii++) {
      allhet[ii]/=(tot/2);
      freq=p[ii]*tot*tot*(1-p[ii])/(tot/2)-(tot-1)*allhet[ii];
      freq/=(4 * ((tot / 2) - 1));
      sumc+= allhet[ii]/2;
      sumb+=freq;
      if(freq+ allhet[ii]>0) {
        freq = freq / (freq+ allhet[ii]/2);
        outfile<<"     "<<left<<setw(3)<<ii<<setw(11)<<" "<<setw(6)<<left<<int(p[ii]*tot+0.5);
        outfile<<setw(7)<<" "<<setw(8)<<right<<p[ii]<<"   "<<right<<setw(7)<<freq<<endl;
      }
    }
    if(sumb+sumc>0) {
      outfile<<"    "<<left<<setw(4)<<"Tot"<<setw(11)<<" "<<setw(6)<<left<<tot;
      outfile<<setw(7)<<" "<<setw(8)<<right<<" "<<"   "<<right<<setw(7)<<sumb/(sumb+sumc);
      // R&H
      U-=tot/2;
      freq = double(tot - 1) * (1 + 2*U / tot) - (tot - allele); // supp ts les alleles sont pr�sents
      freq/= (2. * (tot/2 - 1) * (allele-1));
      outfile<<"  "<<right<<setw(7)<<freq<<endl;
    }
    outfile<<"    -------------------------------------------------------"<<endl;
    outfile<<"\nNormal ending."<<endl;
    outfile.close();
  }
  // delete ptrs allocated in matrice();
  delete[] allhet;
  delete[] p;
  for (unsigned int ii=0;ii<=allele;ii++) delete[] geno[ii];
  delete[] geno;
  noR_cout<<"Normal ending."<<endl;
  noR_cout<<"Edit the file "<<hw_file<<" for results"<<endl;
  if (!perf) ZeGenepopSound();
  if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
  
  
  return 0;
}


//------------------------------------------------------------------------------
//----------------------- Lecture de tous les fichiers diploides ---------------
//------------------------------------------------------------------------------
void fic_lect(){
  using namespace NS_HW;
  using namespace NS_HW4;
  ifstream fichier;
  size_t kfi,allele, k;
  unsigned long int *hom,*het, deuxn;
  string nom_fic,ligne;
  double sumb, sumc;
  float U, freq, ene, Ft, fisOBS;
  double *al,*b,*c,*h;

  kfi = 0;
  for (size_t ifi = 0; ifi < nb_sam; ifi++){
    for(size_t jfi = 0; jfi < nb_locus; jfi++){
      if (fichier_genepop->coding[jfi]>3) {
        kfi++;
        { stringstream locstst;
          locstst<<"P"<<ifi+1<<"_L"<<jfi+1;
          nom_fic=locstst.str();
        }
        fichier.open(nom_fic.c_str());
        if (!fichier.is_open()) {   //Ouverture du fichier dont on vient de construire le nom
          RnoR_cerr << "(!) From fic_lect(): could not open " << nom_fic;
          if (cinGetOnError) cin.get();
          genepop_exit(-1, "(!) From fic_lect(): could not open ");
        }
        getline(fichier,ligne);  //Saute une ligne
        fichier>> allele; //
        
        if (allele > 1) {
          
          hom= new unsigned long int[allele];                  //
          het=new unsigned long int[allele * (allele-1) / 2];    //
          al=new double[allele];
          b=new double[allele];
          c=new double[allele];
          h=new double[allele];
          
          for (size_t i = 0; i < allele; i++) {
            
            hom[i] = 0;
            al[i] = 0;
            b[i] = 0;
            c[i] = 0;
            h[i] = 0;
            
          }
          
          for (size_t i = 0; i < (allele * (allele-1) / 2); i++) het[i] = 0;  //Boucle d'initiallisation des pointeurs
          
          k = 0; deuxn = 0;
          fichier>>hom[0];
          al[0] += 2 * hom[0];
          deuxn += 2 * hom[0];
          
          for(size_t i = 1; i < allele; i++){
            
            for(size_t j = 0; j <= i-1; j++){
              
              fichier>>het[k];
              al[i] += double(het[k]);
              al[j] += double(het[k]);
              h[i] += double(het[k]);
              h[j] += double(het[k]);
              deuxn += 2 * het[k];
              k++;
            }
            fichier>>hom[i];
            deuxn += 2 * hom[i];
            al[i] += 2 * hom[i];
          }
          fichier.close();
          sumb = 0; sumc = 0;
          if (deuxn == 2){
            proba[ifi][jfi][2] = 0;
            fisOBS = -2;
          } else {
            for (size_t i = 0; i < allele; i++){
              h[i] = h[i] / (deuxn / 2);
              c[i] = h[i] / 2;
              b[i] = (al[i] * (deuxn - al[i]) / (deuxn / 2) - (deuxn - 1) * h[i]) / (4 * ((deuxn / 2) - 1));
              sumc += c[i];
              sumb += b[i];
            }
            
            //                        if (sumb + sumc>0)
            fisOBS = float(sumb) / float(sumb + sumc);
            //                        else fisOBS=-2;
            proba[ifi][jfi][2] = fisOBS;
            U = 0;
            
            for (size_t i = 0; i < allele; i++){
              freq = float(al[i]) / float(deuxn);
              if (freq != 0) U = U + float(hom[i]) / freq;
            }
            
            ene = float(deuxn) / 2;
            U = U - ene;
            Ft = (deuxn - 1) * (1 + U / ene) - (deuxn - allele);
            //                        if ((2 * (ene - 1) * (allele - 1))>0)
            Ft = Ft / (2 * (ene - 1) * (allele - 1));
            //                        else Ft=-2;
            proba[ifi][jfi][4] = Ft;
            
          }
          delete[] al;
          delete[] b;
          delete[] c;
          delete[] h;
          delete[] hom;
          delete[] het;
        }else fichier.close(); //if allele >1
      } // if coding
    } //jfi
  } //ifi
}



//------------------------------------------------------------------------------
//------------------------ Fonction analyse_pop() ------------------------------
//------------------------------------------------------------------------------
void analyse_pop(float &chiHW, long int &ddlHW, int &HWinfini, int &HinfINFINI, float &pchi, float &nu, float &chi,string hw_outfile);
void analyse_pop(float &chiHW, long int &ddlHW, int &HWinfini, int &HinfINFINI, float &pchi, float &nu, float &chi,string hw_outfile){
using namespace NS_HW;
using namespace NS_HW4;

long int ddlHWT;
float HWTOT;
int HWinfiniT;
ofstream outfile;
outfile.open(hw_outfile.c_str(),ios::app);
    if  (!outfile.is_open()){
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        cerr<<"analyse_pop() cannot reopen "<<hw_outfile;
        if (cinGetOnError) cin.get();
        #endif

        genepop_exit(-1, "analyse_pop() cannot reopen ");
    }
    outfile<<"\n\n==========================================\n     Results by population\n==========================================\n";
    outfile<<setprecision(4);

    //Analyse par pop

    HWTOT = 0; ddlHWT = 0; HWinfiniT = 0;

    for (size_t j = 0; j < nb_sam; j++){

        HWinfini = 0;  HinfINFINI = 0;
        outfile<<"\n\nPop : "<<fichier_genepop->pops[j]->popName().c_str();
        outfile<<"\n"<<"-----------------------------------------";
        outfile<<"\n                             Fis estimates";
        outfile<<"\n                            ---------------";
        outfile<<"\nlocus       P-val   S.E.    W&C     R&H     Steps ";
        outfile<<"\n----------- ------- ------- ------- ------- ------";
        chiHW = 0; ddlHW = 0;
        outfile.setf(ios_base::fixed,ios_base::floatfield);

        for (size_t i = 0; i < nb_locus; i++){
            if (fichier_genepop->coding[i]>3) {
                outfile<<"\n"<<left<<setw(11)<<fichier_genepop->loci[i]->locusName.substr(0,10)<<" ";
                if (proba[j][i][3]>-0.5) { // test done
                    outfile<<left<<setw(7)<<proba[j][i][0]<<" ";
                    ddlHW += 2;
                    ddlHWT += 2;
                    if (proba[j][i][0] <= numeric_limits<float>::epsilon()) {
                      chiHW+=2.*log(batchnbr*batchlgth); // i.e. *-=* 2 log(bound on freq)
                        HWinfini = 1;
                    } else {
                            chiHW -= 2 * log(proba[j][i][0]);
                    }

                    if (proba[j][i][1]>-numeric_limits<float>::epsilon()) // indic PAS enum complete (moins...)
                       outfile<<left<<setw(7)<<proba[j][i][1]<<" ";
                    else outfile<<"  -     ";

                    //W&C
                    outfile<<internal<<setw(7)<<proba[j][i][2]<<" ";
                    //F de R&Hill
                    outfile<<internal<<setw(7)<<proba[j][i][4]<<" ";
                    //nombre de matrices si enumeration complete
                    if (proba[j][i][1]<-numeric_limits<float>::epsilon())
                       outfile<<setw(6)<<int(proba[j][i][3]+0.5)<<" matrices";
                    else {
                       outfile<<setw(6)<<int(proba[j][i][3]+0.5)<<" switches";
                       if (proba[j][i][3]<lowSwitchNbr) outfile<<" (low!)";
                    }
                } else outfile<<" No information. ";
            } // if coding
        }

        if (probtestbool){

            if (nb_locus != 1 && ddlHW > 2) {
                outfile<<"\n\nAll (Fisher's method):";
                if (HWinfini == 1) HWinfiniT = 1;
                HWTOT += chiHW;
                outfile<<"\n Chi2 :    ";

                if (HWinfini == 1) outfile<<" > ";
                outfile<<chiHW;
                outfile<<"\n Df   :    "<<ddlHW;
                outfile<<"\n Prob :    ";
                nu = ddlHW;
                chi = chiHW;
                chi2(pchi, nu, chi); 
                
                if (pchi == -1) outfile<<"High. sign."; // should not ocur with the new chi2() code
                else {
                  if (pchi > 0.9999) outfile<<"1";
                  else {
                    if (HWinfini == 1) outfile<<" < ";  
                    print_p(pchi,outfile, 6, false) ;
                  }
                }
                
            }
        }
    }

    if (nb_sam != 1 && nb_locus != 1 && ddlHWT > 2) {

        if (probtestbool) {

            outfile<<"\n==========================================";
            outfile<<"\n All locus, all populations ";
            outfile<<"\n==========================================";

            outfile<<"\nAll (Fisher's method) :";
            outfile<<"\n Chi2 :    ";

            if (HWinfiniT == 1) outfile<<" > ";
            outfile<<HWTOT;
            outfile<<"\n Df   :    "<<ddlHWT;
            outfile<<"\n Prob :    ";
            nu = ddlHWT;
            chi = HWTOT;
            chi2(pchi, nu, chi);  
            
            if (pchi == -1) outfile<<"High. sign.";  // should not ocur with the new chi2() code
            else {
              if (pchi > 0.9999) outfile<<"1";
              else {
                if (HWinfiniT == 1) outfile<<" < ";
                print_p(pchi,outfile, 6, false) ;
              }
            }
        }
    }

    outfile<<"\n\nNormal ending\n";
    outfile.close();
//cout<<"ecrit dans "<<hw_outfile.c_str();getchar();

    string renom_fic=gp_file;
    if (probtestbool) renom_fic+=".P"; else if (deficitbool) renom_fic+=".D"; else renom_fic+=".E";
    remove(renom_fic.c_str());
    rename(hw_outfile.c_str(),renom_fic.c_str()); // apparamment cp plutot que mv
    remove(hw_outfile.c_str());


        noR_cout << "\n\nNormal ending.\nEdit the file " << renom_fic << " for informations and global results";
        if (pauseGP) { noR_cout<<"\n(Return) to continue"<<endl; getchar();}



}


//------------------------------------------------------------------------------
//---------------- Ecriture des r�sultats dans le fichier *.HW -----------------
//------------------------------------------------------------------------------
void ecriture_result(string& hw_outfile){
using namespace NS_HW;
using namespace NS_HW4;

int HWinfini=0, HinfINFINI;
long int ddlHW=0;
float chiHW=0.0, nu, chi, pchi;
ofstream outfile;
    outfile.open(hw_outfile.c_str(),ios::app);
    if  (!outfile.is_open()){
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        cerr<<"ecriture_result() cannot open "<<hw_outfile;
        if (cinGetOnError) cin.get();
        #endif

        genepop_exit(-1, "ecriture_result() cannot open ");
    }

    if (probtestbool){
        outfile<<"\nHardy Weinberg: Probability test\n        ************************";
    } else if (deficitbool){
        outfile<<"\nHardy Weinberg test when H1= heterozygote deficit\n                         ************************";
    } else {
        outfile<<"\nHardy Weinberg test when H1= heterozygote excess";
        outfile<<"\n                         ***********************";
    }

    outfile<<"\n\n";
    outfile<<setprecision(4);

    if (nb_sam != 1) {

        outfile<<"\n==========================================\n     Results by locus\n==========================================\n";

        for (size_t i = 0; i < nb_locus; i++){
            if (fichier_genepop->coding[i]>3) {

                HWinfini = 0; HinfINFINI = 0;
                outfile<<"\n\nLocus \""<<fichier_genepop->loci[i]->locusName<<"\"\n";
                outfile<<"-----------------------------------------";
                outfile<<"\n                             Fis estimates";
                outfile<<"\n                            ---------------";
                outfile<<"\nPOP         P-val   S.E.    W&C     R&H     Steps ";
                outfile<<"\n----------- ------- ------- ------- ------- ------";
                chiHW = 0; ddlHW = 0;
                outfile.setf(ios_base::fixed,ios_base::floatfield);

                for (size_t j = 0; j < nb_sam; j++){
                    outfile<<"\n"<<left<<setw(11)<<fichier_genepop->pops[j]->popName().substr(0,10)<<" ";
                    if (proba[j][i][3]>-0.5) { // test done
                        outfile<<left<<setw(7)<<proba[j][i][0]<<" ";
                        ddlHW += 2;
                        if (proba[j][i][0] <= numeric_limits<float>::epsilon()) {
                          HWinfini = 1;
                          proba[j][i][0] = numeric_limits<float>::epsilon();
                        }
                        chiHW -= 2 * log(proba[j][i][0]);
                        
                        if (proba[j][i][1]>-numeric_limits<float>::epsilon()) // indic PAS enum complete (moins...)
                           outfile<<left<<setw(7)<<proba[j][i][1]<<" ";
                        else outfile<<"  -     ";

                        //W&C
                        outfile<<internal<<setw(7)<<proba[j][i][2]<<" ";
                        //F de R&Hill
                        outfile<<internal<<setw(7)<<proba[j][i][4]<<" ";
                        //nombre de matrices si enumeration complete
                        if (proba[j][i][1]<-numeric_limits<float>::epsilon())
                           outfile<<setw(6)<<int(proba[j][i][3]+0.5)<<" matrices";
                        else {
                           outfile<<setw(6)<<int(proba[j][i][3]+0.5)<<" switches";
                           if (proba[j][i][3]<lowSwitchNbr) outfile<<" (low!)";
                        }
                    } else outfile<<" - ";
                } //Next j
            } else {
                outfile<<"\n\nLocus \""<<fichier_genepop->loci[i]->locusName<<"\" not diploid.";
                outfile<<"\n"<<"-----------------------------------------";
            }


            if ( (nb_sam != 1) && (ddlHW > 2) && (probtestbool) ) {

                outfile<<"\n\nAll (Fisher's method):";
                outfile<<"\n Chi2:    ";
                if (HWinfini == 1) outfile<<" > ";
                outfile<<chiHW;
                outfile<<"\n Df   :    "<<ddlHW;
                outfile<<"\n Prob :    ";
                nu = ddlHW; chi = chiHW;
                chi2(pchi, nu, chi);  //Lance la fonction chi2()
                if (pchi == -1) outfile<<"High. sign.";  // should not ocur with the new chi2() code
                else {
                  if (pchi > 0.9999) outfile<<"1";
                  else{
                    if (HWinfini == 1) outfile<<" < ";
                    print_p(pchi,outfile, 6, false) ;
                  }
                }
            }
        } //Next i
    }//Fin du: si (pop != 1)
    outfile.close();
    analyse_pop(chiHW, ddlHW, HWinfini, HinfINFINI, pchi, nu, chi,hw_outfile);


}

void delete_proba() {
using namespace NS_HW;
using namespace NS_HW4;
    for (size_t ii=0;ii<nb_sam;ii++) {
        for (size_t jj=0;jj<nb_locus;jj++) delete[] proba[ii][jj];
        delete[] proba[ii];
    }
    delete[] proba;
}
//------------------------------------------------------------------------------
//------------------------------ main de Global_U (lectures et �critures diverses pour les test globaux; suite dans Wein)------------
//------------------------------------------------------------------------------

int global_U_initialize(vector<vector<bool> >& indic,vector<double>& Uloc,vector<double>& Upop) {
using namespace NS_HW;

//FILE *fichier_in,*fic_upop,*fic_uloc;
ofstream popc,locc,poploc;
ifstream fich_PL;
size_t n,n11,n12,n22,het,deuxn,longueur;
vector <double>f;
double Utot,Uobs;
vector<size_t>al,hom;
double zero;
string nom_fich_PL,buf;
//set_MC_parameters();

longueur = batchlgth * batchnbr ;



//Allocation m�moire
Upop.resize(nb_sam+1,0.0);
Uloc.resize(nb_locus,0.0);
indic.resize(nb_sam+1);
for (size_t i = 0; i < nb_sam+1; i++) {
    indic[i].resize(nb_locus+1,false);
//    for (size_t j=0;j< nb_locus+1; j++) indic[i][j]=false;
}

#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    printf("Hardy Weinberg test: multi-locus and/or multi-population test\n");
#endif

//printf("Computing multi-locus and multi-populations estimates..\n");
//printf("Reading genotypic matrix of pop     and locus    \n");

//for (size_t i=0;i<(nb_sam+1);i++) Upop[i]=0.0;
//for (size_t i=0;i<nb_locus;i++) Uloc[i]=0.0;

//lecture de toutes les matrices genotypiques pour calculer les U des marges
Utot = 0.0;

    for (size_t ii=0;ii<nb_sam; ii++) {
        for (size_t ji=0;ji<nb_locus;ji++) {
            if (fichier_genepop->coding[ji]<4) { //loci haploides
               indic[ii][ji] = false;
               goto nextPL;
            }
            //ELSE
            {
              stringstream locstst;
              locstst<<"P"<<ii+1<<"_L"<<ji+1;
              nom_fich_PL=locstst.str();
            }
            fich_PL.open(nom_fich_PL.c_str());
            if ( ! fich_PL.is_open()) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Error reading "<<nom_fich_PL<<endl;
                if (cinGetOnError) cin.get();
                #endif

                genepop_exit(-1, "Error reading file");
            }
//            printf("Reading file : %s\n",fich_PL);

            //Lecture des fichiers
            getline(fich_PL,buf); // skips first line
            fich_PL>>n;
            if (n <= 1) {
                indic[ii][ji] = false;
                fich_PL.close();
                goto nextPL; // un allele maxi
            }
            if (n == 2) {
                fich_PL>>n11;
                fich_PL>>n12;
                fich_PL>>n22;
                if ((n12 == 1) && ((n11 == 0) || (n22 == 0))) {
                    indic[ii][ji] = false;  //fait double emploi avec effallnbr; et toute cette proc�dure...
                    fich_PL.close();
                    goto nextPL; //deux alleles mais pas d'info
                }
            } //if n <=2
            // si on arrive la c'est soit n==2 et info soit n>2
            indic[ii][ji] = true;
            fich_PL.close(); // close and reopen is a simple way to handle also the case where n=2, where we have read more in the file
            fich_PL.clear();
            fich_PL.open(nom_fich_PL.c_str());
            getline(fich_PL,buf); // skips first line
            fich_PL>>n;
            hom.resize(n);
            f.resize(n);
            al.resize(n);
            for (size_t i=0;i<n;i++) {
                hom[i]=0;
                al[i]=0;
                f[i]=0;
            }
            deuxn = 0;
            fich_PL>>hom[0];
            al[0] = al[0] + (2 * hom[0]);
            deuxn = deuxn + (2 * hom[0]);
            for (size_t i=1;i<n;i++) {
                for (size_t j=0;j<=(i-1);j++) {
                    fich_PL>>het;
                    al[i] = al[i] + het;
                    al[j] = al[j] + het;
                    deuxn = deuxn + (2 * het);
                }
                fich_PL>>hom[i];
                al[i] = al[i] + (2 * hom[i]);
                deuxn = deuxn + (2 * hom[i]);
            }
            fich_PL.close();
            Uobs = 0.0;
            for (size_t i=0;i<n;i++) {
                f[i] = double(al[i])/deuxn;
                Uobs = Uobs + double(hom[i])/f[i];
            }
            Upop[ii] = Upop[ii] + Uobs;
            Uloc[ji] = Uloc[ji] + Uobs;
            Utot = Utot + Uobs;
nextPL: ;
        } //ji ?
    } //ii ?
    Upop[nb_sam] = Utot;
//creation fichiers tempo
    _gotoxy(0,15);

    noR_cout<<"Checking temporary files for multi-locus tests... please wait..  ";



    //creation des fichiers popchain"i" et locchain"j" initialement plein de zero (longueur: chaine)
    // chque valeur sera lu et incr�ment�e selon chaque �tape de chaque chien parcourues
     zero = 0;

     //pour test multi locus
    for (size_t i=0;i<nb_sam;i++) {
      stringstream locstst;
      locstst<<"popc"<<i+1;
      popc.open(locstst.str().c_str(),ios::out|ios::binary);
      if ( !popc.is_open()) {
        RnoR_cerr<<"Error creating "<<locstst.str()<<endl;
        if (cinGetOnError) cin.get();
        genepop_exit(-1, "Error creating file");
      }
      for (size_t k=0;k<longueur;k++) popc.write(reinterpret_cast<char *>(&zero),sizeof zero);
      popc.close();
    }

    //pour test multi pop
    _gotoxy(0,16);

    noR_cout<<"Checking temporary files for multi-pop tests... please wait..  ";

    for (size_t j=0;j<nb_locus;j++) {
      stringstream locstst;
      locstst<<"locc"<<j+1;
      locc.open(locstst.str().c_str(),ios::out|ios::binary);
      if ( !locc.is_open()) {
        RnoR_cerr<<"Error creating "<<locstst.str()<<endl;
        if (cinGetOnError) cin.get();
        genepop_exit(-1, "Error creating file");
      }
      for (size_t k=0;k<longueur;k++) locc.write(reinterpret_cast<char *>(&zero),sizeof zero);
      locc.close();
    }

    //pour test global
    _gotoxy(0,17);

    noR_cout<<"Checking temporary files for global test... please wait..  ";


     poploc.open("poploc",ios::out|ios::binary);
     if ( !poploc.is_open()) {
       RnoR_cerr<<"Error creating poploc";
       if (cinGetOnError) cin.get();
       genepop_exit(-1, "Error creating poploc");
    }

    for (size_t k=0;k<longueur;k++) poploc.write(reinterpret_cast<char *>(&zero),sizeof zero);
    poploc.close();


return 0;
}

vector<double> ChaineD(double plimite, fstream &fichier);
vector<double> ChaineD(double plimite, fstream &fichier) {
  using namespace NS_HW2;
  unsigned long int cas;
  double U ;
  vector<double> val_eca(2);
  double eca=0.0,freq_cas,ffreq_cas=0.0;
  cas = 0;
  plimite-=0.000001; // cf commentaire sur pr�cision dans outpop, outloc, outtot
  for (unsigned long int k=0;k<batchnbr;k++) {
    for (unsigned long int i=0;i<batchlgth;i++) {
      fichier.read(reinterpret_cast<char *>(&U),sizeof U);
      if (plimite <= U) cas++;
      //cout<<plimite<<" "<<U<<" "<<cas;getchar();
    }
    freq_cas = double(cas)/batchlgth;
    ffreq_cas = ffreq_cas + freq_cas;
    eca = eca + freq_cas*freq_cas;
    cas = 0;
  }
  eca = eca - ( ffreq_cas*ffreq_cas/batchnbr);
  ffreq_cas = ffreq_cas /batchnbr;
  eca = eca /  ( (batchnbr) * ( (batchnbr) - 1.0) ) ;
  eca = sqrt(std::max((0.0),eca));
  val_eca[0] = ffreq_cas;
  val_eca[1] = eca;
  return val_eca;
}

vector<double> ChaineE(double plimite, fstream &fichier);
vector<double> ChaineE(double plimite, fstream &fichier) {
  using namespace NS_HW2;
  unsigned long int cas;
  double U ;
  vector<double> val_eca(2);
  double eca=0.0,freq_cas,ffreq_cas=0.0;
  cas = 0;
  plimite+=0.000001;
  for (unsigned long int k=0;k<batchnbr;k++) {
    for (unsigned long int i=0;i<batchlgth;i++) {
      fichier.read(reinterpret_cast<char *>(&U),sizeof U);
      if (plimite >= U) cas++;
    }
    freq_cas = double(cas)/batchlgth;
    ffreq_cas = ffreq_cas + freq_cas;
    eca = eca + freq_cas*freq_cas;
    cas = 0;
  }
  eca = eca - ( ffreq_cas*ffreq_cas/batchnbr);
  ffreq_cas = ffreq_cas /batchnbr;
  eca = eca /  ( (batchnbr) * ( (batchnbr) - 1.0) ) ;
  eca = sqrt(std::max((0.0),eca));
  val_eca[0] = ffreq_cas;
  val_eca[1] = eca;
  return val_eca;
}


// -------------------------------------------------------------------------------
// ------------------ main de HW_G_D/E (dans l'esprit): analyse des fichiers tempos ------------------------------
// -------------------------------------------------------------------------------

int HW_Pvalues_compile(vector<vector<bool> >& indic,vector<double>& Uloc,vector<double>& Upop,string& hw_outfile) {
using namespace NS_HW;
using namespace NS_HW4; 
vector<vector<double> >valoc,valop;
double num;
int denom;
fstream fic_chaine;
ofstream fichier;
ofstream outfile;
//double plimite;
valop.resize(nb_sam+1);
for (size_t i=0;i<nb_sam+1;i++) {
    valop[i].resize(2);
    valop[i][0] = -1;
    valop[i][1] = -1;
}
valoc.resize(nb_locus);
for (size_t i=0;i<nb_locus;i++) {
    valoc[i].resize(2);
    valoc[i][0] = -1;
    valoc[i][1] = -1;
}

for (size_t i=0;i<nb_sam;i++) for (size_t j=0;j<nb_locus;j++) {if (indic[i][j]) indic[i][nb_locus]=true;}

for (size_t j=0;j<nb_locus;j++) {
    for(size_t i=0;i<nb_sam;i++) {if (indic[i][j]) indic[nb_sam][j]=true;}
    if (indic[nb_sam][j]) {indic[nb_sam][nb_locus]=true;}
}

    noR_cout<<endl;



//test sur pop:

if (nb_locus != 1) {
    for (size_t i=0;i<nb_sam;i++) {
        if ( indic[i][nb_locus]) {
          _gotoxy(0,21); //!
          noR_cout<<"Multi-locus test for pop "<<fichier_genepop->pops[i]->popName()<<" ("<<i+1<<" out of "<<nb_sam<<" )\n";
          stringstream locstst;
          locstst<<"popc"<<i+1;
          fic_chaine.open(locstst.str().c_str(), ios::in|ios::binary);
          if (deficitbool) valop[i]=ChaineD(Upop[i],fic_chaine);
          else valop[i]=ChaineE(Upop[i],fic_chaine);
          fic_chaine.close();
        }
    }
}


if (nb_sam != 1) {
   //test sur loc:
   for (size_t i=0;i<nb_locus;i++) {
       if ( indic[nb_sam][i]) {
         _gotoxy(0,22); //!
         noR_cout<<"Multi-pop   test for locus "<<fichier_genepop->loci[i]->locusName<<" ("<<i+1<<" out of "<<nb_locus<<" )\n";
         stringstream locstst;
         locstst<<"locc"<<i+1;
         fic_chaine.open(locstst.str().c_str(), ios::in|ios::binary);
         if (deficitbool) valoc[i]=ChaineD(Uloc[i],fic_chaine);
         else valoc[i]=ChaineE(Uloc[i],fic_chaine);
         fic_chaine.close();
       }
   }

   if (nb_locus != 1) {
      //Test global
      _gotoxy(0,23); //!
    noR_cout<<"Multi-locus and multi-pop test...\n";


      if ( indic[nb_sam][nb_locus]) {
           fic_chaine.open("poploc",ios::in|ios::binary);
           if (deficitbool) valop[nb_sam]=ChaineD(Upop[nb_sam],fic_chaine);
           else valop[nb_sam]=ChaineE(Upop[nb_sam],fic_chaine);
           fic_chaine.close();
      }
   }
}

Upop.resize(0);
Uloc.resize(0);
for (size_t i = 0; i < nb_sam+1; i++) indic[i].resize(0);


outfile.open(hw_outfile.c_str());
if (!outfile.is_open()) {
    #ifdef COMPATIBILITYRCPP
        // Rcpp::R
    #else
    cerr<<"HW_Pvalues_compile() cannot open "<<hw_outfile;
    if (cinGetOnError) cin.get();
    #endif

   genepop_exit(-1, "HW_Pvalues_compile() cannot open ");
}

outfile<<"Global Hardy-Weinberg tests [Score (U) test]\n";
outfile<<"\n";
//d[strlen(d)-1]=0;
outfile<<"File "<<gp_file.c_str()<<endl;
outfile<<"\n";
outfile<<"Number of populations : "<<nb_sam<<"\n";
outfile<<"Number of loci :        "<<nb_locus<<"\n";
outfile<<"\n";
outfile<<"\n";
fichier<<"Unbiased estimates of exact P-value by the Markov chain method. \n";
outfile<<"---------------------------------------------------------------\n";
outfile<<"   Markov chain parameters for all tests :\n";
outfile<<"            Dememorization :       "<<dem<<"\n";
outfile<<"            Batches :              "<<batchnbr<<"\n";
outfile<<"            Iterations per batch : "<<batchlgth<<"\n";
outfile<<"\n";
if (deficitbool)
   outfile<<"Hardy Weinberg test when H1= heterozygote deficit\n";
else
   outfile<<"Hardy Weinberg test when H1= heterozygote excess\n";
outfile<<"                         ************************\n";
outfile<<"\n";

outfile.precision(4);

if (nb_locus >1) {

      outfile<<"=============================================\n";
      outfile<<"   Results by population (test multi-locus) \n";
      outfile<<"=============================================\n";
      outfile<<"\n";
      outfile<<"Population      P-val   S.E.   switches (ave.)\n";
      outfile<<"-----------     ------  ------ --------\n";

      for (size_t ii=0;ii<nb_sam;ii++) {
          outfile<<setw(14)<<left<<fichier_genepop->pops[ii]->popName().substr(0,13);
          outfile<<"  ";
          if (valop[ii][0] == -1) outfile<<" -      \n";
          else {
              num=0.;denom=0;
              for (size_t loc=0;loc<nb_locus;loc++)
                  if (proba[ii][loc][3]>-0.5) {
                      num+=proba[ii][loc][3];
                      denom++;
                  }
              outfile<<fixed<<setw(6)<<valop[ii][0]<<"  ";
              if (valop[ii][1] == -1) outfile<<" -\n";
              else {
                  outfile<<setw(6)<<valop[ii][1];
                  outfile.precision(2);
                  outfile<<right<<" "<<num/denom<<left<<endl;
                  outfile.precision(4);
              }
          }
      }
      outfile<<"\n";
} else outfile<<"\nOnly "<<nb_locus<<" locus, no multi-locus test.\n\n";

if (nb_sam >1) {

      outfile<<"=============================================\n";
      outfile<<"   Results by locus (test multi-population)\n";
      outfile<<"=============================================\n";
      outfile<<"\n";
      outfile<<"Locus           P-val   S.E.   switches (ave.)\n";
      outfile<<"-----------     ------  ------ --------\n";

      for (size_t ii=0;ii<nb_locus;ii++) {
          outfile<<setw(14)<<left<<fichier_genepop->loci[ii]->locusName.substr(0,13);
          outfile<<"  ";
          if (valoc[ii][0] == -1) outfile<<" -      \n";
          else {
              num=0.;denom=0;
              for (size_t pop=0;pop<nb_sam;pop++)
                  if (proba[pop][ii][3]>-0.5) {
                      num+=proba[pop][ii][3];
                      denom++;
                  }
              outfile<<fixed<<setw(6)<<valoc[ii][0]<<"  ";
              if (valoc[ii][1] < -0.5) outfile<<" -    \n";
              else {
                  outfile<<setw(6)<<valoc[ii][1];
                  outfile.precision(2);
                  outfile<<right<<" "<<num/denom<<left<<endl;
                  outfile.precision(4);
              }
          }
      }
} else outfile<<"\nOnly "<<nb_sam<<" population, no multi-population test.\n\n";

if ((nb_sam!=1) && (nb_locus!=1)) {

      outfile<<"\n";
      outfile<<"=============================================\n";
      outfile<<"   Result for all locus and all populations\n";
      outfile<<"=============================================\n";
      outfile<<"\n";
      outfile<<" P-val   S.E.   switches (ave.)\n";
      outfile<<" ------  ------ --------\n";

      if (valop[nb_sam][0] == -1) outfile<<"  -    \n";
      else {
          num=0.;denom=0;
          for (size_t pop=0;pop<nb_sam;pop++)
              for (size_t loc=0;loc<nb_locus;loc++)
                  if (proba[pop][loc][3]>-0.5) {
                      num+=proba[pop][loc][3];
                      denom++;
                  }
          outfile<<" "<<fixed<<setw(6)<<valop[nb_sam][0]<<"  ";
          if (valop[nb_sam][1] < -0.5) outfile<<" -    \n";
          else {
              outfile<<setw(6)<<valop[nb_sam][1];
              outfile.precision(2);
              outfile<<right<<" "<<num/denom<<left<<endl;
              outfile.precision(4);
          }
      }
}

outfile<<"\n\n";
outfile<<"Normal ending.\n";
outfile.close();

string renom_fic=gp_file;
if (deficitbool) renom_fic+=".DG"; else renom_fic+=".EG";
remove(renom_fic.c_str());
rename(hw_outfile.c_str(),renom_fic.c_str()); // apparamment cp plutot que mv
remove(hw_outfile.c_str());

noR_cout<<"\nNormal ending.\nEdit the file "<<renom_fic<<" for informations and global results\n";
if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}


return 0;
}

int Genclean_HW() {
using namespace NS_HW;
// si un remove ne marche pas c'est peut �tre qu'on n'a pas ferm� le fichier...
//     system("dir /p");
     remove("poploc");
     for(size_t jfi=0;jfi<nb_locus;jfi++){
       stringstream locstst;
       locstst<<"locc"<<jfi+1;
       remove(locstst.str().c_str());
     }
     for(size_t ifi=0;ifi<nb_sam;ifi++){
       stringstream locstst;
       locstst<<"popc"<<ifi+1;
       remove(locstst.str().c_str());
     }
//     remove("P*_L*."); ne marche pas...
     for(size_t ifi=0;ifi<nb_sam;ifi++){
        for(size_t jfi=0;jfi<nb_locus;jfi++){
          stringstream locstst;
          locstst<<"P"<<ifi+1<<"_L"<<jfi+1;
          remove(locstst.str().c_str());
        }
     }
return 0;
}

void initializeHWtests() {

  NS_HW::nb_sam = 0;
  NS_HW::nb_locus = 0;
  NS_HW::deficitbool = false;
  NS_HW::probtestbool= false;
  NS_HW::globtestbool= false;
  NS_HW::hwfilebool= false;
  //NS_HW::indic.clear();

  NS_HW2::enumBool = false;
  //NS_HW2::markName.clear();
  //NS_HW2::exacName.clear();
  NS_HW2::compt = 0L;
  NS_HW2::pr=0.0;
  NS_HW2::pech=0.0;
  NS_HW2::ptot=0.0;
  NS_HW2::pU=0.0;
  NS_HW2::pmult=0.0;
  //NS_HW2::effallnbr.clear();

  //NS_HW5::Upop.clear();
  //NS_HW5::Uloc.clear();

  NS_HW3::UInf = 0.f;
  NS_HW3::ecaUinf = 0.f;
  NS_HW3::allele = 0u;
  NS_HW3::Uu=0.0;
  NS_HW3::logLR=0.0;
  //NS_HW3::nom_ficPL.clear();
  NS_HW3::lr=0.0;
  NS_HW3::lr2=0.0;
  NS_HW3::seuil=0.0;
  NS_HW3::seuil2=0.0;
  NS_HW3::switches=0u;

  //NS_HW4::nom_outfile.clear();


}


void cleanHWtests() {

  //NS_HW::indic.clear();

  //NS_HW2::markName.clear();
  //NS_HW2::exacName.clear();
  //NS_HW2::effallnbr.clear();

  //NS_HW4::nom_outfile.clear();

  //NS_HW5::Upop.clear();
  //NS_HW5::Uloc.clear();
}
