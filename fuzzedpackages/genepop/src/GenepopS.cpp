/*
menu ppal : chercher menu()
*/
/***************************************************************************
@ F. Rousset 2006-

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

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream> //stringstream
#include <cstring>
#include <algorithm>

// my header files
#include "GenepopS.h" // menu() !
#include "myutils.h" //commun avec migraine
#include "MersenneTwister.h"

#include "F_est.h"
#include "HW_tests.h" // all HW tests
#include "CT_tests.h"
#include "conversions.h" // ecumenicism... of course
#include "proba.h"
#include "genepop.h"
#include "settings.h"
#include "multimig.h"
#include "bootstrap.h"
#ifdef COMPATIBILITYRCPP
#include "RGenepop.h"
#endif
#include "tools.h"

#ifdef NO_MODULES
// my non-header files...
#include "F_est.cpp" // all Fst and Rho stuff
#include "HW_tests.cpp" // all HW tests
#include "CT_tests.cpp" // contingency table tests
#include "conversions.cpp" // ecumenicism... of course
#include "proba.cpp" // loi normale, chi2
#include "genepop.cpp"
#include "settings.cpp"
#include "multimig.cpp"
#include "bootstrap.cpp"
#include "RGenepop.cpp"
#include "tools.cpp"

#endif

using namespace std;

//--------------------------Variables globales----------------------------------
//------------------------------------------------------------------------------
//string version=" v4.0 (Built on "+datestring+" at "+timestring+").";
//string version="4.6";
std::string getSetting(const std::string which) {
  const std::string version="4.7.5"; // the O N L Y place to story this info.
  if (which.compare("version")==0) return(version);
  if (which.compare("default_settingsfile")==0) return("genepop.txt");
  return("unknown 'which' value");
}


//string settingsfilename="genepop.txt";
vector<vector<int> >MenuOptions;
vector<int> HWfileOptions;
vector<map<int,int> >taille;
bool genicProbaTestBool=false;
bool alleleNbrTestBool=false;
bool geneDivTestBool=false;
vector<int>sequenceGeneDivRanks;
bool identitySettingsBool=true;
bool LDprobaTestBool=false;
bool gp_fileInSettingsBool=false;
bool perf=false;
bool pauseGP=true; // active les getchar() apr?s les messages divers
bool alwaysAskBool=false; //Genepop Mode=Ask
bool HWfileBool=false,strucFileBool=false,isoldeFileBool=false,multiMigFileBool=false; // indic les 2 utilitaires sur fichiers matrice unique
bool estimDiploidBool=true,phylipBool=false,Brookfield96Bool=false,nullIgnoredBool=false,NonNullfailuresBool=false;
string 	gp_file; //in original conception, name of Genepop input file; but fichier_genepop->fileName should be used in revised code (08/2014)
string hw_file,struc_file,isolde_file; // ad hoc single matrix files
unsigned long int alea_seed=67144630; //default value, else through settings
vector<double>ABCweight;
double widthCI=0.95;
string outname; // util dans F_est
char char_tmp[]=".TMP";
char char_iso[]=".ISO";
static char char_mig[]=".MIG";
static bool first_repl; // apparait dans conversion() dans le cas _perf=true... et dans isolde_etc
MTRand alea;

//const double DRIFT_PREC=100.*numeric_limits<double>::epsilon(); //commun ? toutes les MC ? pour U (directionnel)
// afaire: voir theorie sur ces compar... le epsilon est surement trop petit.
size_t dem=1,batchlgth=1;
size_t batchnbr=1;
static char char_num[]=".NUM";
static unsigned int boucle=0; // counts nested calls to menu()
const bool liptakBool=false; // code en r?serve
static bool exit_genepop = false; //ADD by jimmy

namespace NS_GP {
//vector<bool>ploidBool; // pour l'instant pont entre isolde_etc() et F_est.cpp
//static vector<int> 			allMax; //nb alleles par global_loc_nbr
//static vector<string>			nom_locus; //tabLocus;
//static vector<vector<string> >	nom_pop; //tab_pop
//static fstream fichier;
static string fichDATE,fichTIME;
static bool logdist;  // _logdist va etre value soit par les settings, soit par logdist (pour perf), soit manuellement
static ofstream boot_result;}

namespace NS_GP_PERF {
int JobMin=-1,JobMax=-1;
static int isample,nb_sample;
static double compteur=0.;
	string gp_fileRoot;
}


int skipln(ifstream& st) { //appar pas de built-in pour ?a (cf Primer+ p.1003
  while (st.get()!= '\n') continue;
return 0;
}

//---------------- Recupere et affiche la version de Genepop -------------------
//------------------------------------------------------------------------------
void afficher_version() {
//using namespace NS_GP;
 	effacer_ecran();
    noR_cout<<"Genepop version "<< getSetting("version")<<"\n\n";
}

void ZeGenepopSound() {
    noR_cout<<"\a\a\a";
}

int set_MC_parameters(bool pasGlobTest) {
// noter que les 3 variables sont initialis?es ? 1 ? la d?claration et ?vent ensuite modif par les settings
string choix;
unsigned long int setting; // same type as batchnbr
long signcheck;
bool confirmed=false;
        if (!alwaysAskBool && !(dem<100) && !(batchnbr<10) && !(batchlgth<400)) return 0;
        effacer_ecran();
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
            printf("\nSetting Markov chain parameters\n");
        #endif

        if (dem<100) setting=10000; else setting=dem;
        if (!pauseGP) dem=setting; /*prendce qu'il a de mieux en mode batch*/
        if (alwaysAskBool) confirmed=false; else confirmed=true;
        while (!confirmed || dem<100) {
            noR_cout<<"Dememorization number (default = "<<setting<<"): ";
//            cin.getline(choix,100); //fa?on correcte de traiter les possibles return
            getline(cin,choix);
            signcheck=atol(choix.c_str()); // missing number treated as 0
            if (signcheck<0) {
              RnoR_cerr<<"negative values are not allowed"<<endl;
            } else {
              dem=size_t(signcheck);
              if (dem==0){
                dem=setting;
                noR_cout<<dem<<endl;
              }
              if (dem<100){
                RnoR_cerr<<"A value less than 100 is not allowed"<<endl;
              } else confirmed=true;
            }
        }

        if (batchnbr<10) {
           if (pasGlobTest) setting=100; else setting=20;
        } else setting=batchnbr; //test si dans settings ou non
        if (!pauseGP)  batchnbr=setting; /*si on lui donne pas de valeur valable par defaut sous perf*/
        if (alwaysAskBool) confirmed=false; else confirmed=true;
        while (!confirmed || batchnbr<10) {
            noR_cout<<"\nNumber of batches (default = "<<setting<<"): ";
//            cin.getline(choix,100); //fa?on correcte de traiter les possibles return
            getline(cin,choix);
            signcheck=atol(choix.c_str()); // missing number treated as 0
            if (signcheck<0) {
              RnoR_cerr<<"negative values are not allowed"<<endl;
            } else {
              batchnbr=size_t(signcheck);
              if (batchnbr==0){
                batchnbr=setting;
                noR_cout<<batchnbr<<endl;
              }
              if (batchnbr<10){
                RnoR_cerr<<"A value less than 10 is not allowed"<<endl;
              } else confirmed=true;
            }
        }

        if (batchlgth<400) setting=5000; else setting=batchlgth;
        if (!pauseGP)  batchlgth=setting; /*si on lui donne pas de valeur valable par defaut sous perf*/
        if (alwaysAskBool) confirmed=false; else confirmed=true;
        while (!confirmed || batchlgth<400){
            noR_cout<<"\nIterations per batch (default = "<<setting<<"): ";
//            cin.getline(choix,100); //fa?on correcte de traiter les possibles return
            getline(cin,choix);
            signcheck=atol(choix.c_str()); // missing number treated as 0
            if (signcheck<0) {
              RnoR_cerr<<"negative values are not allowed"<<endl;
            } else {
              batchlgth=size_t(signcheck);
              if (batchlgth==0){
                batchlgth=setting;
                noR_cout<<batchlgth<<endl;
              }
              if (batchlgth<400){
                RnoR_cerr<<"A value less than 400 is not allowed"<<endl;
              } else confirmed=true;
            }
        }
        effacer_ecran();
return 0;
}



//---------------- R?cup?re gp_file, Date et time... -------------------
//------------------------------------------------------------------------------
int glance_fichier_in(bool fileCompareBool) {
using namespace NS_GP;
int local_loc_nbr,
	local_pop_nbr;
  	string bidon;
    ifstream locfichier(fichierIn.c_str(),ios::in); // IF locfichier is local, no need to reset the eof bit to reopen it....
    if ( ! locfichier.is_open()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        cerr<<"(!) From glance_fichier_in(): Could not reopen "<<fichierIn<<".\n";
        cerr<<"Exiting...\n";
        if (cinGetOnError) cin.get();
        #endif

        #endif

			genepop_exit(1, "(!) From glance_fichier_in(): Could not reopen " );
		}
    locfichier>>bidon;
    if (gp_file.size()==0) gp_file=bidon;
    // sinon c'est que gp_file a ?t? donn? dans settings
    if (fileCompareBool) // alors on verif que le fichier in correspond au m?me que fichier in Settings
//piegacon: strcmp return O (=false) quand les chaines sont identiques
       if (strcmp(gp_file.c_str(),bidon.c_str())) return 0; // ils sont diff?rents (nouvelles donn?es), on ignore le fichier in
    //else ce sont les m?mes ou il n'y a pas de fichier dans settings
    //on va juste chercher les dates et times
    //(fichier.in ne semble plus servir ? rien d'autre qu'au test prec et ? cet affichage)
    locfichier>>local_pop_nbr>>local_loc_nbr;
    //saut de ligne
    getline(locfichier,bidon);
    for (int i=0; i<local_loc_nbr; i++) getline(locfichier,bidon);
    for (int i=0; i<local_pop_nbr; i++) getline(locfichier,bidon);
    locfichier>>fichDATE;
    locfichier>>fichTIME;
    locfichier.close();
return 0;
}


//---------------- Controle de la saisie dans les menus ------------------------
//------------------------------------------------------------------------------
int controle_choix() {
using namespace NS_GP;
//filtre d'abord leschoix <>0...9 connus, puis traite les choix >0 ...9< sinon retunr -1;
string saisie;
int choix;
	cin>>saisie;
	cin.ignore();
	if (strcmp(saisie.c_str(),"c")==0 || strcmp(saisie.c_str(),"C")==0) return 10;
	if (strcmp(saisie.c_str(),"a")==0) return 11;
	if (strcmp(saisie.c_str(),"e")==0) return 12;
	if (saisie.length() > 1) choix = -1;
	  else {
			choix = atoi(saisie.c_str());
	        if (choix>9) choix=-1;
	  }
return choix;
}




//---------------- nv nom fichier et verif qu'il existe -------------------
//------------------------------------------------------------------------------

int ask_new_gp_file() {
using namespace NS_GP;
    effacer_ecran();
    afficher_version();
    #ifdef COMPATIBILITYRCPP
        // Rcpp::R
    #else
    printf("\nName of the data file ? (press ENTER to quit)\n");
    #endif

        //Entr?e du fichier de donn?es ? traiter
//   	cin>> gp_file; //waits for non-empty line !
    getline(cin,gp_file);  //works with empty lines
//   	cin.ignore(); // vire le \n restant apr?s la lecture dans cin (obsolete ?)
	if (gp_file.size() == 0) {exit_genepop = true; return 0;} //ADD by jimmy
	{string::size_type pos=gp_file.find(".");
//        cin.ignore(); // vire le \n restant apr?s la lecture dans cin (obsolete ?)
		if (pos!=string::npos && strcmp(gp_file.substr(pos+1).c_str(),"txt")!=0) { // first test = did find a "."
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        cerr<<"Genepop input file name shouldn't have an extension different from '.txt'.\n";
//cerr<<gp_file.substr(pos+1).c_str();
        cerr<<"Exiting...\n";
        if (cinGetOnError) cin.get();
        #endif

			genepop_exit(1, "Genepop input file name shouldn't have an extension different from '.txt'.");
		}
	}//block pos
return 0;
}


//---------------- sanitizes, parses, new fichier in et menu -------------------
//------------------------------------------------------------------------------

int check_gp_file_menu(bool verbose) {
using namespace NS_GP;
  noR_cout<<"Current input file: "<<gp_file<<endl;
//    cout<<"Last read at date: "<<fichDATE<<", time: "<<fichTIME<<"\n";
	set_eof_check_EOLtype(gp_file,true); //sanitizes end of file
	fichier_genepop->parseFile();
	if (verbose) {
		fichier_genepop->affiche_nb_alleles();
		if (pauseGP) {
            noR_cout<<"(Return) to continue"<<endl; getchar();
        } else {
        }
	}
	fichier_genepop->createFichierIN();
	glance_fichier_in(false);
  if(exit_genepop) {return 0;}  //ADD by jimmy
	menu();

return(0);
}

int HWexact() {
int choix;
    while (true) {
          if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("Hardy Weinberg tests:\n");
             printf("\n");
             printf("\n");
             printf("HW test for each locus in each population:\n");
             printf("   H1 = Heterozygote deficiency.......1\n");
             printf("   H1 = Heterozygote excess...........2\n");
             printf("   Probability test...................3\n");
             printf("\n");
             printf("Global test:\n");
             printf("   H1 = Heterozygote deficiency.......4\n");
             printf("   H1 = Heterozygote excess...........5\n");
             printf("\n");
             printf("Main menu.............................6\n");
             #endif

             if (((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1)))
                choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();
             switch (choix) {
                    case 1 : HWtest(1);return(0);
                    case 2 : HWtest(2);return(0);
                    case 3 : HWtest(3);return(0);
                    case 4 : HWtest(4);return(0);
                    case 5 : HWtest(5);return(0);
                    case 6 : return 0;
             }
    }
}

int HWfileMenu() {
int choix;

    while (true) {
        if(exit_genepop) {return 0;}  //ADD by jimmy
        effacer_ecran();
        afficher_version();
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
        printf("Hardy Weinberg tests:\n");
        printf("\n");
        printf("\n");
        printf("HW test for each locus in each population:\n");
        printf("   H1 = Heterozygote deficiency .................1\n");
        printf("   H1 = Heterozygote excess .....................2\n");
        printf("   Probability test .............................3\n");
        printf("\n");
        printf("Allele frequencies, expected genotypes, Fis .... 4\n");
        printf("Quit ........................................... 5\n");
        #endif


        if (HWfileOptions.size()>boucle) {
            choix=HWfileOptions[boucle]; //pas encore increment?e
            boucle++;
            if (choix==1) {
                noR_cout<<"\n (1) Heterozygote deficiency, chosen from settings file\n";
            }
            else if (choix==2) {
                noR_cout<<"\n (2) Heterozygote excess, chosen from settings file\n";
            }
            else if (choix==3) {
                noR_cout<<"\n (3) Probability test, chosen from settings file\n";
            }
            else if (choix==4) {
                noR_cout<<"\n (4) Basic information, chosen from settings file\n";
            }
            else if (choix>5) {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                    cerr<<"\nIncorrect choice '"<<choix<<"' given in settings file. I exit.";if (cinGetOnError) cin.get();
                #endif

              genepop_exit(-1, "Incorrect choice given in settings file.");
            }
            if (pauseGP) {
                noR_cout<<"(Return) to continue"<<endl; getchar();
            }
        } else if (!pauseGP) {exit_genepop = true; return 0;} //ADD by jimmy // mode batch et fin de HWfileOptions
        else { // fin de HWfileOptions en mode non batch
           HWfileOptions.resize(0); // vidage par pr?caution peut etre superflue
           choix =  controle_choix();
        }

        switch (choix) {
        case 1 : HWtest(1);HWfileMenu(); return(0);
            case 2 : HWtest(2);HWfileMenu(); return(0);
            case 3 : HWtest(3);HWfileMenu(); return(0);
            case 4 : HWfile_info();HWfileMenu(); return(0);
            case 5 : exit_genepop = true; return 0; //ADD by jimmy
        }
    }
}


int LDexact() {
int choix;
    while (true) {
          if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("\n");
             printf("\n");
             printf("Pairwise associations (haploid and genotypic disequilibrium):\n");
             printf("      Test for each pair of loci in each population ......... 1\n");
             printf("      Only create genotypic contingency tables .............. 2\n");
             printf("\n");
             printf("Menu  ....................................................... 3\n");
             #endif


             if ((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1))
               choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();

             switch (choix) {
                    case 1 : LDtest();return(0);
                    case 2 : LDtables();return(0);
                    case 3 : return 0;
             }
    }
}

int Diffexact() {
int choix;
    while (true) {
            if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("\n");
             printf("\n");
             printf(" Testing population differentiation :\n");
             printf("\n");
             printf("      Genic differentiation:\n");
             printf("           for all populations ........................ 1\n");
             printf("           for all pairs of populations ............... 2\n");
             printf("\n");
             printf("      Genotypic differentiation:\n");
             printf("           for all populations ........................ 3\n");
             printf("           for all pairs of populations ............... 4\n");
             printf("\n");
             printf("      Main menu  ...................................... 5\n");
             #endif


             if ((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1))
               choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();
             switch (choix) {
                    case 1 : Difftest(1);return(0);
                    case 2 : Difftest(2);return(0);
                    case 3 : Difftest(3);return(0);
                    case 4 : Difftest(4);return(0);
                    case 5 : return 0;
             }
    }
}



//---------------- option 1.1 -- 1.5 -------------------
int HWtest(int statindic) {
//statindic=1--5: options du menu HWexact;
using namespace NS_GP;
bool defic=false;
bool prob=false;
bool glob=false;
vector<vector<bool> > indic;
vector<string>markName;
vector<string>exacName;
vector<int> effallnbr;
vector<double> Uloc;
vector<double> Upop;
string hw_outfile;
if (statindic==1 || statindic==4) defic=true;
    if (statindic==3) prob=true;
    if (statindic>3) glob=true;
    if (glob && HWfileBool) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::R
        #else
            cerr<<"(!) No global test on HWFile; check MenuOptions or HWFile settings.\nExiting.";
            if (cinGetOnError) cin.get();
        #endif

        genepop_exit(-1, "(!) No global test on HWFile; check MenuOptions or HWFile settings.");
    }
    hardy1(defic,prob,glob,HWfileBool,hw_outfile); //point d'entree des bool+[create P<>_L<> files si !HWfileBool] // new proba SI !HWfileBool
    if (!HWfileBool && !glob) ecriture_sample_HW(hw_outfile); //from hardy2 with hardly any change; ecriture bidon
    lecture_fich_PL(true,effallnbr); //from hardy2 with hardly any change; lit les fichiers pour d?terminer ceux analysables et comment
    traitement_result_fichiers(markName,exacName,effallnbr,hw_outfile); //from hardy2 with hardly any change; contient appel ? set_MC_parameters()
// et remplit markName etc
if (!glob) enum_test_et_affich(exacName);  //enumeration tests plus affichages divers (~ Fis_D/E/P) et vide exacName
    else global_U_initialize(indic,Uloc,Upop);

    HW_Pvalues_chains(markName); //MCMC tests (~ Wein_D/E/P ou Chain_U)... et vide markName (sinon s'accumule ? chaque appel HWtest)

    if (!HWfileBool) {
       if (glob) HW_Pvalues_compile(indic,Uloc,Upop,hw_outfile); //analyse fichiers des chaines par locus pour test globaux
       else {
            noR_cout << "\n\n...I'm building the output file...";
            fic_lect(); // repasse sur tous lesfichiers pour voir les non info etc
            ecriture_result(hw_outfile); // appelle analyse_pop => messages finaux; renomme fichier sortie
       }
       delete_proba(); // de ptr alloues dans hardy1
       Genclean_HW(); // removes all P<pop>_L<loc> files
    } else {
       //28/04/09: if HWfile then previously  [if (!glob) TRUE] enum_test_et_affich(); => message already displayed
       // maybe the following line is always useless?
       if (glob) {
           noR_cout<<"Edit the file "<<hw_file<<" for results";
          }
       if (!perf) ZeGenepopSound();
       if (pauseGP) {
           noR_cout<<endl<<"(Return) to continue"<<endl; getchar();
           }
    }
return 0;
}

int print_p(double pchi, ostream& fichier_out, int prec, bool endline) {
  streamsize old_prec; 
  if (pchi<1e-04) { 
    old_prec = fichier_out.precision(2);
    fichier_out<<scientific<<pchi<<fixed;
    fichier_out.precision(old_prec);
  } else {
    old_prec = fichier_out.precision(prec);
    fichier_out<<fixed<<pchi;
  }
  if (endline) fichier_out<<endl;
  fichier_out.precision(old_prec);
  return(0);
}

//---------------- option 2.1, 2.3 -------------------
int LDtest() {
using namespace NS_GP;
string nom_fic;
vector<double>testresult(3);
CGenobilocus donnees;
size_t nb_sam=fichier_genepop->pops.size();
size_t nb_locus=fichier_genepop->loci.size();
long int totalNbrTests=nb_sam*nb_locus*(nb_locus-1)/2;
vector<vector<double> > tabF;
size_t pairit;
int ddl=0,ntests=0;
double fisher,liptak;
float pchi=0.0;
bool infini;
time_t _StartTime,_EndTime;
float _secs;
bool affichIntermBool=false;
switchFnPtr=&Cctable::switchSP; // pas de statistique G allelique

    time(&_StartTime);
    tabF.resize(nb_locus*(nb_locus-1)/2);
    for(vector<vector<double> >::iterator i = tabF.begin(); i < tabF.end(); i++) i->resize(nb_sam,0);
    noR_cout<<"Genotypic linkage disequilibrium between each pair of loci"<<endl;
    set_MC_parameters(true);
    nom_fic=gp_file+".DIS";
    ofstream wdisOut(nom_fic.c_str());
    if(!wdisOut.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    wdisOut.setf(ios::left,ios::adjustfield);
    wdisOut<<setprecision(6);
    wdisOut<<"Genepop "<<getSetting("version")<<", Genotypic linkage disequilibrium"<<endl<<endl;
    wdisOut<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    wdisOut<<endl;
    wdisOut<<"Number of populations detected : "<<nb_sam<<endl;
    wdisOut<<"Number of loci detected        : "<<nb_locus<<endl;
    if (nb_locus<2) {
       wdisOut<<"\nOnly "<<nb_locus<<" locus; no test.";
       goto lafin;
    }
    // ELSE
    wdisOut<<endl;
    wdisOut<<"Markov chain parameters"<<endl;
    wdisOut<<"\tDememorisation       : "<<dem<<endl;
    wdisOut<<"\tBatches              : "<<batchnbr<<endl;
    wdisOut<<"\tIterations per batch : "<<batchlgth<<endl;
    wdisOut<<endl;
    wdisOut<<setw(16)<<"Pop"<<setw(9)<<"Locus#1"<<setw(11)<<"Locus#2"<<setw(13)<<"P-Value"<<"S.E.     Switches"<<endl;
    wdisOut<<setw(16)<<"----------"<<setw(9)<<"-------"<<setw(11)<<"-------"<<setw(13)<<"--------"<<"-------- --------"<<endl;

    for (size_t ii=0;ii<nb_sam;ii++) {
        pairit=0;
        for (size_t jj=0;jj<nb_locus;jj++) {
            for (size_t kk=0;kk<jj;kk++) {
                wdisOut<<setw(16)<<fichier_genepop->pops[ii]->popName().substr(0,14)<<setw(9)
                        <<fichier_genepop->loci[kk]->locusName.substr(0,8)<<setw(11)<<fichier_genepop->loci[jj]->locusName.substr(0,8)<<setw(13);
                donnees=fichier_genepop->read_bilocus(ii,kk,jj); //map de (CGenotypes= map+...)
                if (donnees.getMinDim()<1) {
                       wdisOut<<"No data"<<endl;
                       tabF[pairit][ii]=-1;
                } else if (donnees.getMinDim()<2) {
                       wdisOut<<"No contingency table"<<endl;
                       tabF[pairit][ii]=-1;
                } else {
                  Cctable ctable(donnees.tabule()); //Classe pour test CT: cf CT_tests.cpp
                  if (!ctable.verifInfo()) {
                       wdisOut<<"No information"<<endl;
                       tabF[pairit][ii]=-1;
                  } else {
                      if (LDprobaTestBool) testresult=ctable.Proba_test(); else testresult=ctable.G_test();
       		          wdisOut<<setw(13)<<testresult[0]<<setw(8)<<fixed<<testresult[1];
                      wdisOut<<" "<<setw(8)<<right<<int(testresult[2]+0.5)<<left<<endl;
       		          tabF[pairit][ii]=testresult[0];
                  } //verifInfo
                } //getMinDim
                pairit++;
                time(&_EndTime);
                _secs=float(_EndTime- _StartTime);
                if (_secs>1) {
                   affichIntermBool=true;
                   _gotoxy(0,19);
                   noR_cout<<"Already "<<pairit<<" tests done out of "<<totalNbrTests<<"  "<<endl;
                   _StartTime=_EndTime;
                }
            }
        }
    }
    if (affichIntermBool) {
       _gotoxy(0,19);
       noR_cout<<"Already "<<totalNbrTests<<" tests done out of "<<totalNbrTests<<"  "<<endl;
    }
    if (nb_sam>1) {
        wdisOut<<"\nP-value for each locus pair across all populations"<<endl;
        wdisOut<<"(Fisher's method)"<<endl;
        for(int i = 0; i < 53; i++){
         wdisOut<<"-";
        }
        wdisOut<<endl;
        wdisOut<<setw(30)<<"Locus pair"<<setw(10)<<"Chi2"<<setw(5)<<"df"<<"P-Value"<<endl;
        wdisOut<<setw(30)<<"--------------------"<<setw(10)<<"--------"<<setw(5)<<"---"<<"--------"<<endl;
        pairit = 0;
        for(size_t ll = 0; ll < nb_locus; ll++){
            for(size_t mm = 0; mm < ll; mm++){
            fisher = 0;
            ddl = 0;
            infini = false;
            for(size_t p = 0; p < nb_sam; p++){
               if(tabF[pairit][p] == -1){
               //next p
               }else{
                  if(tabF[pairit][p] < numeric_limits<float>::epsilon()){
                    infini = true;
                    fisher-=2.*log(numeric_limits<float>::epsilon());
                  }else{
                     fisher-=2.*log(tabF[pairit][p]);
                  }
                  ddl+=2;
               }
            } // next p
            wdisOut<<setw(13)<<fichier_genepop->loci[mm]->locusName.substr(0,13)<<" & "<<setw(14)<<fichier_genepop->loci[ll]->locusName.substr(0,13);
            if(ddl > 0){
              if (infini) {
                wdisOut<<">"<<setw(9)<<fisher;
              } else {
                wdisOut<<setw(10)<<fisher;
              }
              chi2(pchi,ddl,float(fisher));
              wdisOut<<setw(5)<<ddl;
              if(pchi != -1){
                if (infini) {
                  wdisOut<<"<"<< setw(7);
                } else {
                  wdisOut<<setw(8);
                }
                print_p(pchi, wdisOut, 6, true);
              } else{
                wdisOut<<"Highly sign."<<endl;
              } // if infini et pchi
            }else{
               wdisOut<<"Not possible"<<endl;
            } // if nu
            pairit+=1; // kk = (ll-1)*ll/2+mm-ll+1
            } // next mm
        } // next ll
    } // if nb_sam>1

// DEBUT CODAGE Z TEST ; NDTRI EST L INVERSE DE NDTR -> QUANTILE
if (liptakBool) // false pour l'instant, mais garder le code
    if (nb_sam>1) {
        wdisOut<<"\nP-value for each locus pair across all populations"<<endl;
        wdisOut<<"(Liptak's method = Cochran's Z test)"<<endl;
        for(int i = 0; i < 53; i++){
         wdisOut<<"-";
        }
        wdisOut<<endl;
        wdisOut<<setw(30)<<"Locus pair"<<setw(10)<<"Z"<<setw(5)<<"#"<<"P-value"<<endl;
        wdisOut<<setw(30)<<"--------------------"<<setw(10)<<"--------"<<setw(5)<<"---"<<"--------"<<endl;
        pairit = 0;
        for(size_t ll = 0; ll < nb_locus; ll++){
            for(size_t mm = 0; mm < ll; mm++){
            liptak = 0;
            ntests=0;
            infini = false;
            for(size_t p = 0; p < nb_sam; p++){
               if(tabF[pairit][p] == -1){
               //next p
               }else{
                 ntests+=1;
                 if(tabF[pairit][p] < numeric_limits<float>::epsilon()){
                     noR_cout<<"\nQuick patch in Z test for estimated p=0\n";
                     liptak+=ndtri(numeric_limits<float>::epsilon()); //QUICK PATCH quantile at (epsilon)~est.pvalue
                     infini = true;
                     // next p
                  }else if(tabF[pairit][p] == 1){
                    noR_cout<<"\nQuick patch in Z test for estimated p=1\n"; // combining the two patches is not perfect...
                    liptak+=ndtri(1.-1./(batchnbr*batchlgth)); //QUICK PATCH
                  }else{
                     liptak+=ndtri(tabF[pairit][p]);
                  }
               }
            } // next p
            wdisOut<<setw(13)<<fichier_genepop->loci[mm]->locusName.substr(0,13)<<" & "<<setw(14)<<fichier_genepop->loci[ll]->locusName.substr(0,13);
            if(ntests > 0){
              if (infini) {
                wdisOut<<"<"<<setw(9)<<liptak;
              } else {
                wdisOut<<setw(10)<<liptak;
              }
              liptak=ndtr(liptak/sqrt(ntests));
              wdisOut<<setw(5)<<ntests;
              if(liptak != 0.0){
                if(infini){
                  wdisOut<<"<"<<setw(7);
                }else{
                  wdisOut<<setw(8);
                }
                print_p(liptak, wdisOut, 6, true);
              } else{
                wdisOut<<"Highly sign."<<endl;
              } // if infini et pchi
            }else{
               wdisOut<<"Not possible"<<endl;
            } // if nu
            pairit+=1; // kk = (ll-1)*ll/2+mm-ll+1
            } // next mm
        } // next ll
    } // if nb_sam>1
// fin if liptakBool
lafin:
    wdisOut<<"\nNormal ending."<<endl;
    wdisOut.close();
    noR_cout<<"Normal ending."<<endl;
    noR_cout<<"Edit the file "<<gp_file.c_str()<<".DIS for results"<<endl;
        if (!perf) ZeGenepopSound();
        if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}

//---------------- option 2.2: only write tables -------------------
int LDtables() {
using namespace NS_GP;
    string nom_fic;
    size_t nb_sam=fichier_genepop->pops.size();
    size_t nb_locus=fichier_genepop->loci.size();
    CGenobilocus donnees;
    CGenotypes ligne;
    time_t begTime,endTime;
    unsigned champ,gnourf,llong;
    ssize_t typ1,typ2;
    int fac1,fac2;
    size_t rowsum,maxCellCount;
    stringstream stst;
    long int compteur=0;
    size_t total=(nb_sam)*nb_locus*(nb_locus-1)/2;
    nom_fic=gp_file+".TAB";
    ofstream wdisOut(nom_fic.c_str());
    if(!wdisOut.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    wdisOut<<"Genepop "<<getSetting("version")<<", Contingency tables for genotypic disequilibrium\n\n";
    wdisOut<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")\n\n";
    wdisOut<<"Number of populations detected : "<<nb_sam<<endl;
    wdisOut<<"Number of loci detected        : "<<nb_locus<<endl;
    wdisOut<<endl;
    if (nb_locus<2) {
       wdisOut<<"\nOnly "<<nb_locus<<" locus; no table.";
       goto lafin;
    }
    time(&begTime);
    for (size_t ii=0;ii<nb_sam;ii++) {
        for (size_t jj=0;jj<nb_locus;jj++) {
            for (size_t kk=0;kk<jj;kk++) {
                donnees=fichier_genepop->read_bilocus(ii,kk,jj); //map de (CGenotypes= map+...)
                if (donnees.coding1==3 || donnees.coding1==6) fac1=1000; else fac1=100; //kk
                if (donnees.coding2==3 || donnees.coding2==6) fac2=1000; else fac2=100; //jj
                wdisOut<<"-----next table------------------\n";
                wdisOut<<"Table with "<<donnees.getligNbr()<<" lines and "<<donnees.getcolNbr()<<" columns\n";
                wdisOut<<"Pop: "<<fichier_genepop->pops[ii]->popName();
                wdisOut<<", loci: "<<fichier_genepop->loci[jj]->locusName<<" and "<<fichier_genepop->loci[kk]->locusName;
                if (donnees.getMinDim()<2) {
                    wdisOut<<" --> Missing data\n";
               } else {
                    {Cctable unetable(donnees.tabule());
                     maxCellCount=unetable.maxCellCount();
                     champ=0;
                    }
                    wdisOut<<"\n\n";
                    wdisOut<<"             "<<fichier_genepop->loci[jj]->locusName<<"\n";
                    wdisOut<<setw(13)<<" ";
                    donnees.marginal.resetIterator();
                    if (donnees.coding2>4) { //diploide 3-digits
                    // premier passage pour determiner la largeur des champs
                        donnees.marginal.resetIterator();
                        while ((typ2=donnees.marginal.getNext()) >= 0) {
                            stst<<typ2/fac2;
                            gnourf=unsigned(stst.str().size());
                            stst.str("");
                            stst.clear();
                            if (gnourf>champ) champ=gnourf;
                        }
                        donnees.marginal.resetIterator();
                        while ((typ2=donnees.marginal.getNext()) >= 0) {
                            stst<<typ2%fac2;
                            gnourf=unsigned(stst.str().size());
                            stst.str("");
                            stst.clear();
                            if (gnourf>champ) champ=gnourf;
                        }
                        gnourf=champ; //swap tempo
                        stst<<maxCellCount;
                        champ=unsigned(stst.str().size()); // int() bc setw(int n)
                        stst.str("");
                        stst.clear();
                        if (gnourf>2) {// real 3-digits : on ?crit les genos sur 2 lignes
                           champ=std::max(champ+1,gnourf+2); //espaces justifi?s par l'esth?tique probable
                           //impression premi?re ligne des genotypes du 2e locus
                           donnees.marginal.resetIterator();
                           while ((typ2=donnees.marginal.getNext()) >= 0) {
                              wdisOut<<setw(int(champ))<<left<<typ2/fac2;
                           }
                           wdisOut<<endl;
                           //impression deuxieme ligne des genotypes du 2e locus
                           wdisOut<<setw(13)<<" ";
                           donnees.marginal.resetIterator();
                           while ((typ2=donnees.marginal.getNext()) >= 0)
                              wdisOut<<setw(int(champ))<<left<<typ2%fac2;
                        } else { // on peut compresser en moins de 3 digits
                            //on recalcule la largeur
                            //champ est la largeur de maxCellCount ? ce point
                            donnees.marginal.resetIterator();
                            while ((typ2=donnees.marginal.getNext()) >= 0) {
                                stst<<typ2/fac2<<"."<<typ2%fac2<<" ";
                                gnourf=unsigned(stst.str().size());
                                stst.str("");
                                stst.clear();
                                if (gnourf>champ) champ=gnourf;
                            }
                            champ+=1; //espaces justifi?s par l'esth?tique probable
                            // impression genotypes du 2e locus
                            donnees.marginal.resetIterator();
                            while ((typ2=donnees.marginal.getNext()) >= 0) {
                                    stst<<typ2/fac2<<"."<<typ2%fac2;
                                    wdisOut<<setw(int(champ))<<left<<stst.str();
                                    stst.str("");
                                    stst.clear();
                            }
                        }
                    } else if (donnees.coding2==4) { //diploide 2-digits
                        //on calcule la largeur
                        donnees.marginal.resetIterator();
                        while ((typ2=donnees.marginal.getNext()) >= 0) {
                            stst<<typ2/fac2<<"."<<typ2%fac2;
                            gnourf=unsigned(stst.str().size());
                            stst.str("");
                            stst.clear();
                            if (gnourf>champ) champ=gnourf;
                        }
                        gnourf=champ;
                        stst<<maxCellCount;
                        champ=unsigned(stst.str().size());
                        stst.str("");
                        stst.clear();
                        champ=std::max(champ+1,gnourf+2);
                        // impression genotypes du 2e locus
                        donnees.marginal.resetIterator();
                        while ((typ2=donnees.marginal.getNext()) >= 0) {
                                stst<<typ2/fac2<<"."<<typ2%fac2;
                                wdisOut<<setw(int(champ))<<left<<stst.str();
                                stst.str("");
                                stst.clear();
                        }
                    } else { //haploid
                        donnees.marginal.resetIterator();
                        while ((typ2=donnees.marginal.getNext()) >= 0) {
                            stst<<typ2;
                            gnourf=unsigned(stst.str().size());
                            stst.str("");
                            stst.clear();
                            if (gnourf>champ) champ=gnourf;
                        }
                        gnourf=champ;
                        stst<<maxCellCount;
                        champ=unsigned(stst.str().size());
                        stst.str("");
                        stst.clear();
                        champ=std::max(champ+1,gnourf+2);
                        donnees.marginal.resetIterator();
                        while ((typ2=donnees.marginal.getNext()) >= 0)
                            wdisOut<<setw(int(champ))<<left<<typ2;
                    }
                    wdisOut<<endl;
                    if ((champ * donnees.marginal.getNumber()) < 239) llong = unsigned(champ * donnees.marginal.getNumber()); else llong = 239;
                    wdisOut<<setw(13)<<left<<fichier_genepop->loci[kk]->locusName.substr(0,11);
                    for (unsigned k=0;k<(llong);k++) wdisOut<<"_";
                    wdisOut<<endl;
                    donnees.resetIterator();
                    while ((typ1=donnees.getNext()) >= 0) {
                        if (fac1>99) wdisOut<<setw(4)<<right<<typ1/fac1<<"."<<setw(8)<<left<<typ1%fac1;
                        else wdisOut<<setw(4)<<right<<typ1<<setw(9)<<" "<<left;
                        donnees.marginal.resetIterator();
                        ligne=donnees.mapmap[size_t(typ1)];
                        rowsum=0;
                        while ((typ2=donnees.marginal.getNext()) >= 0) {
                            wdisOut<<setw(int(champ))<<ligne.getEffective(typ2);
                            rowsum+=ligne.getEffective(typ2);
                        }
                        wdisOut<<"  "<<setw(int(champ))<<rowsum<<endl;
                    }
                    wdisOut<<setw(13)<<" ";
                    for (unsigned k=0;k<(llong);k++) wdisOut<<"_";
                    wdisOut<<endl<<setw(13)<<" ";
                    donnees.marginal.resetIterator();
                    while ((typ2=donnees.marginal.getNext()) >= 0)
                        wdisOut<<setw(int(champ))<<donnees.marginal.getEffective(typ2);
                    wdisOut<<"  "<<setw(int(champ))<<donnees.marginal.getSum()<<endl;
               }
               compteur++;
        	   time(&endTime);
               if (float(endTime- begTime) >1) {
                  _gotoxy(0,12);
                  noR_cout<<"Already "<<compteur<<" tables written out of "<<total<<"  ";
                  begTime=endTime;
               }
            }  //kk
        }   //jj
    }   //ii
    wdisOut<<"-----End of tables---------------";
    _gotoxy(0,12);
    noR_cout<<"Already "<<compteur<<" tables written out of "<<total<<"  ";
lafin:
    wdisOut<<"\nNormal ending."<<endl;
    wdisOut.close();
    noR_cout<<"\nNormal ending."<<endl;
    noR_cout<<"Edit the file "<<gp_file.c_str()<<".TAB for results"<<endl;
    if (!perf) ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}


//---------------- option 3 -------------------
int Difftest(int statindic) {
using namespace NS_GP;
string nom_fic;
vector<double>testresult(3);
CGenobilocus donnees;
size_t nb_sam=fichier_genepop->pops.size();
size_t nb_locus=fichier_genepop->loci.size();
vector<vector<double> > tabF;
int ddl=0,ntests=0;
double fisher,liptak;
float pchi=0.0;
bool infini;
string statstring;
time_t _StartTime,_EndTime;
float _secs;
bool affichIntermBool=false;

    time(&_StartTime);
    set_MC_parameters(true);
    nom_fic=gp_file;
    if (statindic==1) {
       tabF.resize(1);
       statstring="Genic differentiation";
       if (genicProbaTestBool) nom_fic+=".PR";
       else nom_fic+=".GE";
    } else if (statindic==2) {
       tabF.resize(nb_sam*(nb_sam-1)/2);
       statstring="Genic differentiation for each population pair";
       if (genicProbaTestBool) nom_fic+=".PR2";
       else nom_fic+=".GE2";
    } else if (statindic==3) {
       tabF.resize(1);
       statstring="Genotypic differentiation";
       if (alleleNbrTestBool) nom_fic+=".AL";        //private option for EI
       else if (geneDivTestBool) nom_fic+=".GD";        //private option for EI
       else nom_fic+=".G";
       genicProbaTestBool=false;
    } else if (statindic==4) {
       tabF.resize(nb_sam*(nb_sam-1)/2);
       statstring="Genotypic differentiation for each population pair";
       nom_fic+=".2G2";
       genicProbaTestBool=false;
    }
    for(vector<vector<double> >::iterator i = tabF.begin(); i < tabF.end(); i++) i->resize(nb_locus);
    if (genicProbaTestBool) statstring+=" (Fisher's exact Probability test)";
    else if (alleleNbrTestBool) statstring+=" (test based on trend on number of alleles)";        //private option for EI
    else if (geneDivTestBool) statstring+=" (test based on trend on gene diversity)";        //private option for EI
    else statstring+=" (exact G test)";
    noR_cout<<statstring<<endl;
    ofstream fichier_out(nom_fic.c_str());
    if(!fichier_out.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    fichier_out.setf(ios::left,ios::adjustfield);
    if ((statindic%2)==0) fichier_out<<setprecision(5); // probablement ? generaliser aux autres cas
    else fichier_out<<setprecision(6);
    fichier_out<<"Genepop "<<getSetting("version")<<", "<<statstring<<endl<<endl;
    fichier_out<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    fichier_out<<endl;
    fichier_out<<"Number of populations detected : "<<nb_sam<<endl;
    fichier_out<<"Number of loci detected        : "<<nb_locus<<endl;
    if (nb_sam<2) {
       fichier_out<<"\nOnly "<<nb_sam<<" populations, no differentiation test.\n\n";
       goto lafin;
    }
    fichier_out<<endl;
    fichier_out<<"Markov chain parameters"<<endl;
    fichier_out<<"\tDememorisation       : "<<dem<<endl;
    fichier_out<<"\tBatches              : "<<batchnbr<<endl;
    fichier_out<<"\tIterations per batch : "<<batchlgth<<endl;
    fichier_out<<endl;
// BOUCLE PRINCIPALE
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
        fichier_out<<endl<<"Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
        fichier_out<<"=================================\n";
        if (fichier_genepop->coding[iLoc]<4 && statindic>2) {
           fichier_out<<"Genotypic test not applicable on haploid data"<<endl<<endl;
           for(vector<vector<double> >::iterator i = tabF.begin(); i < tabF.end(); i++) (*i)[iLoc]=-2;
        } else crunchLocTable(statindic,iLoc,fichier_out,&tabF); //version tabF double

        time(&_EndTime);
        _secs=float(_EndTime- _StartTime);
        if (_secs>1) {
           affichIntermBool=true;
           _gotoxy(0,19);
           noR_cout<<"Already "<<iLoc+1<<" loci analyzed out of "<<nb_locus<<"  "<<endl;
           _StartTime=_EndTime;
        }
    }
    if (affichIntermBool) {
       _gotoxy(0,19);
       noR_cout<<"Already "<<nb_locus<<" loci analyzed out of "<<nb_locus<<"  "<<endl;

    }
//TESTS MULTIPLES
    if ((statindic%2)==0) { //sur paires
      size_t pairit;
        if (nb_sam>1) {
            fichier_out<<"\n=================================\n";
            fichier_out<<"\nP-value for each population pair across all loci"<<endl;
            fichier_out<<"(Fisher's method)"<<endl;
            for(int i = 0; i < 53; i++){
             fichier_out<<"-";
            }
            fichier_out<<endl;
            fichier_out<<setw(30)<<"Population pair"<<setw(11)<<"Chi2"<<setw(6)<<"df"<<"P-Value"<<endl;
            fichier_out<<setw(30)<<"--------------------"<<setw(11)<<"---------"<<setw(6)<<"-----"<<"---------"<<endl;
            pairit = 0;
            for(size_t ll = 0; ll < nb_sam; ll++){
                for(size_t mm = 0; mm < ll; mm++){
                fisher = 0;
                ddl = 0;
                infini = false;
                for(size_t p = 0; p < nb_locus; p++){
                   if(tabF[pairit][p] <= -1){ // -1: crunchloctable indic donnees pas testables; -2: test genotypic appele sur donn hap
                   //next p
                   }else{
                      if(tabF[pairit][p] == 0){
                         infini = true;
                         fisher+=2.*log(batchnbr*batchlgth); // i.e. *-=* 2 log(bound on freq)
                         ddl+=2;
                         // next p
                      }else{
                         ddl+=2;
                         fisher-=2.*log(tabF[pairit][p]);
                      }
                   }
                } // next p
                fichier_out<<setw(13)<<fichier_genepop->pops[mm]->popName().substr(0,13)<<" & "<<setw(14)<<fichier_genepop->pops[ll]->popName().substr(0,13);
                if(ddl > 0){
                  if (infini) {
                    fichier_out<<">"<<setw(9)<<fisher;
                  } else {
                    fichier_out<<setw(10)<<fisher;
                  }
                  chi2(pchi,ddl,float(fisher));
                  fichier_out<<" "<<setw(6)<<ddl;
                  if(pchi != -1){
                    if (infini) {
                      fichier_out<<"<"<< setw(7);
                    } else {
                      fichier_out<<setw(8);
                    }
                    print_p(pchi, fichier_out, 6, true);
                  } else{
                    fichier_out<<"Highly sign."<<endl;
                  }
                }else{
                   fichier_out<<"Not possible"<<endl;
                } // if nu
                pairit+=1; // kk = (ll-1)*ll/2+mm-ll+1
                } // next mm
            } // next ll
        } // if nb_sam>1
    if (liptakBool) // la aussi false...
        if (nb_sam>1) {
            fichier_out<<"\nP-value for each population pair across all loci"<<endl;
            fichier_out<<"(Liptak's method = Cochran's Z test)"<<endl;
            for(int i = 0; i < 53; i++){
             fichier_out<<"-";
            }
            fichier_out<<endl;
            fichier_out<<setw(30)<<"Population pair"<<setw(11)<<"Z"<<setw(6)<<"#"<<"P-value"<<endl;
            fichier_out<<setw(30)<<"--------------------"<<setw(11)<<"--------"<<setw(6)<<"-----"<<"--------"<<endl;
            pairit = 0;
            for(size_t ll = 0; ll < nb_sam; ll++){
                for(size_t mm = 0; mm < ll; mm++){
                liptak = 0;
                ntests=0;
                infini = false;
                for(size_t p = 0; p < nb_locus; p++){
                   if(tabF[pairit][p] <= -1){
                   //next p
                   }else{
                     ntests+=1;
                     if(tabF[pairit][p] == 0){
                        noR_cout<<"\nQuick patch in Z test for estimated p=0\n";
                        liptak+=ndtri(numeric_limits<float>::epsilon()); //QUICK PATCH quantile at (epsilon)~est.pvalue
                        infini = true;
                         // next p
                      }else if(tabF[pairit][p] == 1){
                         noR_cout<<"\nQuick patch in Z test for estimated p=1\n"; // combining the two patches is not perfect...
                         liptak+=ndtri(1.-1./(batchnbr*batchlgth)); //QUICK PATCH quantile at (1-epsilon)~est.pvalue
                         // next p
                      }else{
                         liptak+=ndtri(tabF[pairit][p]); // quantile at p-value
                      }
                   }
                } // next p
                fichier_out<<setw(13)<<fichier_genepop->pops[mm]->popName().substr(0,13)<<" & "<<setw(14)<<fichier_genepop->pops[ll]->popName().substr(0,13);
                if(ntests > 0){
                  if (infini) {
                    fichier_out<<"<"<<setw(9)<<liptak;
                  } else {
                      fichier_out<<setw(10)<<liptak;
                   }
                  liptak=ndtr(liptak/sqrt(ntests));
                  fichier_out<<" "<<setw(6)<<ntests;
                  if(liptak != 0.0){
                    if(infini){
                      fichier_out<<"<"<<setw(7);
                    } else{
                      fichier_out<<setw(8);
                    }
                    print_p(liptak, fichier_out, 6, true);
                   } else {
                      fichier_out<<"Highly sign."<<endl;
                   } // if infini et pchi
                }else{
                   fichier_out<<"Not possible"<<endl;
                } // if nu
                pairit+=1; // kk = (ll-1)*ll/2+mm-ll+1
                } // next mm
            } // next ll
        } // fi nb_sam()>1
    // fin liptakBool
    } else { // pas sur paires
        if (nb_sam>1) {
            fichier_out<<"\n=================================\n";
            fichier_out<<"\nP-value across all loci"<<endl;
            fichier_out<<"(Fisher's method)"<<endl;
            for(int i = 0; i < 53; i++){
             fichier_out<<"-";
            }
            fichier_out<<endl;
            fichier_out<<setw(15)<<"Locus"<<setw(9)<<"P-Value"<<endl;
            fichier_out<<setw(15)<<"-------------"<<setw(9)<<"--------"<<endl;
            for (size_t ii=0;ii<nb_locus;ii++)
               if (tabF[0][ii]>=0) fichier_out<<setw(15)<<fichier_genepop->loci[ii]->locusName.substr(0,13)<<setw(9)<<tabF[0][ii]<<endl;
            fichier_out<<endl;
            fisher = 0;
            ddl = 0;
            infini = false;
            for(size_t p = 0; p < nb_locus; p++){
               if(tabF[0][p] <= -1){
               //next p
               }else{
                 ddl+=2;
                 if(tabF[0][p] < numeric_limits<float>::epsilon()){
                   fisher-=2.*log(numeric_limits<float>::epsilon());
                   infini = true;
                  }else{
                     fisher-=2.*log(tabF[0][p]);
                  }
               }
            } // next p
            if(ddl > 0){
              if (infini) {
                fichier_out<<"All: Chi2< ";
              } else {
                fichier_out<<"All: Chi2= ";
              }
              print_p(fisher, fichier_out, 4, false);
              fichier_out<<" (df= "<<ddl<<")";
              chi2(pchi,ddl,float(fisher));
              fichier_out.precision(6);
               if(pchi != -1){
                 if(infini){
                   fichier_out<<", P-value< ";
                 } else {
                   fichier_out<<", P-value= ";
                 }
                 print_p(pchi, fichier_out, 6, true);
               } else{
                  fichier_out<<", highly significant"<<endl;
               }
            }else{
               fichier_out<<"All: not possible"<<endl;
            } // if nu
    } // fi tabF.size()>1
    if (liptakBool) // toujours false
        if (nb_sam>1) {
            fichier_out<<"\nP-value across all loci"<<endl;
            fichier_out<<"(Liptak's method = Cochran's Z test)"<<endl;
            for(int i = 0; i < 53; i++){
             fichier_out<<"-";
            }
            fichier_out<<endl;
            fichier_out<<setw(15)<<"Locus"<<setw(9)<<"P-Value"<<endl;
            fichier_out<<setw(15)<<"-------------"<<setw(9)<<"--------"<<endl;
            for (size_t ii=0;ii<nb_locus;ii++)
               if (tabF[0][ii]>=0) fichier_out<<setw(15)<<fichier_genepop->loci[ii]->locusName.substr(0,13)<<setw(9)<<tabF[0][ii]<<endl;
            fichier_out<<endl;
            liptak = 0;
            ntests=0;
            infini = false;
            for(size_t p = 0; p < nb_locus; p++){
               if(tabF[0][p] <= -1){
               //next p
               }else{
                 ntests+=1;
                 if(tabF[0][p] < numeric_limits<float>::epsilon()){
                   noR_cout<<"\nQuick patch in Z test for estimated p=0\n";
                   liptak+=ndtri(numeric_limits<float>::epsilon());
                     infini = true;
                     // next p
                  }else if(tabF[0][p] == 1){
                    noR_cout<<"\nQuick patch in Z test for estimated p=1\n";
                    liptak+=ndtri(1.-1./(batchnbr*batchlgth)); //QUICK PATCH
                     // next p
                  }else{
                     liptak+=ndtri(tabF[0][p]);
                  }
               }
            } // next p
            if(ntests > 0){
              if (infini) {
                fichier_out<<"All: Z< "<<liptak;
              } else {
                fichier_out<<"All: Z= "<<liptak;
              }
              fichier_out<<" ("<<ntests<<" tests)";
              liptak=ndtr(liptak/sqrt(ntests));
              if(liptak != 0.0){
                if(infini){
                  fichier_out<<", P-value< ";
                } else{
                  fichier_out<<", P-value= ";
                }
                print_p(liptak, fichier_out, 6, true);
              } else{
                fichier_out<<"Highly sign."<<endl;
              } // if infini et pchi
            }else{
               fichier_out<<"All: Not possible"<<endl;
            } // if nu
        } // fi nb_sam()>1
    // fin liptakBool
    }
lafin:
    fichier_out<<"\nNormal ending."<<endl;
	fichier_out.close();
    noR_cout<<"Normal ending."<<endl;
    noR_cout<<"Edit the file "<<nom_fic<<" for results"<<endl;
    if (!perf) ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}

//---------------- option 4: private alleles -------------------
int BartonS86() {
using namespace NS_GP;
    if(exit_genepop) {return 0;} //ADD by jimmy
    string nom_fic;
    vector<unsigned long int>dummyvec;
	vector<vector<unsigned long int> >toutable(0);
    int allele;
    long int nbpriv;
    double privfreq,a10 = -.489,b10 = -.951,a25 = -.576,b25 = -1.11,a50 = -.612,b50 = -1.21,nm,nm10,nm25,nm50;
    vector<double>ssize(2,0);
	size_t nb_sam=fichier_genepop->pops.size();
	size_t nb_locus=fichier_genepop->loci.size();
    nom_fic=gp_file+".PRI";
    ofstream wdisOut(nom_fic.c_str());
    if(!wdisOut.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    wdisOut<<"Genepop "<<getSetting("version")<<", Number of migrants using private alleles\n(see Barton & Slatkin, Heredity (1986),56:409-415)\n\n";
    wdisOut<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")\n\n";
    wdisOut<<"Number of populations detected : "<<nb_sam<<endl;
    wdisOut<<"Number of loci detected        : "<<nb_locus<<endl;
    wdisOut<<endl;
    privfreq=0;nbpriv=0;
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
        toutable.resize(0);
        for(vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
            dummyvec.resize(0);
            if (fichier_genepop->coding[iLoc]>3 && estimDiploidBool) {// donnes diploides mais table genique
                fichier_genepop->loci[iLoc]->resetgIterator();
                while ((allele = fichier_genepop->loci[iLoc]->getgNext()) >= 0)
                      dummyvec.push_back((*ii)->loci[iLoc]->getgEffective(allele));
            } else if (fichier_genepop->coding[iLoc]<4 && !estimDiploidBool) { // donnees haploides
                fichier_genepop->loci[iLoc]->resetIterator();
                while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0)
                      dummyvec.push_back((*ii)->loci[iLoc]->getEffective(allele));
            }
            toutable.push_back(dummyvec);
        }
        Cctable ctable(toutable); // structure pour test CT
        if (ctable.purgeZeros(false)) {
            ctable.cumul(privfreq,nbpriv,ssize);
        }
    }
    if (nbpriv>0) {
        if (estimDiploidBool) ssize[0]/=(2*ssize[1]); //diploid !
        else ssize[0]/=ssize[1];
        privfreq/=nbpriv;
        wdisOut<<"Mean sample size: "<<ssize[0]<<endl;
        wdisOut<<"Mean frequency of private alleles p(1)= "<<privfreq<<endl<<endl;
        nm = log(privfreq) / log(10); //hmmm mindless translation of BASIC code
        nm10 = pow(10.,((nm - b10) / a10));
        nm25 = pow(10.,((nm - b25) / a25));
        nm50 = pow(10.,((nm - b50) / a50));
        if (ssize[0] <= 17.5) {nm = nm10 * 10 / ssize[0];}
        else if (ssize[0] > 37.5) {nm = nm50 * 50 / ssize[0];}
        else nm = nm25 * 25 / ssize[0];
        wdisOut<<"Number of migrants for mean N=10: "<<nm10<<endl;
        wdisOut<<"Number of migrants for mean N=25: "<<nm25<<endl;
        wdisOut<<"Number of migrants for mean N=50: "<<nm50<<endl;
        wdisOut<<"Number of migrants after correction for size= "<<nm<<endl;
    } else wdisOut<<"No private alleles."<<endl;
        wdisOut<<"\nNormal ending."<<endl;
    wdisOut.close();
    noR_cout<<"\nNormal ending."<<endl;
    noR_cout<<"Edit the file "<<nom_fic<<" for results"<<endl;
    if (!perf) ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}



//---------------- option 5.1: basic info, frequencies, Fis -------------------
int basic_info() {
using namespace NS_GP;
    string nom_fic;
	size_t nb_sam=fichier_genepop->pops.size();
	size_t nb_locus=fichier_genepop->loci.size();
	vector<CPopulation *>::iterator p; // it?rateur sur les populations
	map<int,CAllele* >::iterator mit,mjt;
	CGenotypes genos; // g?notypes pour chaque population au locus courant
	int longueur;
	size_t nb_all,popit;
	unsigned long int count,tot,heterobs,homobs;
	size_t invcoding,s_alli,s_allj; // max=1000
	char coding;
	double freq,hetero,homo,sumc,sumb,U;
    vector<double>allhet;
    nom_fic=gp_file+".INF";
    ofstream wdisOut(nom_fic.c_str());
    if(!wdisOut.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    wdisOut.setf(ios_base::fixed,ios_base::floatfield);
    wdisOut.precision(4);
    wdisOut<<"Genepop "<<getSetting("version")<<", Basic data for each locus in each population\n\n";
    wdisOut<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")\n\n";
    wdisOut<<"Number of populations detected : "<<nb_sam<<endl;
    wdisOut<<"Number of loci detected        : "<<nb_locus<<endl;
    wdisOut<<endl;

    wdisOut<<"'Expected' numbers of homozygotes or heterozygotes\nare computed using Levene's correction\n\n";
    wdisOut<<"Fis: computed as in Weir & Cockerham (1984);\nalso as in Robertson & Hill (1984).\n";
    wdisOut<<"=============================================\n  Detailed analyses\n=============================================\n\n";

	// instanciation des structures de stockage des g?notypes
	// on traite pop par pop (boucle externe) ;
	popit=0;
	for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
	    for (size_t iLoc = 0; iLoc < fichier_genepop->loci.size(); iLoc ++) {
		// it?rations sur les populations pour le locus courant
		    coding=fichier_genepop->coding[iLoc];
	      if (coding==6) {
	        invcoding=1000;
	      } else invcoding=100;
	      allhet.resize(0,0);
	      allhet.resize(invcoding); //size= 100 or 1000 : inelegant but simple
	      wdisOut<<" Pop: "<<(*p)->popName()<<"   Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
    		wdisOut<<"-----------------------------------------\n";
      	    nb_all=(*p)->loci[iLoc]->alleles.size();
       		if (nb_all<1) wdisOut<<"    No data.\n";
		    else if (coding>3) {
		      genos.clear(); // nettoyage
		      genos.fillGenotypes(iLoc, *p,coding);
		      wdisOut<<"\n     Genotypic matrix:\n      ";
		      for (mit=(*p)->loci[iLoc]->alleles.begin();
             mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
		        wdisOut.width(4);wdisOut<<mit->second->_identif;
		      }
		      wdisOut<<endl<<"        ";
		      if (nb_all>62) longueur = 62 * 4; else longueur = int(4 * nb_all);
		      for (int k=0;k<longueur;k++) wdisOut<<"_";
		      wdisOut<<endl;
		      tot=0;homo=0;hetero=0;heterobs=0;homobs=0;U=0;
		      for (mit=(*p)->loci[iLoc]->alleles.begin();
             mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
		        wdisOut<<"    "<<setw(2)<<mit->second->_identif;
		        s_alli=size_t(mit->second->_identif);
		        for (mjt=(*p)->loci[iLoc]->alleles.begin();mjt!=mit;mjt++) {
		          //cout<<invcoding<<" "<<mit->second->_identif;getchar();
		          s_allj=size_t(mjt->second->_identif);
		          count=genos.getEffective(invcoding*s_alli+s_allj);
		          wdisOut.width(4);wdisOut<<count;
		          tot+=count;
		          heterobs+=count;
		        }
		        count=genos.getEffective(invcoding*s_alli+s_alli);
		        wdisOut.width(4);wdisOut<<count;
		        tot+=count;
		        homobs+=count;
		        wdisOut<<endl;
		      }
		      tot*=2;
		      if (tot>2) {
		        wdisOut<<endl<<"    Genotypes  Obs.      Expected\n";
		        for (mit=(*p)->loci[iLoc]->alleles.begin();
               mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
		          s_alli=size_t(mit->second->_identif);
		          for (mjt=(*p)->loci[iLoc]->alleles.begin();mjt!=mit;mjt++) {
		            wdisOut<<right<<setw(6)<<mit->second->_identif<<" , "<<left<<setw(7)<<mjt->second->_identif;
		            s_allj=size_t(mjt->second->_identif);
		            count=genos.getEffective(invcoding*s_alli+s_allj);
		            allhet[s_alli]+=count;
		            allhet[s_allj]+=count;
		            wdisOut<<left<<setw(4)<<count;
		            freq=double(mit->second->getEffective()) * mjt->second->getEffective() / (tot - 1);
		            wdisOut<<right<<setw(11)<<freq<<endl;
		            hetero+=freq;
		          }
		          wdisOut<<right<<setw(6)<<mit->second->_identif<<" , "<<left<<setw(7)<<mit->second->_identif;
		          count=genos.getEffective(invcoding*s_alli+s_alli);
		          wdisOut<<left<<setw(4)<<count;
		          U+=double(count)*tot/mit->second->getEffective();
		          freq=double(mit->second->getEffective()) * (mit->second->getEffective() - 1) / (2*(tot - 1));
		          wdisOut<<right<<setw(11)<<freq<<endl;
		          homo+=freq;
		        }
		        wdisOut<<endl<<endl<<"    Expected number of homozygotes  : "<<homo<<endl;
		        wdisOut<<"    Observed number of homozygotes  : "<<homobs<<endl;
		        wdisOut<<"    Expected number of heterozygotes: "<<hetero<<endl;
		        wdisOut<<"    Observed number of heterozygotes: "<<heterobs<<endl<<endl;

		        wdisOut<<endl<<endl<<"    Allele frequencies and Fis:"<<endl;
		        wdisOut<<"    -------------------------------------------------------"<<endl;
		        wdisOut<<"                                           Fis"<<endl;
		        wdisOut<<"                                           ----------------"<<endl;
		        wdisOut<<"    Allele     Sample count     Frequency   W&C      R&H"<<endl;
		        sumb = 0; sumc = 0;
		        for (mit=(*p)->loci[iLoc]->alleles.begin();
               mit!=(*p)->loci[iLoc]->alleles.end();mit++) {
		          allhet[size_t(mit->second->_identif)]/=(tot/2);
		          freq=double(mit->second->getEffective())*(tot-mit->second->getEffective())/(tot/2)-(tot-1)*allhet[size_t(mit->second->_identif)];
		          freq/=(4 * ((tot / 2) - 1));
		          sumc+= double(allhet[size_t(mit->second->_identif)])/2;
		          sumb+=freq;
		          if(freq+ allhet[size_t(mit->second->_identif)]>0) {
		            freq = freq / (freq+ double(allhet[size_t(mit->second->_identif)])/2);
		            wdisOut<<"     "<<left<<setw(3)<<mit->second->_identif<<setw(11)<<" "<<setw(6)<<left<<mit->second->getEffective();
		            wdisOut<<setw(7)<<" "<<setw(8)<<right<<double(mit->second->getEffective())/tot<<"   "<<right<<setw(7)<<freq<<endl;
		          }
		        }
		        if(sumb+sumc>0) {
		          wdisOut<<"    "<<left<<setw(4)<<"Tot"<<setw(11)<<" "<<setw(6)<<left<<tot;
		          wdisOut<<setw(7)<<" "<<setw(8)<<right<<" "<<"   "<<right<<setw(7)<<sumb/(sumb+sumc);
		          // R&H
		          U-=tot/2;
		          freq = double(tot - 1) * (1 + 2*U / tot) - (tot - (*p)->loci[iLoc]->getNumber());
		          freq/= (2. * (tot/2 - 1) * ((*p)->loci[iLoc]->getNumber()-1));
		          wdisOut<<"  "<<right<<setw(7)<<freq<<endl;
		        }
		        wdisOut<<"    -------------------------------------------------------"<<endl;
		      } //tot>2
		    } else {
                wdisOut<<"    Haploid data, no genotypic matrix.\n";
            }
        wdisOut<<endl;
		}
		popit++;
    }

//puis boucle externe sur les loc
    wdisOut.precision(3);
    wdisOut<<"Tables of allelic frequencies for each locus:\n\n";
    for (size_t iLoc = 0; iLoc < fichier_genepop->loci.size(); iLoc ++) {
        nb_all=fichier_genepop->loci[iLoc]->alleles.size();
   		wdisOut<<" Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
   		wdisOut<<"------------------\n";
   		if (nb_all==0) wdisOut<<"   No data\n";
   		else {
            if (nb_all==1) wdisOut<<"   Pop    Allele Genes\n          ------ -----\n";
            else {
               wdisOut<<"   Pop     Alleles"<<setw(int(6*(nb_all-1)-1))<<" "<<"Genes\n           ";
               wdisOut<<setw(int(6*nb_all-1))<<setfill('-')<<""<<setfill(' ')<<" -----\n";
            }
// allele indices
        	wdisOut<<"           ";
            for (mit=fichier_genepop->loci[iLoc]->alleles.begin();
           		mit!=fichier_genepop->loci[iLoc]->alleles.end();mit++) {
            				wdisOut<<left<<setw(6)<<mit->second->_identif;
            }
            wdisOut<<endl;
            for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
                wdisOut<<left<<setw(11)<<(*p)->popName().substr(0,10);
                count=0;
                for (mit=fichier_genepop->loci[iLoc]->alleles.begin();
           		    mit!=fichier_genepop->loci[iLoc]->alleles.end();mit++) {
                       count+=(*p)->loci[iLoc]->getEffective(mit->second->_identif);
                }
                for (mit=fichier_genepop->loci[iLoc]->alleles.begin();
           		    mit!=fichier_genepop->loci[iLoc]->alleles.end();mit++) {
                    if (count>0) wdisOut<<left<<setw(6)<<double((*p)->loci[iLoc]->getEffective(mit->second->_identif))/count;
                    else wdisOut<<" -    ";
                }
                wdisOut<<left<<setw(6)<<count<<endl;
            }
        }
    wdisOut<<endl;
    }
    wdisOut<<"\nNormal ending."<<endl;
    wdisOut.close();
    noR_cout<<"Normal ending."<<endl;
    noR_cout<<"Edit the file "<<nom_fic<<" for results"<<endl;
    if (!perf) ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}

//---------------- option 5.2 and 5.3 -------------------
int Fis_Div(bool identitybool,bool warnOnExit=true) {
using namespace NS_GP;
	size_t nb_sam=fichier_genepop->pops.size();
	size_t nb_locus=fichier_genepop->loci.size();
	vector<vector<unsigned long int> >tabF(nb_sam+2);  //va contenir les effectifs genos et les types alleliques (r??crite ? chaque locus)
	std::vector<double> MSgTotHez;
	vector<int>NLoc;
	vector<int>NLocHez;
	vector<int>NtotLoc;
	vector<double>SSiTotLoc;
	NLocHez.resize(nb_sam); // pour heterosigo unweighted over loci for each pop: Nombre loci avec Ntot>0
	vector<double>SSgTotLoc;
	MSgTotHez.resize(nb_sam); // idem, la somme des MSg non pond?r?e
    NLoc.resize(nb_sam);
    NtotLoc.resize(nb_sam);
    SSiTotLoc.resize(nb_sam);
    SSgTotLoc.resize(nb_sam);
    fill(NLocHez.begin(),NLocHez.end(),0);
    fill(MSgTotHez.begin(),MSgTotHez.end(),0);
    fill(NLoc.begin(),NLoc.end(),0);
    fill(NtotLoc.begin(),NtotLoc.end(),0);
    fill(SSiTotLoc.begin(),SSiTotLoc.end(),0);
    fill(SSgTotLoc.begin(),SSgTotLoc.end(),0);
    string nom_fic;
    nom_fic=gp_file;
    if (identitybool) nom_fic+=".DIV"; else nom_fic+=".MSD";
    ofstream fichier_out(nom_fic.c_str());
    if(!fichier_out.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    fichier_out.setf(ios::left,ios::adjustfield);
    if (identitybool) {
        fichier_out<<"Genepop "<<getSetting("version")<<"\nAllele frequency based correlation (Fis) and diversity ('heterozygosity')\n\n";
        fichier_out<<"One locus estimates following standard ANOVA as in Weir and Cockerham (1984)\n\n";
    } else {
        fichier_out<<"Genepop "<<getSetting("version")<<"\nAllele size based covariance (rhoIS) and diversity (Fis) and diversity\n (MSD: mean squared allele size difference)\n\n";
        fichier_out<<"One locus estimates following standard ANOVA analogous to Weir and Cockerham (1984)\n";
        fichier_out<<"See e.g. Michalakis and Excoffier (1996) for formulas\n";
    }
    fichier_out<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    fichier_out<<endl;
    fichier_out<<"Number of populations detected : "<<nb_sam<<endl;
    fichier_out<<"Number of loci detected        : "<<nb_locus<<endl;
    fichier_out<<endl;
    if (identitybool) {
       fichier_out<<fixed<<setprecision(4);
    } else {
       fichier_out<<setprecision(5); //  <<fixed mettrait 5 decim apr?s point
    }
// BOUCLE PRINCIPALE
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
        fichier_out<<endl<<"Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
        fichier_out<<"=================================\n";
        if (!identitybool) {
           if (taille.size()<=iLoc || taille[iLoc].size()==0)
              fichier_out<<" Allele sizes = allele names\n";
           else {
              fichier_out<<" Allele names: ";
              for (map<int,int>::iterator ptr=taille[iLoc].begin();ptr!=taille[iLoc].end();ptr++)
                  fichier_out<<setw(4)<<ptr->first;
              fichier_out<<"\n Allele sizes: ";
              for (map<int,int>::iterator ptr=taille[iLoc].begin();ptr!=taille[iLoc].end();ptr++)
                  fichier_out<<setw(4)<<ptr->second;
              fichier_out<<"\n(for unspecified alleles, size = name)\n";
           }
           fichier_out<<"---------------------------------\n";
        }
        if (fichier_genepop->coding[iLoc]<4) {
           fichier_out<<"Fis computation not applicable for haploid data"<<endl<<endl;
           crunchLocTable(1,iLoc,fichier_out,&tabF); //version tabF entier dans CT_tests
           tabFtotabM(&tabF); //contient des new : delete ci dessous
           if (identitybool) fichier_out<<"\n          1-Qinter      \n-------------------\n";
           else fichier_out<<"\n          MSDinter    \n-------------------\n";
           FisParPop(identitybool,iLoc,fichier_out,MSgTotHez, NLocHez,NLoc, NtotLoc,SSiTotLoc,SSgTotLoc); //version tabF entier dans F_est: iLoc devient F_est_locIt
           delete_tabM_tabCode();
       	} else {
           if  (fichier_genepop->loci[iLoc]->getNumber()>0) {
               crunchLocTable(2,iLoc,fichier_out,&tabF); //version tabF entier dans CT_tests
               tabFtotabM(&tabF); //contient des new : delete ci dessous
    // ici on doit etre prete pour l'equi de lecture_Paires->lecture->remplit sfreqs etc.
               if (identitybool) fichier_out<<"\n          1-Qintra   1-Qinter      Fis\n-------------------------------------------\n";
               else fichier_out<<"\n           MSDintra   MSDinter    Rho(is)\n-------------------------------------------\n";
               FisParPop(identitybool,iLoc,fichier_out,MSgTotHez, NLocHez,NLoc, NtotLoc,SSiTotLoc,SSgTotLoc); // dans F_est: iLoc devient F_est_locIt
               delete_tabM_tabCode();
           } else fichier_out<<"No complete diploid genotypes\n";
        }
    } /*de l'iteration sur les locus*/
    fichier_out<<"\n\n============================================\n";
    fichier_out<<"Statistics per sample over all loci with at least two individuals typed:\n";
    if (estimDiploidBool) fichier_out<<"(DIPLOID loci only)\n"; else fichier_out<<"(HAPLOID loci only)\n";
    if (identitybool) fichier_out<<"Sample    1-Qintra   1-Qinter      Fis\n--------------------------------------------";
    else fichier_out<<"Sample     MSDintra   MSDinter    Rho(is)\n--------------------------------------------";
    for (size_t popit=0;popit<nb_sam;popit++) {
//1-Qintra
          if (NtotLoc[popit]>0 && estimDiploidBool)
              fichier_out<<endl<<left<<setw(11)<< fichier_genepop->pops[popit]->popName().substr(0,10)<<setw(11)<<SSgTotLoc[popit]/NtotLoc[popit];
          else
              fichier_out<<endl<<left<<setw(11)<< fichier_genepop->pops[popit]->popName().substr(0,10)<<"     -     ";
//1-Qinter
          if (NtotLoc[popit]-NLoc[popit]>0)
             if (estimDiploidBool) fichier_out<<setw(11)<<(SSgTotLoc[popit]+SSiTotLoc[popit])/(2.*NtotLoc[popit]);
             else fichier_out<<setw(11)<<SSiTotLoc[popit]/NtotLoc[popit];
          else
              fichier_out<<setw(11)<< "     -     ";
//Fis
          if (SSgTotLoc[popit]+SSiTotLoc[popit]>0 && estimDiploidBool) {
             ios_base::fmtflags old=fichier_out.setf(ios_base::fixed,ios_base::floatfield);
             fichier_out.precision(4); // appar il *faut* ca pour activer la bonne combin... le precision() prec ne joue plus ?
             fichier_out<<std::internal<<setw(7)<<(SSiTotLoc[popit] - SSgTotLoc[popit])/(SSiTotLoc[popit] + SSgTotLoc[popit]);
             fichier_out.setf(old,ios_base::floatfield); // necess pour SSD
             if (identitybool) fichier_out.precision(4); else fichier_out.precision(5);
          } else
             fichier_out<< "     -     ";
    }
    fichier_out<<"\n============================================\n\n\n======================";
    fichier_out<<"\n1-Qintra per sample over all loci with at least one diploid individual typed:\n";
    if (identitybool) fichier_out<<"Sample    1-Qintra       (unweighted average)\n----------------------";
    else fichier_out<<"Sample     MSDintra       (unweighted average)\n----------------------";
    for (size_t popit=0;popit<nb_sam;popit++) {
          if (NLocHez[popit]>0)
              fichier_out<<endl<<left<<setw(11)<< fichier_genepop->pops[popit]->popName().substr(0,10)<<setw(11)<<MSgTotHez[popit]/NLoc[popit];
          else
              fichier_out<<endl<<left<<setw(11)<< fichier_genepop->pops[popit]->popName().substr(0,10)<<"     -     ";
    }
    fichier_out<<"\n======================\n";
    fichier_out<<"\nNormal ending."<<endl;
	fichier_out.close();
    noR_cout<<"Normal ending."<<endl;
    noR_cout<<"Edit the file "<<nom_fic<<" for results"<<endl;
    if (warnOnExit) {
      if (!perf) ZeGenepopSound();
      if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
    }
return 0;
}

//---------------- option 6.1 -- 6.4 -------------------
//------------------------------------------------------------------------------
int Fstat(bool identitybool,int Indic) {
// Indic= 2 or 0 (pairwise or not)
using namespace NS_GP;
using namespace NS_F_est;
using namespace NS_FFF_slmt;
	nb_sam=fichier_genepop->pops.size();
	nb_locus=fichier_genepop->loci.size();
	double deno,Fis=0,Fst;
	vector<vector<unsigned long int> >tabF(nb_sam+2);  //va contenir les effectifs genos et les types alleliques (r??crite ? chaque locus)
	vector<vector<double> >FFF;  //table des F par locus, ou des num et denum des F par paires, selon !pairwise ou si
    SSiTot=0;SSgTot=0;MSi2P=0;MSg2P=0;MSp2P=0;MSi2Pw=0;MSg2Pw=0;
    ofstream fichier_out,fichier_mig;
    string nom_fic,nom_mig;
    nom_fic=gp_file;
    if (Indic==0) FFF.resize(nb_locus); else {
       FFF.resize(nb_sam);
       for (vector<vector<double> >::iterator ii=FFF.begin();ii<FFF.end(); ii++) {
           ii->resize(nb_sam);
           fill(ii->begin(),ii->end(),0);
       }
    }
    if (Indic==0) {if (identitybool) nom_fic+=".FST"; else nom_fic+=".RHO";}
    else {nom_fic+=".ST2";}
    fichier_out.open(nom_fic.c_str());
    if(!fichier_out.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    fichier_out.setf(ios::left,ios::adjustfield);
    if (Indic==2) {
        nom_mig=gp_file+char_mig;
        fichier_mig.open(nom_mig.c_str());
        if(!fichier_mig.is_open()){
            noR_cout<<"Error while opening file "<<nom_mig<<endl;
           genepop_exit(-1, "Error while opening file ");
        }
        fichier_mig.setf(ios::left,ios::adjustfield);
        fichier_mig.precision(4);
    }
    if (identitybool) {
        if (Indic==2) fichier_out<<"Genepop "<<getSetting("version")<<"\nPairwise Fst's\n\n";
        else fichier_out<<"Genepop "<<getSetting("version")<<"\nAllele frequency-based correlation (Fis, Fst, Fit)\n\n";
        fichier_out<<"One locus estimates following standard ANOVA as in Weir and Cockerham (1984)\n\n";
    } else {
        if (Indic==2) fichier_out<<"Genepop "<<getSetting("version")<<"\nPairwise Rhost's\n\n";
        else fichier_out<<"Genepop "<<getSetting("version")<<"\nAllele size based covariance (rhoIS, rhoST, rhoIT)\n\n";
        fichier_out<<"One locus estimates following standard ANOVA analogous to Weir and Cockerham (1984)\n";
        fichier_out<<"See e.g. Michalakis and Excoffier (1996) for formulas\n";
    }
    if (Indic==2) {
       fichier_mig<<"From file: "<<fichier_genepop->fileName<<endl<<nb_sam<<" populations\nGenetic statistic ";
       if (identitybool) fichier_mig<<"(Fst):\n"; else fichier_mig<<"(Rhost):\n";
    }
    fichier_out<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    fichier_out<<endl;
    fichier_out<<"Number of populations detected : "<<nb_sam<<endl;
    fichier_out<<"Number of loci detected        : "<<nb_locus<<endl;
    fichier_out<<endl;
    fichier_out.setf(ios_base::fixed,ios_base::floatfield);
    fichier_out.precision(6);

    if (Indic==2) {
       fichier_out<<"Indices for populations:\n----     -------------\n";
       {int bidon=1;
     	 for (vector<CPopulation *>::iterator p = fichier_genepop->pops.begin();p<fichier_genepop->pops.end();p++) {
    		fichier_out <<setw(9)<<bidon<< (*p)->popName().substr(0,8)<<endl;
    		bidon++;
         }
        } //bidon local
        fichier_out<<"----------------------\n\nEstimates for each locus:\n========================";
    }

// BOUCLE PRINCIPALE
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
        fichier_out<<endl<<"  Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
        fichier_out<<"---------------------------------\n";
        if (!identitybool) {
           if (taille.size()<=iLoc || taille[iLoc].size()==0)
              fichier_out<<" Allele sizes = allele names\n";
           else {
              fichier_out<<" Allele names: ";
              for (map<int,int>::iterator ptr=taille[iLoc].begin();ptr!=taille[iLoc].end();ptr++)
                  fichier_out<<setw(4)<<ptr->first;
              fichier_out<<"\n Allele sizes: ";
              for (map<int,int>::iterator ptr=taille[iLoc].begin();ptr!=taille[iLoc].end();ptr++)
                  fichier_out<<setw(4)<<ptr->second;
              fichier_out<<"\n(for unspecified alleles, size = name)\n";
           }
           fichier_out<<"---------------------------------\n";
        }
        if (Indic==2) {
           fichier_out<<"pop      ";
           for (size_t ii=1;ii<nb_sam;ii++) fichier_out<<left<<setw(8)<<ii;
           fichier_out<<endl;
        }
        if (fichier_genepop->coding[iLoc]<4) {// haploid
           if (Indic==0) crunchLocTable(1,iLoc,fichier_out,&tabF); //version tabF entier dans CT_tests
           else { //si pairwise, n'?crit pas la matrice des effectifs dans fichier_out
              ostringstream poubelle; //va ?ponger toutes les ?critures de crunchLocTable
              crunchLocTable(1,iLoc,poubelle,&tabF); //version tabF entier
           }
       	} else { //diploid
           if (Indic==0) crunchLocTable(2,iLoc,fichier_out,&tabF); //version tabF entier dans CT_tests
           else { //si pairwise, n'?crit pas la matrice des effectifs  dans fichier_out
              ostringstream poubelle; //va ?ponger toutes les ?critures de crunchLocTable
              crunchLocTable(2,iLoc,poubelle,&tabF); //version tabF entier
           }
        }
        tabFtotabM(&tabF); //contient des new : delete ci dessous
// ici on doit etre prete pour l'equi de lecture_Paires->lecture->remplit sfreqs etc.
        hierFstat(identitybool,Indic,iLoc,fichier_out,&FFF); // dans F_est: iLoc devient F_est_locIt
        delete_tabM_tabCode();
    } /*de l'iteration sur les locus*/
    fichier_out.precision(4);
    if (Indic==0) { ////////// NOT PAIRWISE ////////////////////
       if (estimDiploidBool) { //diploid multilocus
           fichier_out<<"Multilocus estimates for diploid data\n";
           if (identitybool) fichier_out<<"Locus           Fwc(is)     Fwc(st)     Fwc(it)\n";
           else fichier_out<<"Locus           Rho(is)     Rho(st)     Rho(it)\n";
           fichier_out<<"------------    -------     -------     -------\n";
           for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
               if (fichier_genepop->coding[iLoc]>3) {
                   fichier_out<<left<<setw(12)<<fichier_genepop->loci[iLoc]->locusName.substr(0,11)<<"    "<<std::internal;
                   fichier_out<<setw(7);
                   if (FFF[iLoc][0]>-665) fichier_out<<FFF[iLoc][0]; else fichier_out<<"   -   ";
                   fichier_out<<"     "<<setw(7);
                   if (FFF[iLoc][1]>-665) fichier_out<<FFF[iLoc][1]; else fichier_out<<"   -   ";
                   fichier_out<<"     "<<setw(7);
                   if (FFF[iLoc][2]>-665) fichier_out<<FFF[iLoc][2]; else fichier_out<<"   -   ";
                   fichier_out<<endl;
               }
           }
           fichier_out<<setw(12)<<"           All: "<<std::internal;
           fichier_out<<setw(7);
           if (SSiTot+SSgTot>0) {
              Fis=(SSiTot-SSgTot)/(SSiTot+SSgTot);
              fichier_out<<Fis;
           } else fichier_out<<"   -   ";
           fichier_out<<"     "<<setw(7);
           deno=MSp2P+MSi2Pw+MSg2Pw;
           if (deno>0) {
              Fst=(MSp2P-MSi2P)/deno;
              fichier_out<<Fst;
              fichier_out<<"     "<<setw(7)<<1-(1-Fis)*(1-Fst);
           } else fichier_out<<"   -           -";
           fichier_out<<"\n-----------------------------------------------";
       } else { // haploid multilocus
           fichier_out<<"Multilocus estimate for haploid data\n";
           if (identitybool) fichier_out<<"Locus           Fwc(st)\n";
           else fichier_out<<"Locus           Rho(st)\n";
           fichier_out<<"------------    -------\n";
           for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
               if (fichier_genepop->coding[iLoc]<4) {
                   fichier_out<<left<<setw(12)<<fichier_genepop->loci[iLoc]->locusName.substr(0,11)<<"    "<<std::internal;
                   fichier_out<<setw(7);
                   if (FFF[iLoc][1]>-665) fichier_out<<FFF[iLoc][1]; else fichier_out<<"   -   ";
                   fichier_out<<endl;
               }
           }
           fichier_out<<setw(12)<<"           All: "<<std::internal;
           fichier_out<<setw(7);
           deno=MSp2P+MSi2Pw;
           if (deno>0) {
              Fst=(MSp2P-MSi2P)/deno;
              fichier_out<<Fst;
           } else fichier_out<<"   -           -";
           fichier_out<<"\n-----------------------------------------------";
       }
    } else { ///////// PAIRWISE ////////////////////
       if (nb_sam<2) {
          fichier_out<<"\nOnly "<< nb_sam<<" population. No pairwise estimation.\n";
          goto lafin;
       }
       //ELSE
       if (estimDiploidBool) //diploid multilocus
          fichier_out<<"\nEstimates for all loci (diploid):\n=========================\n";
       else
          fichier_out<<"\nEstimates for all loci (haploid):\n=========================\n";
       fichier_out<<"pop      ";
       for (size_t ii=1;ii<nb_sam;ii++) fichier_out<<left<<setw(8)<<ii;
       fichier_out<<endl;
       for (global_pop_it=1;global_pop_it<nb_sam;global_pop_it++) {
            fichier_out<<setw(6)<<left<<global_pop_it+1<<std::internal;
        	for (global_pop2_it=0;global_pop2_it<global_pop_it;global_pop2_it++) {
                if (FFF[global_pop2_it][global_pop_it]>0) {
                   fichier_out<<setw(7)<<FFF[global_pop_it][global_pop2_it]/FFF[global_pop2_it][global_pop_it]<<" ";
                   fichier_mig<<setw(7)<<FFF[global_pop_it][global_pop2_it]/FFF[global_pop2_it][global_pop_it]<<" ";
                } else {
                  fichier_out<< "   -    ";
                  fichier_mig<< "   -    ";
                }
            }
            fichier_out<<endl;
            fichier_mig<<endl;
       }
       fichier_out<<"\nThe file "<<nom_mig<<" contains the matrix ready for further analysis";
    } ////////////////////////////////////////////////////////////
lafin:
    fichier_out<<"\nNormal ending."<<endl;
	fichier_out.close();
	if (Indic==2) {
       fichier_mig<<"Geographic distances:\n";
       fichier_mig.close();
    }
    noR_cout<<"Normal ending."<<endl;
    noR_cout<<"Edit the file "<<nom_fic<<" for results"<<endl;
    if (!perf) ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}



//---------------- option 6.5--6.8, isoldeFile, multiMigFile -------------------
//------------------------------------------------------------------------------
/*ca calcule le diverses statistiques sur les moments (Fst, a_r, e_r...)
genotip2() est juste un remplissage de tableau de donn?es
pairwMS() est un calcul de MS ? la Cockerham
Les diff?rents cas jouent dans
create_matrices()
o? joue ausi la pond?ration pour le bootstrap
*/

namespace varForBootstrapGenepop {
    vector<size_t> idxPloid;
   string nom_fich_mig;
   string nom_fich_tmp;
}

vector<double> creatMat_isolde(std::vector<double> inputwei) {
  /** inputwei as computed by bootstrapOverLoci should be a nbPloidBool-sized vector
  this is mapped to a nb_locus sized vector patched with 'zero' weights.
  Then create_matrices and isoldeproc are called.   **/
  using namespace varForBootstrapGenepop;
  using namespace NS_F_est;
  string nom_fich;
  bool was_first=_first_of_repl; // because the orginal value is needed after it is changed by isoldeproc
  for(size_t loc=0;loc<ABCweight.size();loc++) {ABCweight[loc]=0.;}
  for(size_t loc=0;loc<idxPloid.size();loc++) {ABCweight[idxPloid[loc]]=inputwei[loc];}
  vector<double> tt(3);
  if( _first_of_repl ) {
    nom_fich= nom_fich_mig;
  } else {
    nom_fich= nom_fich_tmp;
  }
  int indicnodist=create_matrices(nom_fich.c_str()); // this writes a .mig file in all cases
  if (indicnodist==-1) {
    delete_ptrs(); //deletes houla[][][] ETC;
    noR_cout<<"\nNo coordinates or equal coordinates for all samples;\n";
    noR_cout<<"No further analysis of isolation by distance.\n";
    noR_cout<<"(Return) to continue";
    if (pauseGP) cin.get();
    tt[0]=numeric_limits<double>::quiet_NaN();
    tt[1]=numeric_limits<double>::quiet_NaN();
    tt[2]=numeric_limits<double>::quiet_NaN();
    return(tt);
  } else {
    tt=isoldeproc(nom_fich.c_str()); //this reads a .mig file in all cases
    if( ! was_first ) {
      remove(nom_fich.c_str());
    }
  }
  return(tt);
}

//FR->FR a single function would be more efficient, but would require redefining the interface of bootstrapOverLoci
double slope_from_creatMat_isolde(std::vector<double> inputwei) {
    vector<double> tt=creatMat_isolde(inputwei);
    return(tt[1]);
}

double intercept_from_creatMat_isolde(std::vector<double> inputwei) {
    vector<double> tt=creatMat_isolde(inputwei);
    return(tt[0]);
}

double mean_from_creatMat_isolde(std::vector<double> inputwei) {
  vector<double> tt=creatMat_isolde(inputwei);
  return(tt[2]);
}

int isolde_etc(bool indiv) {
  using namespace varForBootstrapGenepop;
  using namespace NS_GP;
  using namespace NS_F_est;
  string nom_fich_CI;
  char ch='\0';
  vector<bool>ploidBool;
  set_phylipmatrix(phylipBool);
  // e_stat must be set before set_first_of_repl_ptrs() (case indiv==true below)
  if ( ! indiv) {
    _a_stat=false;_e_stat=false;
    // here singleGeneDiv could be true or false
  } else /** indiv**/ if ( ! estimDiploidBool) { /** haploid **/
  singleGeneDiv=true; // mandatory on haploid individual data
    if (_e_stat) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::R
#else
      cerr<<"\nRequesting '^e' estimator with haploid data is not meaningful. Check settings\n";
      if (pauseGP) getchar();
#endif

      genepop_exit(-1, "Requesting '^e' estimator with haploid data is not meaningful.");
    }
    // but performance=a... should work
    if (!identitySettingsBool) {
      effacer_ecran();
      noR_cout<<"Allele size-based analysis on individuals is NOT advised.\n";
      noR_cout<<"This combination of options has not been allowed.\nIdentity-based analysis will be performed.\n";
      identitySettingsBool=true;
      if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
    }
  }
  if (first_repl) {
    set_options(perf,indiv,identitySettingsBool); // indiv est argument de isolde_etc;
  }
  if (multiMigFileBool) { /** ad hoc multiple (per locus) matrices file; on sort de la fonction ? la fin de ce bloc **/
    effacer_ecran();
    _gotoxy(0,0);
    noR_cout<<"\nReading 'MultiMigFile'...";
    readMultiMigFile(isolde_file.c_str());
    noR_cout<<"\n...done.\n";
    if (datamatrix::nb_sam_migf<3) {
      noR_cout<<"\nOnly "<< datamatrix::nb_sam_migf<<" populations. No further analysis.\n";
      if (pauseGP) getchar();
      genepop_exit(-1, "Only populations. No further analysis.");
    } /*else*/
    outname=isolde_file; // 'outname' used by calcwritecorw called within bootstrapOverLoci
    CABCbootstrap ABC(multimig::nb_loc_migf);
    _first_of_repl=true;  // set to false in the first call of calcwritecorw()
    if (meanDiffBool) {
      ABC.bootstrapOverLoci(&mean," for INTERCEPT",isolde_file,false);
    } else {
      ABC.bootstrapOverLoci(&slope," for SLOPE",isolde_file);
      /// _first_of_repl can be set to true only once bc when rewrites the output file. But may not be a problem for multiMigFileBool case
      _first_of_repl=true;  // set to false in the first call of calcwritecorw()
      readMultiMigFile(isolde_file.c_str()); // it's better to reset the geo distances...
      ABC.bootstrapOverLoci(&intercept," for INTERCEPT",isolde_file,false);
    }
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    printf("\n results are stored in file %s\n",isolde_file.c_str());
#endif

    if (pauseGP) {
      noR_cout<<"\n(Return) to continue"<<endl; getchar();
    }
    return 0;
  } // on est sorti de isolde_etc
  /** ELSE **/
  if (isoldeFileBool) { // ad hoc single matrix file; on sort de la fonction ? la fin de ce bloc
    noR_cout<<"\nReading 'isolationFile'...\n";
    _first_of_repl=true; //FER 3/2012
    readGGFile(isolde_file.c_str()); // no conversion ! Note this is called again by isoldeproc below, whcih does not seem useful for an isoldeFile.
    _first_of_repl=false; //FER 3/2012
    noR_cout<<"\n...done.\n";
    if (datamatrix::nb_sam_migf<3) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::R
#else
      cerr<<"\nOnly "<< datamatrix::nb_sam_migf<<" populations. No further analysis.\n";
      if (pauseGP) getchar();
#endif

      genepop_exit(-1, "Only populations. No further analysis.");
    } /*else*/
  if (!phylipBool) {
    noR_cout<<"Analysis of Isolation by distance\n"; // for MultiMigFile
  }
  outname=isolde_file; // utilis? par calcwritecorw appel? par la ligne suivante
  // geo distances were read but not converted when the multimig file was read
  conversionFst(); // conversion since data have been read afresh
  if (phylipBool) {
    writepma(); // only uses genetic matrix
    return 0; //retour au menu
  }
  /** ELSE **/
  conversionGeo(); // note that this operates typeSelection on population (habitats) types
  writeGraOnly(isolde_file.c_str());
  isoldeproc(isolde_file.c_str()); //matrices -> regression estimates; calls conversion() for geo dist.
  if (cmp_nocase(typeSelection,"inter_all_types")!=0
      && cmp_nocase(typeSelection,"intra_all_types")!=0) mantelTest(true,mantelRankBool); // calcule Mantel aussi  //FR->FR small feature: in contrast to isolde_etc(), this is called even if there is no information.
  //delete_ptr_all_repl(); //deletes houla[][][] ETC;  // FR 030111: but set_first_of_repl_ptrs has not been called...
#ifdef COMPATIBILITYRCPP
  // Rcpp::R
#else
  printf("\n results are stored in file %s\n",isolde_file.c_str());
#endif

  if (pauseGP) {
    noR_cout<<"\n(Return) to continue"<<endl; getchar();
  }
  return 0;
  } // on est sorti de isolde_etc
  /***** ELSE : Genepop input file ******/
  { /**********************************************  bloc for clarity *****************************************/
  vector<double>t0(2);
    nb_sam=fichier_genepop->pops.size();
    nb_pair_sam=(nb_sam)*(nb_sam-1)/2;  // non transparent code as this is used rather deep in the inner 'dolls'
    nb_locus=fichier_genepop->loci.size();
    ABCweight.resize(nb_locus);
    /** output file names **/
    nom_fich_mig=gp_file;
    nom_fich_tmp=gp_file+char_tmp;
    nom_fich_CI=nom_fich_mig+char_iso; // whether perf or not;
    //performance output goes is result.CI through boot_result.open("result.CI"); in performance_main();
    nom_fich_mig+=char_mig;
    outname=nom_fich_CI; // 'outname' is where calcwritecorw (called within bootstrapOverLoci) writes
    /** identify loci with correct ploidy **/
    ploidBool.resize(0);
    idxPloid.resize(0);
    {   
      bool boule;
      size_t idx=0;
      for (vector<char>::iterator ii=fichier_genepop->coding.begin();ii<fichier_genepop->coding.end();ii++,idx++) {
        boule=((estimDiploidBool && (*ii)>3) ||(!estimDiploidBool && (*ii)<4));
        ploidBool.push_back(boule);
        if (boule) idxPloid.push_back(idx); // the aim is to replace the table of booleans by a table of indices
      }
    }
    size_t nbPloidLoc=idxPloid.size();
    /** checking appropriateness of data **/
    if (nbPloidLoc==0) {
      if (pauseGP) {
        if (estimDiploidBool) noR_cout<<"\n\nNo diploid locus (use estimationPloidy=Haploid setting?).\n\n";
        else noR_cout<<"\nNo haploid locus.\n";
        noR_cout<<"(Return) to continue"<<endl; getchar();
      }
      ofstream bootOut(nom_fich_CI.c_str(),ios::app);
      if ( ! bootOut.is_open()) {
#ifdef COMPATIBILITYRCPP
        // Rcpp::R
#else
        cerr<<"(!) From isolde_etc(): error while opening file "<<nom_fich_CI<<endl;
        if (cinGetOnError) cin.get();
#endif

        genepop_exit(-1, "(!) From isolde_etc(): error while opening file ");
      }
      if (estimDiploidBool) bootOut<<"\nNo diploid locus.\n";
      else bootOut<<"\nNo haploid locus.\n";
      bootOut.close();
      return -1;
    }
    if (nb_sam<3) {
      noR_cout<<"\nOnly "<< nb_sam<<" populations. No further analysis.\n";
      if (pauseGP) getchar();
      return 0;
    } /*else*/
    if (indiv) {
      for (vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();
           ii<fichier_genepop->pops.end();ii++) {
        if ((*ii)->inds.size() > 1) {
          noR_cout<<"\nSome 'pop'ulation contains more than one individual.";
          noR_cout<<"\nIndividual-based analyses cannot be performed. I exit...";
          noR_cout<<"\nIf you really want to perform this analysis on these data,";
          noR_cout<<"\nPut each individual in its own 'pop' (see menu option 8.4).";
          if (pauseGP) getchar();
          genepop_exit(-1, "Some 'pop'ulation contains more than one individual.");
        }
      }
    }
    /** settings not previously set **/
    // e_stat must be set before set_first_of_repl_ptrs()
    // !perf doit etre redondant avec pauseGP ?
    // for haploid there's never any choice
    if ((pauseGP) && (!perf) && (!IsolBDstatInSettingsBool)) {
      if (estimDiploidBool) {
        if (indiv) {
#ifdef COMPATIBILITYRCPP
          // Rcpp::R
#else
          noR_cout<<"\a\n Fit to (a) '^a' or to (e) '^e' statistic?\n";
          noR_cout<<"\n Enter 'a' or  'e':\n";
          cin>>ch;
          cin.ignore();
#endif

          if (ch=='A' || ch=='a') {_a_stat=true; _e_stat=false;}
          else {_a_stat=false; _e_stat=true;}
        } else { // diploid on pops
#ifdef COMPATIBILITYRCPP
          // Rcpp::R
#else
          noR_cout<<"\a\n (F) Only use pairwise Fst/(1-Fst) or (C) use a common denominator for all pairs?\n";
          noR_cout<<"\n Enter 'F' or  'C' ('C' may be more advisable, see Documentation):\n";
          cin>>ch;
          cin.ignore();
#endif

          if (ch=='C' || ch=='c') {singleGeneDiv=true;}
          else {singleGeneDiv=false;}
        }
      } else { //haploid
        if (indiv) { // no choice here
        } else { //haploid on pops
#ifdef COMPATIBILITYRCPP
          // Rcpp::R
#else
          noR_cout<<"\a\n (F) Only use pairwise Fst/(1-Fst) or (C) use a common denominator for all pairs?\n";
          noR_cout<<"\n Enter 'F' or  'C' ('C' may be more advisable, see Documentation):\n";
          cin>>ch;
          cin.ignore();
#endif

          if (ch=='C' || ch=='c') {singleGeneDiv=true;}
          else {singleGeneDiv=false;}
        }
      }
      //cout<<_a_stat<<" "<<_e_stat;
    }
    if (estimDiploidBool) {
      if (_a_stat) statname="'^a' statistic"; else if (_e_stat) statname="'^e' statistic"; else {
        if (singleGeneDiv) statname="'F/(1-F)'-like with common denominator"; else statname="F/(1-F)";
      }
    } else {
      if (_a_stat) statname="'^a'-like statistic"; else if (_e_stat) statname="'^e'-like statistic"; else {
        if (singleGeneDiv) statname="'F/(1-F)'-like with common denominator"; else statname="F/(1-F)";
      }
    }
    /** allocating pointers, the old way... **/
    set_ptrs(); /******************  ! ************************/ // Allocates houla[][][]   etc.
    /** fills per locus tables**/
    genotip2();
    pairwMS(ploidBool); //calcul MS analyse de variance. la prochaine etape est create_matrice
    {  // remove LOCUS files created by genotip2 and read by pairwMS (2017/02/17)
      for(size_t jfi=0;jfi<nb_locus;jfi++){
        stringstream stst;
        stst<<"LOCUS"<<jfi+1;
        remove(stst.str().c_str());
      }
    }

    /** geo coordinates **/
    if ( ! geoDistFromGeoFile) fichier_genepop->extract_coord_pop();
    /** phylip output case**/
    if (phylipBool) {
      // giving equal weights to all loci:
      for(size_t loc=0;loc<nb_locus;loc++) {
        if (ploidBool[loc]) ABCweight[loc]=1.0/nbPloidLoc; else ABCweight[loc]=0;
        //cout<<ploidBool[loc]<<" "<<" "<<nbPloidLoc<<" "<<ABCweight[loc];getchar();
      }
      //create matrix appelle ecriture pop_tot qui ecrit le .mig  (nom_fich_mig cr?? pr?c.)
      create_matrices(nom_fich_mig.c_str()); //->genetic matrix; geog matrix
      writepma(); // only uses genetic matrix
      remove(nom_fich_mig.c_str());
      noR_cout<<"Normal ending."<<endl;
      noR_cout<<"Edit the file "<<gp_file.c_str()<<".PMA for results"<<endl;
      if (!perf) ZeGenepopSound();
      if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
      delete_ptrs(); //deletes houla[][][] ETC;
      return 0; //retour au menu
    }
    /** ELSE **/
    /** point estimate and bootstrap **/
    effacer_ecran();
    _gotoxy(0,0);
    noR_cout<<"Analysis of Isolation by distance\n"; // for general case
    if (perf) remove(nom_fich_mig.c_str());
    CABCbootstrap ABC(idxPloid.size());
    _first_of_repl=true;  // set to false in the first call of calcwritecorw()
    /// _first_of_repl can be set to true only once bc it rewrites the output file
    /// and _first_of_repl must be run independently for 'mean' and BD analyses => they are exclusive alternatives
    if (meanDiffBool) {
      noR_cout<<"Analysis of mean differentiation\n";
      ABC.bootstrapOverLoci(&mean_from_creatMat_isolde," for mean",nom_fich_CI,false);
    } else {
      noR_cout<<"Analysis of Isolation by distance\n"; // for general case
      ABC.bootstrapOverLoci(&slope_from_creatMat_isolde," for SLOPE",nom_fich_CI); // modifies datamatrix::data; using it as scratch!
      if (perf) {	 //ecriture des resultats
        _gotoxy(5,9);
        noR_cout<<"Slope= "<<ABC.t0<<" ["<<ABC.tinf<<" , "<<ABC.tsup<<"]       ";  //screen
        if(NS_GP_PERF::isample==NS_GP_PERF::JobMin) { // comments at beginning of multisample output file
          if (_a_stat) statname="'^a' slope   "; else if (_e_stat) statname="'^e' slope   "; else statname="F/(1-F) slope";
          boot_result<<"# "<<statname;
          boot_result<<"       Lower bound (slope)  Upper bound (slope)";
          if ( ! std::isnan(testPointslope)) boot_result<<" Pvalue for slope="<<testPointslope;
          boot_result<<endl;
          boot_result.precision(15);
        }
        boot_result<<ABC.t0<<"   "<<ABC.tinf<<"   "<<ABC.tsup; //file (order sup/inf inverted 01/2011)
        if ( ! std::isnan(testPointslope)) boot_result<<"    "<<ABC.testPointPvalue;
        boot_result<<endl;
      } // else non-perf file output made within bootstrapOverLoci
      /** Mantel uses namespace datamatrix hence must be called before bootstrap for intercept**/
      if ( ! std::isnan(ABC.t0) 
             && cmp_nocase(typeSelection,"inter_all_types")!=0 
             && cmp_nocase(typeSelection,"intra_all_types")!=0 ) mantelTest(false,mantelRankBool); // calcule Mantel aussi
      /** **/
      ABC.bootstrapOverLoci(&intercept_from_creatMat_isolde," for INTERCEPT",nom_fich_CI,false);
    }
    //_first_of_repl=true;  // set to false in the first call of calcwritecorw()
    delete_ptrs(); //deletes houla[][][] ETC;
    if (!perf) {
      noR_cout<<"\n**** normal ending****";
    }
    noR_cout<<"\nResults are stored in file "<<nom_fich_CI;
    if (pauseGP) { noR_cout<<"\n(Return) to continue"<<endl; getchar();}
    return(0);
  } /** END bloc GenepopFile for clarity**/
}


vector<double> bootstrapNullAllele(CGenotypes *popGenos,const size_t iLoc, const size_t iPop,
                                   map<ssize_t,double>& genocopy, // copy locale en <double> de la map de popGenos;
                                   // ssize_t to avoid conversions from getNext() that can be -1
                                     double& genocopySum, //differs from original sample for delta-perturbed samples !
                                     int& typeMax, // Max allelic type over pops at given locus
                                     string& failure,
                                     double& gauss_inf,
                                     double& gauss_sup,
                                     ofstream& fichier_out) {
// direct weighing of genotypic counts (not most straightforward but more efficient)
//    string nom_fic;
    bool ABCind_is_new_type;
    double *delta, *dt, *ddt;
    size_t nbInd=popGenos->getSum();
    ssize_t geno; // mismatch sur type map key est fatal
    size_t genocumul;
	delta=new double[nbInd];
	dt=new double[nbInd];
	ddt=new double[nbInd];
	double t0;
	vector<double>estimates(2);
	vector<double> tminput(nbInd);
	vector<double> tpinput(nbInd);
    double epsn,tm,tp,cq,bhat,curv,tinf,tsup,bidul25,bidul975,machin,z;//ABC bootstrap variables
	double sigmahat=0.0;
	double ahat=0.0;
    double epsn_value=0.001/nbInd;
    failure=""; //string from NS_Null storing cause of failure
//    _gotoxy(0,1);
//    cout<<"CI[ "<<gauss_inf<<"- "<<gauss_sup<<" ] will be computed";

    if (nbInd<2) {
       estimates.resize(0); // indicator failed computation
       failure="Fewer than two individuals     ";
	   return estimates;
    }

     //point estimate for original data (ABCweigth constant)
    genocopy.clear();
    genocopySum=0;
    popGenos->resetIterator();
    while ((geno=popGenos->getNext()) >=0) {
        genocopy[geno]=popGenos->getEffective(geno);
        genocopySum+=genocopy[geno];
        //        cout<<geno<<" "<<genocopy[geno]<<endl;
    }
    //ATTENTION si RHS a size()=1 LHS est reduit ? size=1
    estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out);
    t0=estimates[0];
//    _gotoxy(0,8);
//	cout<<" Computing confidence interval... beginning              ";
    typeMax=fichier_genepop->loci[iLoc]->alleleMax; //! max sur ttes les pops = NULL allele


	//pour borne inf
	epsn=-epsn_value;
	for(size_t ABCind=0;ABCind<nbInd;ABCind++) {
        genocopy.clear();
        popGenos->resetIterator();
        genocumul=0;
        ABCind_is_new_type=false; 
        for(size_t ind=0;ind<nbInd;ind++) { // ind does not refer to the lines in the GP file
          /*
           * Suppose we have 10 individual and genotypes in counts 4 3 2 1. We run over 'ind'
           * and thus 4 times over the first genotyep, 3 over the second... 
           * We meet a new (geno)type for ind = 0 4 7 9 as tested by (ind>=genocumul)
           * Then (for each new type: ind = 0 4 7 9) we change its count to getEffective(geno)*(1.-epsn)  
           * Further (ind==ABCind) is evaluated : 
           *  * conditionally on ind= 0 4 7 9 and thus ABCind_is_new_type is set to true only for ABCind= 0 4 7 9. 
           *     Since ABCind 0 1 2 3 have the same genotype by def, there is no need to recompute 'estimates' for _ABCind_ = 1 2 3
           *     Since ind runs on all values for each ABCind, all genocopy[geno] are reduced in this way.
           *  * unconditionnally, for all ind and then genocopy[geno] is increased by genocopy[geno]. But this affects estimates 
           *     only if (ABCind_is_new_type)...
           *     
           * Here the weigths are the genotype counts and sum to nbInd they are allready by a *factor* (1.-epsn) 
           * hence nbInd*epsn must be added on the genotype of the ABCind     
           */
            if (ind>=genocumul) { // new geno on indiv ind;  // could test == instead of >= ?
              genocumul+=popGenos->getEffective((geno=popGenos->getNext()));
              genocopy[geno]=popGenos->getEffective(geno)*(1.-epsn); // only once for each geno
              if (ind==ABCind) { 
                ABCind_is_new_type=true; 
                genocopy[geno]+=epsn*nbInd; // only once in the ind loop, after genocopy[<geno of ABCind>] has been initialized 
              } // could test ABCind=genocumul IF first line, before updating genocumul ?
            } 
        } // end loop on each ind
        if (ABCind_is_new_type) estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out); 
        // elee keeps estimate for previous ABCind
        tminput[ABCind]=estimates[0];
    }

	//pour borne sup
	epsn=epsn_value;
	for(size_t ABCind=0;ABCind<nbInd;ABCind++) {
        genocopy.clear();
        popGenos->resetIterator();
        genocumul=0;
        ABCind_is_new_type=false;
        for(size_t ind=0;ind<nbInd;ind++) {
            if (ind>=genocumul) { // individual -> genotype count
                genocumul+=popGenos->getEffective((geno=popGenos->getNext()));
                genocopy[geno]=popGenos->getEffective(geno)*(1.-epsn);
                if (ind==ABCind) {
                  ABCind_is_new_type=true;
                  genocopy[geno]+=epsn*nbInd;
                }
            }
        }
       if (ABCind_is_new_type) estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out);
       tpinput[ABCind]=estimates[0];
    }
	for(size_t ind=0;ind<nbInd;ind++) {
	  dt[ind]=(tpinput[ind]-tminput[ind])/(2*epsn);
		ddt[ind]=(tpinput[ind]+tminput[ind]-2*t0)/pow(epsn,2);
	}
	for(size_t ind=0;ind<nbInd;ind++) {
        sigmahat+=pow(dt[ind],2);
    }
	sigmahat=sqrt(sigmahat)/nbInd; // low if no info, or quite a lot such that dt for each indiv low ?
    if (sigmahat<((1e-7)/sqrt(nbInd))) { // likely cause being all dt =0 to precision error
       estimates.resize(0); // indicator failed computation
       failure="sigmahat suspiciously low                   ";
	   return estimates;
    }
	for(size_t ind=0;ind<nbInd;ind++) ahat+=pow(dt[ind],3);
	ahat=ahat/(6*pow(double(nbInd),3)*pow(sigmahat,3));
    if (std::isnan(ahat)) { // likely cause being all dt =0 should now not pass previous test
       estimates.resize(0); // indicator failed computation
       failure="ahat is not a number                   ";
	   return estimates;
    }
	for(size_t ind=0;ind<nbInd;ind++) {
        delta[ind]=dt[ind]/(pow(double(nbInd),2)*sigmahat);
	}
	epsn=-epsn_value;

    genocopy.clear();
    genocumul=0;
    genocopySum=0; // must be recomputed for delta-perturbed sample.
    popGenos->resetIterator();
    for(size_t ind=0;ind<nbInd;ind++)
        if (ind>=genocumul) { // individual -> genotype count
            genocumul+=popGenos->getEffective((geno=popGenos->getNext()));
            // w0 = 1/nbInd a 1/nBind + delta[ind]*epsn donc chaque indiv a poids relatif 1-epsn*delta*nbInd
            genocopy[geno]=popGenos->getEffective(geno)*(1.+delta[ind]*epsn*nbInd);
            genocopySum+=genocopy[geno];
        }
    estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out);
    tp=estimates[0];
//cout<<estimates[0];getchar();
    genocopy.clear();
    genocumul=0;
    genocopySum=0;
    popGenos->resetIterator();
    for(size_t ind=0;ind<nbInd;ind++)
        if (ind>=genocumul) { // individual -> genotype count
            genocumul+=popGenos->getEffective((geno=popGenos->getNext()));
            genocopy[geno]=popGenos->getEffective(geno)*(1.-delta[ind]*epsn*nbInd);
            genocopySum+=genocopy[geno];
        }
    estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out);
    tm=estimates[0];
//cout<<estimates[0];getchar();

//	printf("\r Computing confidence interval... about %.2g %% done",100*(2*(double)nb_locus+2)/(2*(double)nb_locus+5));
	cq=(tp+tm-2*t0)/(2*sigmahat*pow(epsn,2));
	bhat=0.0;
	for(size_t ind=0;ind<nbInd;ind++) bhat+=ddt[ind];
	bhat/=(2*pow(double(nbInd),2));
	curv=(bhat/sigmahat)-cq;
	if (fabs(curv)>37) {
	   	 estimates.resize(0); // indicator failed computation
         failure="curv is too large                  ";
	   	 return estimates; // not exit() !
	} else {
       machin=2*ndtr(ahat)*ndtr(-curv);
	   if((machin>=1)||(machin<=0.0)) {//z=INFINITY; a priori machin est une proba...
	   	 estimates.resize(0); // indicator failed computation
         failure="ndtri argument is not within ]0,1[";
	   	 return estimates; // not exit() !
       } else z=ndtri(machin);
    }
	bidul25=(z+ndtri(gauss_inf))/(pow(1-ahat*(z+ndtri(gauss_inf)),2));
	bidul975=(z+ndtri(gauss_sup))/(pow(1-ahat*(z+ndtri(gauss_sup)),2));

    genocopy.clear();
    genocumul=0;
    popGenos->resetIterator();
    genocopy.clear();
    genocumul=0;
    genocopySum=0;
    popGenos->resetIterator();
    for(size_t ind=0;ind<nbInd;ind++)
        if (ind>=genocumul) { // individual -> genotype count
            genocumul+=popGenos->getEffective((geno=popGenos->getNext()));
// this (and below) is where the w0 : nbInd stuff matters
            genocopy[geno]=popGenos->getEffective(geno)*(1.+delta[ind]*bidul25*nbInd);
            genocopySum+=genocopy[geno];
            if (genocopy[geno]<0) {
               estimates.resize(0); // indicator failed computation
               failure="genocopy<0 in tinf computation             ";
        	   return estimates;
            }
    }
    estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out);
    tinf=estimates[0];

//	printf("\r Computing confidence interval... about %.2g %% done",100*(2*(double)nbInd+4)/(2*(double)nbInd+5));

    genocopy.clear();
    genocumul=0;
    genocopySum=0;
    popGenos->resetIterator();
    for(size_t ind=0;ind<nbInd;ind++)
        if (ind>=genocumul) { // individual -> genotype count
            genocumul+=popGenos->getEffective((geno=popGenos->getNext()));
            genocopy[geno]=popGenos->getEffective(geno)*(1.+delta[ind]*bidul975*nbInd);
            genocopySum+=genocopy[geno];
            if (genocopy[geno]<0) {
               estimates.resize(0); // indicator failed computation
               failure="genocopy<0 in tsup computation";
        	   return estimates;
            }
        }
    estimates=estimNullLocPop(iLoc,iPop,false,genocopy,genocopySum,typeMax,fichier_out);
    tsup=estimates[0];
//    if (tsup<tinf) {tsup=tinf;tinf=estimates[0];} //!
    if (tinf<0 || tsup>1) {
       estimates.resize(0); // indicator failed computation
       failure="tinf<0 or tsup>0                      ";
	   return estimates;
    }
    delete[] delta;
    delete[] dt;
    delete[] ddt;
//pcq'il a pu etre reduit par les =estimNull...
    estimates.resize(2);
    estimates[0]=tinf;
    estimates[1]=tsup;
return estimates;
}
/////////////////////////////////////////



vector<double> estimNullLocPop(const size_t iLoc, const size_t iPop, const bool printBool,
                                 map<ssize_t,double>& genocopy, // copy locale en <double> de la map de popGenos;
                                 double& genocopySum, //differs from original sample for delta-perturbed samples !
                                 int& typeMax,
                                 ofstream& fichier_out) {
   vector<double>results(0);
   int homtype,heztype,nulltype;
   map<int,CAllele*> * allelesPtr;
   double Broo_A,Broo_B,nullhom,Ho,He,beta,betanext;
   double diff,numer,denom,fauxnull;
   map<int,double> estP,estPnext;
   //map<int,double>::iterator copyiter=genocopy.begin();
   allelesPtr=&(fichier_genepop->pops[iPop]->loci[iLoc]->alleles); //! alleles visibles dans la pop iPop
   double genosum=genocopySum; // on va en retirer les visible null hez.
   if (nullIgnoredBool) { // method Chakraborty et al 92 ... secret option!
   //observed freq of obvious hez
      Ho=0.;
      for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++) {
           if (fichier_genepop->coding[iLoc]==4) {
              homtype=101*jj->first;
           } else {
              homtype=1001*jj->first;
           }
           if (jj->first==typeMax) {//on retire les visible null homoz
              genosum-=genocopy[homtype];
           } else Ho+=genocopy[homtype];
      }
      //"expected" freq of obvious hez
      He=0;
      for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++)
           for (map<int,CAllele*>::iterator kk=(*allelesPtr).begin();kk!=jj;kk++) {
               if (jj->first!=typeMax) He+=(jj->second->getEffective())*(kk->second->getEffective());
           }
      if (genosum>0 && He+Ho>0) {
          He/=(2.*genosum*genosum); //2 pq but divided by nb genotypes
          Ho=1.-Ho/genosum;
          results.push_back((He-Ho)/(He+Ho));
      } else {
         results.push_back(-1);
      }
   } else if (Brookfield96Bool) { // method Brookfield 96
      //null
      if (fichier_genepop->coding[iLoc]==4) {
          homtype=101*typeMax;
      } else {
          homtype=1001*typeMax;
      }
      nullhom=genocopy[homtype]/genocopySum;
//if (iLoc==2 && iPop==3) {cout<<" "<<genocopy[101]<<" "<<genocopy[202]<<" "<<nullhom;getchar();}
      //"expected" freq of obvious hez
      He=0;
      for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++)
           for (map<int,CAllele*>::iterator kk=(*allelesPtr).begin();kk!=jj;kk++) {
               if (jj->first!=typeMax) He+=(jj->second->getEffective())*(kk->second->getEffective());
//cout<<He;getchar();
           }
      He/=(2.*genocopySum*genocopySum); //2 pq but divided by nb genotypes
      He/=((1.-nullhom)*(1.-nullhom)); // avant derni?re de l'article
      //observed freq of obvious hez
      Ho=0.;
      for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++) {
           if (fichier_genepop->coding[iLoc]==4) {
              homtype=101*jj->first;
           } else {
              homtype=1001*jj->first;
           }
           Ho+=genocopy[homtype];
//if (iLoc==2 && iPop==3) {cout<<Ho<<" "<<homtype<<" "<<genocopy[homtype]<<" "<<genocopySum;getchar();}
      }
      Ho=1.-Ho/genocopySum;
      Broo_A=He*(1.+nullhom)-Ho;
      Broo_B=4.*nullhom*(1.-He*He);
      estP[typeMax]=(Broo_A+sqrt(Broo_A*Broo_A+Broo_B))/(2.*(1.+He));
//if (iLoc==2 && iPop==3) {cout<<estP[typeMax]<<" "<<He<<" "<<Ho<<" "<<nullhom;getchar();}
      results.push_back(estP[typeMax]);
   } else { // full ML using EM algo
       if (fichier_genepop->pops[iPop]->loci[iLoc]->alleleExists(typeMax)) // null allele in allele map
          for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++)
              estP[jj->first]=1./(*allelesPtr).size();
       else { // cas ou null allele pas le max *de cette pop l?* (not in Genepop->3.4 included)
          for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++)
              estP[jj->first]=1./((*allelesPtr).size()+1);
          estP[typeMax]=1./((*allelesPtr).size()+1); //(otherwise no convergence for zero initial freq)
       }
       if (NonNullfailuresBool && (*allelesPtr).size()>2) beta=0.05; else beta=0;
       diff=1.;
       while (diff>1e-10) {
//cout<<diff;
//Expectation step: expected # of true hoz and hez among apparent hoz
           for(map<int,double>::iterator jj=estP.begin();jj!=estP.end();jj++) {
               if(jj->first!=typeMax) {
                   denom=jj->second*jj->second+2.*jj->second*estP[typeMax]; //p_i^2+2p_ip_NULL
// ordre vient de la def de fillGenoypes:... declareGenotype(...getMaxAllele(geno) * 100 + ...getMinAllele(geno));
                   if (fichier_genepop->coding[iLoc]==4) {
                      homtype=101*jj->first;
                      heztype=100*typeMax+jj->first;
                   } else {
                      homtype=1001*jj->first;
                      heztype=1000*typeMax+jj->first;
                   }
                   numer=genocopy[homtype]+genocopy[heztype]; //geno_ii+geno_iNULL (the data)
                   genocopy[homtype]=numer*jj->second*jj->second/denom; //numer*p_i^2/denom
                   genocopy[heztype]=numer*2.*jj->second*estP[typeMax]/denom; //numer*2p_ip_NULL/denom
               } //else null homoz, genocopy is not changed
           } //boucle sur alleles
//Maximisation step:= new allele frequencies by simple gene counting since this is ML for complete data
//note that when apparent null homoz are discarded as unreliable, ML differs from gene counting
           for(map<int,double>::iterator jj=estP.begin();jj!=estP.end();jj++) estPnext[jj->first]=0;
           for(map<ssize_t,double>::iterator jj=genocopy.begin();jj!=genocopy.end();jj++) { //allelic from genotypic
               estPnext[int(minAllele(jj->first,fichier_genepop->coding[iLoc]))]+=jj->second;
               estPnext[int(maxAllele(jj->first,fichier_genepop->coding[iLoc]))]+=jj->second;
           }
           if (fichier_genepop->coding[iLoc]==4) {
              nulltype=101*typeMax;
           } else {
              nulltype=1001*typeMax;
           }
           denom=beta+(1-beta)*estP[typeMax]*estP[typeMax];
           if (denom>0) {
               for (map<int,CAllele*>::iterator jj=(*allelesPtr).begin();jj!=(*allelesPtr).end();jj++) {
                   if (jj->first!=typeMax) {
                      fauxnull=2*genocopy[nulltype]*estP[jj->first]*beta/denom;
                      estPnext[jj->first]+=fauxnull;
                      estPnext[typeMax]-=fauxnull;
                   }
               }
               betanext=beta/denom*genocopy[nulltype]/genocopySum;
           } else {betanext=0;} //no apparent null hom
           for(map<int,double>::iterator jj=estPnext.begin();jj!=estPnext.end();jj++) jj->second/=(2.*genocopySum);
//Convergence criterion
           diff=0;
           for (map<int,double>::iterator jj=estP.begin();jj!=estP.end();jj++) {
               diff+=fabs(estPnext[jj->first]-jj->second);
               jj->second=estPnext[jj->first];
           }
           diff+=fabs(betanext-beta);
           beta=betanext;
       } // fin while
       results.push_back(estP[typeMax]);
       if ((*allelesPtr).size()>2) results.push_back(beta); else results.push_back(-1);
       if (printBool) {
           fichier_out<<"Population: "<<fichier_genepop->pops[iPop]->popName()<<"\n-----------------------------------------\n";
           fichier_out<<"Allele   Freq.     Homoz.     Null Heter.\n";
           for (map<int,double>::iterator jj=estP.begin();jj!=estP.end();jj++) {
               if (fichier_genepop->coding[iLoc]==4) {
                  homtype=101*jj->first;
                  heztype=100*typeMax+jj->first;
               } else {
                  homtype=1001*jj->first;
                  heztype=1000*typeMax+jj->first;
               }
               if(jj->first!=typeMax) {
                  fichier_out<<"  "<<setw(7)<<jj->first;
                  fichier_out<<setw(7)<<jj->second<<"   "<<setw(7)<<genocopy[homtype];
                  fichier_out<<"    "<<setw(7)<<genocopy[heztype]<<endl;
               } else fichier_out<<" Null    "<<setw(7)<<jj->second<<endl;
           }
           if (NonNullfailuresBool) fichier_out<<"Genotyping failure rate= "<<beta<<endl;
       } //printBool
   }
//cout<<results[0];getchar();
return results;
}

//---------------- option 8.1 NEW-------------------
int estimNull() {
  bool locdebug=false;
using namespace NS_GP;
 vector<vector<double> >tabFnull; // tableau frequence null allele...
 vector<vector<double> >tabBeta; //
 char coding;
 vector<vector<vector<double> > >bounds; // tableau frequence null allele...
 ofstream fichier_out;
 map<ssize_t,double> genocopy; // copy locale en <double> de la map de popGenos;
 double genocopySum; //differs from original sample for delta-perturbed samples !
 int typeMax; // Max allelic type over pops at given locus
 string failure;
 double gauss_inf;
 double gauss_sup;

    gauss_inf=(1.0-widthCI)/2.0;
    gauss_sup=1.0-gauss_inf;
	size_t nb_sam=fichier_genepop->pops.size();
	size_t nb_locus=fichier_genepop->loci.size();
	ssize_t genotype;
    const vector<double>bidon(0);
	vector<double>estimates;
    CGenotypes popGenotypes;
    string nom_fic;
    nom_fic=gp_file+".NUL";
    fichier_out.open(nom_fic.c_str());
    if(!fichier_out.is_open()){
        noR_cout<<"Error while opening file "<<nom_fic<<endl;
       genepop_exit(-1, "Error while opening file ");
    }
    fichier_out.setf(ios::left,ios::adjustfield);
    if (nullIgnoredBool) //secret option
       fichier_out<<"Genepop "<<getSetting("version")<<"\nEstimation of null allele frequency by Chakraborty et al's 1992 method.\n\n";
    else if (Brookfield96Bool)
       fichier_out<<"Genepop "<<getSetting("version")<<"\nEstimation of null allele frequency by Brookfield's 1996 method.\n\n";
    else {
         fichier_out<<"Genepop "<<getSetting("version")<<"\nMaximum likelihood estimation of null allele frequency\n\n";
         fichier_out<<"EM algorithm (Dempster, Laird and Rubin, 1977)\n\n";
    }
    fichier_out<<"File: "<<fichier_genepop->fileName<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    fichier_out<<endl;
    fichier_out<<"Number of populations detected : "<<nb_sam<<endl;
    fichier_out<<"Number of loci detected        : "<<nb_locus<<endl;
    fichier_out<<endl;
    fichier_out.setf(ios_base::fixed,ios_base::floatfield);
    fichier_out<<setprecision(4);



    _gotoxy(0,11);
    noR_cout<<"Computing confidence intervals... Locus        Pop       \n";
    tabFnull.resize(nb_locus);
    tabBeta.resize(nb_locus);
    bounds.resize(nb_locus);
// BOUCLE PRINCIPALE
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
      if (locdebug) RnoR_cerr<<endl<<iLoc<<" ";
        tabFnull[iLoc].resize(0);
        tabBeta[iLoc].resize(0);
        bounds[iLoc].resize(0);
        _gotoxy(40,11);
        noR_cout<<iLoc+1;
        if (!nullIgnoredBool && !Brookfield96Bool) {
           fichier_out<<endl<<"Locus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
           fichier_out<<"=================================\n";
        }
        if ((coding=fichier_genepop->coding[iLoc])<4) {
            if (!nullIgnoredBool && !Brookfield96Bool)
                 fichier_out<<"Computation not applicable on haploid data"<<endl<<endl;
       	} else {
           typeMax=fichier_genepop->loci[iLoc]->alleleMax; //! max sur ttes les pops = NULL allele
       	  if (locdebug) RnoR_cerr<<typeMax<<" ";
       	  for (size_t iPop=0;iPop<fichier_genepop->pops.size();iPop++) {
               _gotoxy(51,11);
               noR_cout<<iPop+1<<"     \n"<<flush;
               popGenotypes.clear();
               popGenotypes.fillGenotypes(iLoc, fichier_genepop->pops[iPop],coding);
               if (popGenotypes.getSum()==0)  {
                  fichier_out<<"Population: "<<fichier_genepop->pops[iPop]->popName();
                  fichier_out<<"\n-----------------------------------------\n";
                  fichier_out<<"No information for this locus in this population\n";
                  fichier_out<<"(No full genotype)\n";
                  tabFnull[iLoc].push_back(-1);
                  tabBeta[iLoc].push_back(-1);
                  bounds[iLoc].push_back(bidon);
               } else if (NonNullfailuresBool && fichier_genepop->pops[iPop]->loci[iLoc]->getNumber()<3) {
// not enough alleles for joint estimation of null allele freq and beta
                  fichier_out<<"Population: "<<fichier_genepop->pops[iPop]->popName();
                  fichier_out<<"\n-----------------------------------------\n";
                  fichier_out<<"No information for this locus in this population\n";
                  fichier_out<<"(Less than 3 alleles)\n";
                  tabFnull[iLoc].push_back(-3);
                  tabBeta[iLoc].push_back(-1);
                  bounds[iLoc].push_back(bidon);
               } else if (fichier_genepop->pops[iPop]->loci[iLoc]->getNumber()<2) {
// not enough alleles for null allele freq estimation
                  fichier_out<<"Population: "<<fichier_genepop->pops[iPop]->popName();
                  fichier_out<<"\n-----------------------------------------\n";
                  fichier_out<<"No information for this locus in this population\n";
                  fichier_out<<"(Less than 2 alleles)\n";
                  tabFnull[iLoc].push_back(-2);
                  tabBeta[iLoc].push_back(-1);
                  bounds[iLoc].push_back(bidon);
               } else {
                 popGenotypes.resetIterator();
                   genocopy.clear(); // clears the map
                   genocopySum=popGenotypes.getSum();
                   while ((genotype = popGenotypes.getNext()) >= 0) //sur les genos complets seulement
                         genocopy[genotype]=popGenotypes.getEffective(genotype); //copy en <double>
                   estimates.resize(2);
                   estimates=estimNullLocPop(iLoc,iPop,true,genocopy,genocopySum,typeMax,fichier_out);
                   tabFnull[iLoc].push_back(estimates[0]);
                   if (NonNullfailuresBool) tabBeta[iLoc].push_back(estimates[1]);
                   estimates=bootstrapNullAllele(&popGenotypes,iLoc,iPop,
                                                   genocopy, // copy locale en <double> de la map de popGenos;
                                                   genocopySum, //differs from original sample for delta-perturbed samples !
                                                   typeMax, // Max allelic type over pops at given locus
                                                   failure,
                                                   gauss_inf,
                                                   gauss_sup,fichier_out);
                   bounds[iLoc].push_back(estimates); // estimates a longueur nulle si failed
               }
               if (!nullIgnoredBool && !Brookfield96Bool) fichier_out<<endl;
//if (failure.size()>0) {cout<<failure;getchar();}
           } //iter su pops
        } // test si locus diploide
    } /*de l'iteration sur les locus*/
    fichier_out<<"\n\n(Locus by population) table of estimated null allele frequencies\n";
    fichier_out<<"================================================================\n";
    fichier_out<<"Locus:     Populations (! names truncated to 6 characters):\n           ";
    for (vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++)
        fichier_out<<setw(7)<<(*ii)->popName().substr(0,6);
    fichier_out<<"\n           -----------------------------------------------------\n";
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
        fichier_out<<setw(11)<<fichier_genepop->loci[iLoc]->locusName.substr(0,10);
        if (tabFnull[iLoc].size()==0) fichier_out<<"No information for this locus\n";
        else for(vector<double>::iterator jj=tabFnull[iLoc].begin();jj<tabFnull[iLoc].end();jj++)
             if ((*jj)>-0.5) fichier_out<<setw(7)<<(*jj); else fichier_out<<"No inf ";
        fichier_out<<endl;
    }
    fichier_out<<"================================================================\n";
    if (NonNullfailuresBool) {
        fichier_out<<"\n\n(Locus by population) table of estimated genotyping failure rate (beta).\n";
        fichier_out<<"================================================================\n";
        fichier_out<<"Locus:     Populations (! names truncated to 6 characters):\n           ";
        for (vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++)
            fichier_out<<setw(7)<<(*ii)->popName().substr(0,6);
        fichier_out<<"\n           -----------------------------------------------------\n";
        for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
            fichier_out<<setw(11)<<fichier_genepop->loci[iLoc]->locusName.substr(0,10);
            if (tabBeta[iLoc].size()==0) fichier_out<<"No information for this locus\n";
            else for(vector<double>::iterator jj=tabBeta[iLoc].begin();jj<tabBeta[iLoc].end();jj++)
                 if ((*jj)>-0.5) fichier_out<<setw(7)<<(*jj); else fichier_out<<"No inf ";
            fichier_out<<endl;
        }
        fichier_out<<"================================================================\n";
    }
    fichier_out<<"\n\nConfidence intervals for null allele frequencies\n";
    fichier_out<<"=================================================\n";
    fichier_out<<"                       Frequency   "<<setw(7)<<gauss_inf<<"  "<<setw(7)<<gauss_sup<<"\n";
    fichier_out<<"Locus      Population   estimate   bound    bound";
    fichier_out<<"\n-------------------------------------------------\n";
    for (size_t iLoc=0;iLoc<nb_locus;iLoc++) {
        fichier_out<<setw(11)<<fichier_genepop->loci[iLoc]->locusName.substr(0,10);
        if (tabFnull[iLoc].size()==0) fichier_out<<"No information for this locus\n";
        else {
           vector<vector<double> >::iterator bb=bounds[iLoc].begin();
           vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();
           for(vector<double>::iterator jj=tabFnull[iLoc].begin();jj<tabFnull[iLoc].end();jj++,bb++,ii++) {
             if (jj!=tabFnull[iLoc].begin()) fichier_out<<setw(11)<<" ";  //additional lines for same locus
             fichier_out<<setw(12)<<(*ii)->popName().substr(0,11);
             if ((*jj)>-0.5) { //frequency estimate
                fichier_out<<setw(11)<<(*jj);
                if ((*bb).size()==2) {
                   fichier_out<<setw(8)<<(*bb)[0]<<setw(8)<<(*bb)[1];
                   if ((*bb)[0]<0.0000000001 && (*bb)[1]<0.0000000001) {
                      fichier_out<<"(suspicious CI: heterozygote excess?)";
                   } else if (fabs((*bb)[0]-(*bb)[1])<0.0000000001) {
                      fichier_out<<"(suspicious CI: no heterozygote?)";
                   }
                   fichier_out<<endl;
                }
// add "<<failure" if you want to know why this occurred *in the last instance it did*:
                else fichier_out<<"(No info for CI)"<<endl;
             } else fichier_out<<"No information.\n";
           }
        }
    }
    fichier_out<<"=================================================\n";
    fichier_out<<"\nNormal ending."<<endl;
	fichier_out.close();
    tabFnull.resize(0); // important to flush out anything for next sample
    tabBeta.resize(0);
    bounds.resize(0);
    noR_cout<<"\nNormal ending."<<endl;
    noR_cout<<"Edit the file "<<nom_fic<<" for results"<<endl;
    if (!perf) ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
return 0;
}

//---------------- option 5 -------------------
//------------------------------------------------------------------------------
int descriptif() {
using namespace NS_GP;
int choix=0;

   while (choix != 4) {
            if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("\n");
             printf("\n");
             printf("--->  Allele frequencies, various Fis and gene diversities\n");
             printf("\n");
             printf("      Allele and genotype frequencies per locus and per sample .. 1\n");
             printf("\n");
             printf("      Gene diversities & Fis :\n");
             printf("                                  Using allele identity ......... 2\n");
             printf("                                  Using allele size ............. 3\n");
             printf("\n");
             printf("      Main menu ................................................. 4\n");
             #endif


             if ((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1))
               choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();

             switch (choix) {
                    case 1 : basic_info();return 0;
                    case 2 : Fis_Div(true);return 0;
                    case 3 : Fis_Div(false);return 0;
                    case 4 : return 0;
             }
   }
    return 0;
}




//---------------- option 6 -------------------
//------------------------------------------------------------------------------
int FstIBD() {
using namespace NS_GP;
int choix=0;

if (fichier_genepop->pops.size()==1) {
  RnoR_cerr<<"(!) Only one 'pop' in input file: no information for genetic differentiation."<<endl;
}

while (choix != 8) { // FR 2017/08/18: why not 7 ?
             if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("\n");
             printf("\n");
             printf(" Estimating spatial structure:\n");
             printf("\n");
             printf(" The information considered is :\n");
             printf("      --> Allele identity (F-statistics)\n");
             printf("                For all populations ............ 1\n");
             printf("                For all population pairs ....... 2\n");
             printf("      --> Allele size (Rho-statistics)\n");
             printf("                For all populations ............ 3\n");
             printf("                For all population pairs ....... 4\n");
             printf("\n");
             printf(" Isolation by distance  \n");
             printf("                between individuals ............ 5\n");
             printf("                between groups.................. 6\n");
             printf("\n");
             printf("    Main menu  ................................. 7\n");
             #endif


             if ((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1))
               choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();

             switch (choix) {
                    case 1 : Fstat(true,0);return(0);
                    case 2 : Fstat(true,2);return(0);
                    case 3 : Fstat(false,0);return(0);
                    case 4 : Fstat(false,2);return(0);
                    case 5 : { // individuals
                                    first_repl=true;
                                    isolde_etc(true);
                                    return(0);
                                  }
                    case 6 : { //groups
                                    first_repl=true;
                                    isolde_etc(false);
                                    return(0);
                                  }
                    case 7 : return 0;
             }
             }
    return 0;
}

//---------------- option 7 -------------------
//------------------------------------------------------------------------------
int conversions() {
using namespace NS_GP;
int choix=0;
    while (choix != 5) {
            if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("\n");
             printf("\n");
             printf(" File conversion (diploid data, 2-digits coding only):\n");
             printf("\n");
             printf("      GENEPOP --> FSTAT (F statistics) ........................ 1\n");
             printf("      GENEPOP --> BIOSYS (letter code) ........................ 2\n");
             printf("      GENEPOP --> BIOSYS (number code) ........................ 3\n");
             printf("      GENEPOP --> LINKDOS (D statistics) ...................... 4\n");
             printf("\n");
             printf("      Main menu  .............................................. 5\n");
             #endif


             if ((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1))
               choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();

             switch (choix) {
                    case 1 : conversion(1);return 0;
                    case 2 : conversion(2);return 0;
                    case 3 : conversion(3);return 0;
                    case 4 : conversion(4);return 0;
                    case 5 : return 0;
             }
    }
    return 0;
}

//---------------- option 7 -------------------
//------------------------------------------------------------------------------
int misc() {
using namespace NS_GP;
int choix=0;
while (choix != 5) {
             if(exit_genepop) {return 0;}  //ADD by jimmy
             effacer_ecran();
             afficher_version();
             #ifdef COMPATIBILITYRCPP
                 // Rcpp::R
             #else
             printf("\n");
             printf("\n");
             printf(" Miscellaneous :\n");
             printf("    Null allele: estimates of allele frequencies .......... 1\n");
             printf("    Diploidisation of haploid data ........................ 2\n");
             printf("    Relabeling alleles .................................... 3\n");
             printf("    Conversion to individual data with population names ... 4\n");
             printf("    Conversion to individual data with individual names ... 5\n");
             printf("    Random sampling of haploid genotypes from diploid ones  6\n");
             printf("\n");
             printf("    Main Menu   ........................................... 7\n");
             #endif


             if ((MenuOptions.size()>boucle-1) && (MenuOptions[boucle-1].size()>1))
               choix=MenuOptions[boucle-1][1];
             else choix = controle_choix();

             switch (choix) {
                    case 1 : estimNull();return 0;
                    case 2 : conversion(5);return 0;
                    case 3 : conversion(6);return 0;
                    case 4 : conversion(7);return 0;
                    case 5 : conversion(8);return 0; // added 03/2014 = version 4.2.3
                    case 6 : conversion(9);return 0; // added 03/2014 = version 4.2.3
                    case 7 : return 0;
             }
             }
    return 0;
}

//------------------------------- Menu principal -------------------------------
//------------------------------------------------------------------------------
int menu(){
using namespace NS_GP;

int choix;

while (true) {
    if(exit_genepop) {return 0;}  //ADD by jimmy
    effacer_ecran();
    afficher_version();
    noR_cout<<"Current input file: "<<gp_file<<endl;
    noR_cout<<"Last read at date: "<<fichDATE<<", time: "<<fichTIME<<"\n";
    #ifdef COMPATIBILITYRCPP
        // Rcpp::R
    #else
    printf("\n");
    printf("-------> Change Data ................... C\n");
    printf("\n");
    printf("Testing :\n");
    printf("    Hardy-Weinberg exact tests (several options) ...................... 1\n");
    printf("    Exact tests for genotypic disequilibrium (several options) ........ 2\n");
    printf("    Exact tests for population differentiation (several options) ...... 3\n");
    printf("Estimating:\n");
    printf("    Nm estimates (private allele method) .............................. 4\n");
    printf("    Allele frequencies, various Fis and gene diversities .............. 5\n");
    printf("    Fst & other correlations, isolation by distance (several options).. 6\n");
    printf("Ecumenicism and various utilities:\n");
    printf("    Ecumenicism: file conversion (several options) .................... 7\n");
    printf("    Null alleles and miscellaneous input file utilities ............... 8\n\n");
    printf("QUIT Genepop .......................................................... 9\n\n");
    printf("\n");
    printf("Your choice? : ");
    #endif

//cout<<boucle<<"! ";getchar();
    if (MenuOptions.size()>boucle) {
       choix=MenuOptions[boucle][0];
       boucle++; // attention faire appel a boucle-1 dans les sous options
    } else if (perf) { // perf (donc mode batch) et fin de MenuOptions
       if (MenuOptions.size()==0) {
           #ifdef COMPATIBILITYRCPP
               // Rcpp::R
           #else
                cerr<<"\n(!) Suspect call of performance evaluation without any explicit analysis specified.";
                cerr<<"\n    Check settings. I exit.";
           #endif

            /*** either we had some Performance= call (=> explicitPerf=true and performance_menu is called rather than menu)
                 or we had an 'implicit' performance call without MenuOptions and we arrive here.
            ***/
            genepop_exit(-1, "(!) Suspect call of performance evaluation without any explicit analysis specified.");
            //return -1;
       }
       return 0; // on sort de menu();
    } else { // end of MenuOptions reached and not perf
        if ( ! pauseGP) {
            if (MenuOptions.size()==0 && Mode=="Batch") {
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<endl<<endl<<"(!) Suspect combination of options:  "<<endl;
        		cerr<<"    no apparent MenuOptions, and Mode is Batch "<<endl;
                cerr<<"    (as if Performance setting has been used), "<<endl;
                cerr<<"    but a single input file is considered."<<endl;
                cerr<<"    It may be worth reconsidering the settings."<<endl<<endl;
                #endif

            }
            noR_cout<<"Normal exit; running Mode was "<<Mode<<"."<<endl;
            // mode batch non perf et fin de MenuOptions
            exit_genepop = true;
            return 0; // ADD by jimmy
        } else { // fin de menuOptions en mode non batch
            MenuOptions.resize(0); // vidage par precaution peut etre superflue
            choix =  controle_choix();
        }
    }
    //S?lection des options suivant le choix r?alis?.
//cout<<choix;getchar();

    switch(choix){
                  case 1 : HWexact(); menu(); break;
                  case 2 : LDexact(); menu(); break;
                  case 3 : Diffexact(); menu(); break;
                  case 4 : BartonS86(); menu(); break;
                  case 5 : descriptif(); menu(); break;
                  case 6 : FstIBD(); menu(); break;
                  case 7 : conversions(); menu(); break;
                  case 8 : misc(); menu(); break;
                  case 9 : exit_genepop = true; return 0;//ADD by jimmy //Sortie du programme //ADD by jimmy
                  case 10 : //Changement de donn?es
							remove(fichierIn.c_str());
							delete fichier_genepop;
							ask_new_gp_file(); //exits if Enter...
// eviter le new dans check_... sinon risques de fuites
							fichier_genepop=new CFichier_genepop(gp_file);
							check_gp_file_menu(true); // boucle sur menu possible
                  			return 0;
    }
}
}




//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
using namespace NS_GP;
exit_genepop = false; //ADD by jimmy
#ifdef COMPATIBILITYRCPP
    // Rcpp::R
#else
    cout.setf(std::ios_base::boolalpha); //displays booleans as text
#endif
    fstream fichier;
string cmdlinefilename;
string settingsfilename=getSetting("default_settingsfile");
	effacer_ecran();
	afficher_version();
    if (argc>1) {
// to give inline the name of the file in which command line is written ; or at random...
    	string buf(argv[1]);
    	string::size_type pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
    	string var=buf.substr(0,pos).c_str();
    	if(cmp_nocase(var,"CmdlineFileName")==0) cmdlinefilename=buf.substr(pos+1);
        else cmdlinefilename="cmdline.txt";
        ofstream cmdline(cmdlinefilename.c_str(),ios::out);
    	for (int it=0;it<argc;it++) cmdline<<argv[it]<<endl;
    	cmdline<<endl;
    	cmdline.close();
    // seeks optional SettingsFile in cmdline
    	seeks_settings_file_name(cmdlinefilename,settingsfilename);
    }
    read_settings_file(settingsfilename.c_str()); //cf settings.cpp... READS and SETS...
    if (argc>1) read_settings_file(cmdlinefilename.c_str());
    alea.seed(alea_seed);
    if (perf) {
        /*pauseGP=false;*/ /*removed 30/12/2009. Already the default, we must be able to overcome it.
        indeed some interactive code is present at the beginning of performance_main,
        after which pauseGP is again set to false*/
    	performance_main();
    	return 0;
   	}
    // ELSE	GENEPOP classique
    if (HWfileBool) {HWfileMenu(); if(exit_genepop) {return 0;}} //ADD by jimmy
    else if (isoldeFileBool || multiMigFileBool) {isolde_etc(false);}
    else if (strucFileBool) struc();
	else {
        if (gp_fileInSettingsBool) {
    		fichier.open(fichierIn.c_str(),ios::in);
    	  	if (!fichier.is_open()) fichier.clear();//if no fichier.in, then must be new data
    		else {
			fichier.close(); // slmt pour la clart? du code
			glance_fichier_in(true); //will check it's the same gpFile as in the settings and if so display last read time (seulement affichage, pas d'effet sur la suite)
			}
   			fichier_genepop=new CFichier_genepop(gp_file);
    		check_gp_file_menu(true);
        if(exit_genepop) {return 0;}  //ADD by jimmy
        } else {    // SINON: nom du fichier a lire dans fichier in ou ? donner
    		fichier.open(fichierIn.c_str(),ios::in); //recherche fichier.in
    	  	if (!fichier.is_open()) {//if no fichier.in, then must ask for a file name
    			fichier.clear();
    			ask_new_gp_file();
    			fichier_genepop=new CFichier_genepop(gp_file);
    			check_gp_file_menu(true);
          if(exit_genepop) {return 0;}  //ADD by jimmy
    		}
    		else {
			    fichier.close(); // slmt pour la clart? du code
    			glance_fichier_in(false); // gets gp_file from fichier in
    			fichier_genepop=new CFichier_genepop(gp_file);
    			check_gp_file_menu(false);
          if(exit_genepop) {return 0;}  //ADD by jimmy
    		}
    	}
    }
return 0;
}//main
/****************FIN de MAIN*********************************************************/


int performance_main ( void ) {
	using namespace NS_GP;
	using namespace NS_GP_PERF;
	string zutst;
	if (pauseGP) {// occurs here only if mode has been specified after performance setting
	   if (gp_fileRoot.size()==0 || alwaysAskBool) {
           noR_cout<<"\nName root of the files to analyse? ";
           cin>>gp_fileRoot;
       } // else a root name has been specified
    } else if (gp_fileRoot.size()==0) gp_fileRoot="GP"; //default in batch mode
	if (pauseGP) {
        if (JobMin<0 || alwaysAskBool) {
            noR_cout<< "Index of the first file to analyse? ";
            cin>>JobMin;
        }
    } else if (JobMin<0) JobMin=1; //default in batch mode
	if (pauseGP) {
        if (JobMax<0 || alwaysAskBool) {
            noR_cout<< "\nIndex of the last file to analyse? ";
        	cin>>JobMax;
        }
    } else if (JobMax<0) JobMax=1; //default in batch mode
    alwaysAskBool=false; /** pauseGP not modified: Ask => Default (Default=>Default, Batch=>Batch)**/
	nb_sample=JobMax-JobMin+1;
	effacer_ecran();
	_gotoxy(0,2);
    noR_cout<<"sample: ";
	//***debut de la boucle sur les fichiers (isample membre NS_GP
	for(isample=JobMin;isample<=JobMax;isample++) {
		if(isample==JobMin) {
			first_repl=true; //
			boot_result.open("result.CI");
            if(!boot_result.is_open()){
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Error while opening file result.CI. Exiting."<<endl;
                if (cinGetOnError) cin.get();
                #endif

               genepop_exit(-1, "Error while opening file result.CI.");
            }
        } else {
			first_repl=false; //
			boot_result.open("result.CI",ios::app);
            if(!boot_result.is_open()){
                #ifdef COMPATIBILITYRCPP
                    // Rcpp::R
                #else
                cerr<<"Error while opening file result.CI. Exiting."<<endl;
                if (cinGetOnError) cin.get();
                #endif

               genepop_exit(-1, "Error while opening file result.CI.");
            }
		}
		//cr?? le nom du fichier a analyser
		stringstream stst;
		stst<<gp_fileRoot<<isample;		//nom:Nare->nom:Nare1
		stst>>gp_file;  //de stst vers string
        _gotoxy(9,2);
// affichage subliminal dans le cas du menu(), superpos? au menu qui affiche le fichier pr?c?dent:
        noR_cout<<gp_file<<endl;
		fichier_genepop=new CFichier_genepop(gp_file);
		fichier_genepop->checkName();
		set_eof_check_EOLtype(fichier_genepop->fileName,true); //sanitizes end of file
		fichier_genepop->parseFile();
		fichier_genepop->createFichierIN();
		glance_fichier_in(false);
// ici on doit pouvoir faire un appel ? menu qui utilisera les MenuOptions!
        if (MenuOptions.size()>0) {
            boucle=0;
            menu();
            if(exit_genepop) {return 0;} //ADD by jimmy
        } // MenuOptions pas vide
		else if (explicitPerf) {// options explicites valides sur ligne performance=...
           isolde_etc(indivBool);
        } else {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::R
            #else
            cerr<<"\n(!) Suspect call of performance evaluation without any explicit analysis specified.";
            cerr<<"\n    I exit.";
            if (cinGetOnError) cin.get();
            #endif

            genepop_exit(1, "(!) Suspect call of performance evaluation without any explicit analysis specified.");
        }
        boot_result.close();
	    //if(fclose(boot_result)) printf("result.CI file close error\n");
		delete fichier_genepop;
        noR_cout<<"\a";
	}//fin boucle sur isample
	ZeGenepopSound();
    if (pauseGP) { noR_cout<<"(Return) to continue"<<endl; getchar();}
	return 0;
}
/****************FIN de PERFORMANCE_MAIN*********************************************************/


/*========== Add by Jimmy lopez ==========*/

int mainJimmy(int argc, string argv[]) {
  using namespace NS_GP;
  fstream fichier;
  initialize_for_R();
  string cmdlinefilename;
  string settingsfilename=getSetting("default_settingsfile");
  exit_genepop = false; //ADD by jimmy
	effacer_ecran();
	afficher_version();
  if (argc>1) {
// to give inline the name of the file in which command line is written ; or at random...
    string buf(argv[1]);
  	string::size_type pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
  	string var=buf.substr(0,pos).c_str();
  	if(cmp_nocase(var,"CmdlineFileName")==0) cmdlinefilename=buf.substr(pos+1);
      else cmdlinefilename="cmdline.txt";
      ofstream cmdline(cmdlinefilename.c_str(),ios::out);
  	for (int it=0;it<argc;it++) cmdline<<argv[it]<<endl;
  	cmdline<<endl;
  	cmdline.close();
  // seeks optional SettingsFile in cmdline
  	seeks_settings_file_name(cmdlinefilename,settingsfilename);
  }
  read_settings_file(settingsfilename.c_str()); //cf settings.cpp... READS and SETS...
  // ...(cmdlinefilename.c_str()) finds the gp_file :
  if (argc>1) read_settings_file(cmdlinefilename.c_str());
  alea.seed(alea_seed);
#ifdef COMPATIBILITYRCPP
  // wrong idea. Use std::mt19937 instead
  // set_seed(alea_seed); // set the seed in  the base environment of R
#endif
  if (perf) {
      /*pauseGP=false;*/ /*removed 30/12/2009. Already the default, we must be able to overcome it.
      indeed some interactive code is present at the beginning of performance_main,
      after which pauseGP is again set to false*/
  	performance_main();
  	return 0;
 	}
  // ELSE	GENEPOP classique
  if (HWfileBool) {
    HWfileMenu();
    if(exit_genepop) { clean(false); return 0; }
  } //ADD by jimmy
    else if (isoldeFileBool || multiMigFileBool) {
      isolde_etc(false);
  } else if (strucFileBool) {
    struc();
  } else {
    if (gp_fileInSettingsBool) {
      fichier.open(fichierIn.c_str(),ios::in);
	  	if (!fichier.is_open()) {
	  	  fichier.clear();//if no fichier.in, then must be new data
	  	} else {
	      fichier.close(); // slmt pour la clart? du code
	      glance_fichier_in(true); //will check it's the same gpFile as in the settings and if so display last read time (seulement affichage, pas d'effet sur la suite)
	    }
 		  fichier_genepop=new CFichier_genepop(gp_file);
	  	check_gp_file_menu(true);
    } else {    // SINON: nom du fichier a lire dans fichier in ou ? donner
		  fichier.open(fichierIn.c_str(),ios::in); //recherche fichier.in
	  	if (!fichier.is_open()) {//if no fichier.in, then must ask for a file name
			  fichier.clear();
			  ask_new_gp_file();
			  fichier_genepop=new CFichier_genepop(gp_file);
			  check_gp_file_menu(true);
		  } else {
	      fichier.close(); // slmt pour la clart? du code
			  glance_fichier_in(false); // gets gp_file from fichier in
			  fichier_genepop=new CFichier_genepop(gp_file);
			  check_gp_file_menu(false);
		  }
	  }
    clean(); // all cases with a genepop file return here (whether exit_genepop or not)
    return 0;
  }
  clean(false); // only cases without a genepop file can reach this point
  return 0;
}//mainJimmy

void initialize_for_R() {
  reinitializeGenepopS(); // reset all global variables to their default values
  initializegenepop(); // has cinGetOnError=false;
                       // =>  no need to wrap cinGetOnError use in (! COMPATIBILITYRCPP)
  initializeFest();
  initializeCTtests();
  initializeHWtests();
  initializeMultimig();
  initializeSetting();
  #ifdef COMPATIBILITYRCPP
    initRGenepop();
  #endif
}

void reinitializeGenepopS() {
    // init var
  //version="4.6";
  //settingsfilename="genepop.txt";
  MenuOptions.clear();
  HWfileOptions.clear();
  taille.clear();
  genicProbaTestBool=false;
  alleleNbrTestBool=false;
  geneDivTestBool=false;
  sequenceGeneDivRanks.clear();
  identitySettingsBool=true;
  LDprobaTestBool=false;
  gp_fileInSettingsBool=false;
  perf=false;
  pauseGP=true;
  alwaysAskBool=false;
  HWfileBool=false;
  strucFileBool=false;
  isoldeFileBool=false;
  multiMigFileBool=false;
  estimDiploidBool=true;
  phylipBool=false;
  Brookfield96Bool=false;
  nullIgnoredBool=false;
  NonNullfailuresBool=false;
  gp_file.clear();
  hw_file.clear();
  struc_file.clear();
  isolde_file.clear();
  alea_seed=67144630;
  ABCweight.clear();
  widthCI=0.95;
  outname.clear();
  strcpy(char_tmp, ".TMP");
  strcpy(char_iso, ".ISO");
  strcpy(char_mig, ".MIG");
  first_repl = false;
  dem=1;
  batchlgth=1;
  batchnbr=1;
  strcpy(char_num, ".NUM");
  boucle=0;
  exit_genepop = false;

  //NS_GP::allMax.clear();
  //NS_GP::nom_locus.clear();
  //NS_GP::nom_pop.clear();
  //NS_GP::ploidBool.clear();
  //NS_GP::fichier.close();
  NS_GP::fichDATE.clear();
  NS_GP::fichTIME.clear();
  NS_GP::logdist = false;
  NS_GP::boot_result.close();
  //NS_GP::perfbreak=false;

  NS_GP_PERF::JobMin=-1;
  NS_GP_PERF::isample = 0;
  NS_GP_PERF::JobMax=-1;
  NS_GP_PERF::nb_sample = 0;
  NS_GP_PERF::compteur=0.;
  //fclose(NS_GP_PERF::f_in);
  NS_GP_PERF::gp_fileRoot.clear();

  //NS_Null::coding = 0;
  //NS_Null::tabFnull.clear();
  //NS_Null::tabBeta.clear();
  //NS_Null::bounds.clear();
  //NS_Null::fichier_out.close();
  //NS_Null::genocopy.clear();
  //NS_Null::genocopySum = 0.0;
  //NS_Null::typeMax = 0;
  //NS_Null::failure.clear();
  //NS_Null::gauss_inf = 0.0;
  //NS_Null::gauss_sup = 0.0;
}

void clean(bool fich_genepop) {
  cleanGenepopS();
  if (fich_genepop) cleangenepop();
  cleanFest();
  //cleanCTtests();
  cleanHWtests();
  cleanMultimig();
  cleanSetting();
  #ifdef COMPATIBILITYRCPP
  cleanRGenepop();
  #endif
}

void cleanGenepopS() {
  MenuOptions.clear();
  HWfileOptions.clear();
  taille.clear();
  sequenceGeneDivRanks.clear();
  ABCweight.clear();
  memset(char_iso, 0, sizeof char_iso);
  memset(char_mig, 0, sizeof char_mig);
  memset(char_num, 0, sizeof char_num);

  //NS_GP::allMax.clear();
  //NS_GP::nom_locus.clear();
  //NS_GP::nom_pop.clear();
  //NS_GP::ploidBool.clear();
  //NS_GP::fichier.close();
  NS_GP::boot_result.close();

  //fclose(NS_GP_PERF::f_in);
  NS_GP_PERF::gp_fileRoot.clear();

  //NS_Null::tabFnull.clear();
  //NS_Null::tabBeta.clear();
  //NS_Null::bounds.clear();
  //NS_Null::fichier_out.close();
  //NS_Null::genocopy.clear();

}
