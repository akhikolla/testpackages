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
#ifndef GENEPOPS_H
#define GENEPOPS_H

#include <fstream>
#include <map>
#include <vector>
#include "MersenneTwister.h"
#include "genepop.h"

#define BIGTABLES
const double DRIFT_PREC=0.00000001; //commun ? toutes les MC ? pour U (directionnel)
const double COMP_PREC=0.000000000001; //commun ? toutes les MC ? pour logLR (maximum)
const double HW_PR_PREC=0.000000000000001;
const double HW_U_PREC=.00000000000001;
//const std::string datestring=__DATE__;
//const std::string timestring=__TIME__;

extern std::string getSetting(const std::string which);



extern std::vector<std::map<int,int> >taille;
extern std::vector<std::vector<int> >MenuOptions;
extern std::vector<int> HWfileOptions;
extern std::string outname,gp_file,hw_file,struc_file,isolde_file;
//extern std::string settingsfilename;
extern unsigned long int alea_seed;
extern MTRand alea;
extern size_t dem,batchlgth;
extern size_t batchnbr;
extern char char_iso[];
extern std::vector<double>ABCweight;
extern double widthCI;
extern bool identitySettingsBool,genicProbaTestBool,alleleNbrTestBool,geneDivTestBool,estimDiploidBool,LDprobaTestBool,
       gp_fileInSettingsBool,HWfileBool,strucFileBool,isoldeFileBool,multiMigFileBool,pauseGP,phylipBool,
       alwaysAskBool,perf,GeometryInSettingsBool,Brookfield96Bool,nullIgnoredBool,NonNullfailuresBool;
extern std::vector<int>sequenceGeneDivRanks;

namespace NS_GP {
    //extern std::vector<bool>ploidBool;
}

namespace varForBootstrapGenepop {
  extern std::vector<size_t> idxPloid;
  extern std::string nom_fich_mig;
}

namespace NS_GP_PERF {
  extern int JobMin,JobMax;
	extern std::string gp_fileRoot;
}



int skipln(std::ifstream& st);
//int set_eof_check_EOLtype(); // FR: no def!
int controle_choix();
int option(int choix); // should become obsolete
int check_gp_file_menu(bool verbose);
int ask_new_gp_file();
int menu();
int glance_fichier_in(bool fileCompareBool);
int HWexact();
int isolde_etc(bool indiv);
int FstIBD();
int HWtest(int statindic);
int LDtest();
int LDtables();
int Difftest(int statindic);
void ZeGenepopSound();
void afficher_version();
int set_MC_parameters(bool pasGlobTest);
int performance_main ( void );
int mainJimmy(int argc, std::string argv[]); //ADD by Jimmy
void initialize_for_R();
void clean(bool fich_genepop=true);
void reinitializeGenepopS();
void cleanGenepopS();
int print_p(double pchi, std::ostream& fichier_out, int prec, bool endline);

std::vector<double> estimNullLocPop(const size_t iLoc,const size_t iPop,const bool printBool,
                                    std::map<ssize_t,double>& genocopy, // copy locale en <double> de la map de popGenos;
                                    double& genocopySum, //differs from original sample for delta-perturbed samples !
                                    int& typeMax,
                                    std::ofstream& fichier_out);


int misc();
int conversions();
int descriptif();
int estimNull();
std::vector<double> bootstrapNullAllele(CGenotypes *popGenos,const size_t iLoc, const size_t iPop,
                                          std::map<ssize_t,double>& genocopy, // copy locale en <double> de la map de popGenos;
                                          double& genocopySum, //differs from original sample for delta-perturbed samples !
                                          int& typeMax, // Max allelic type over pops at given locus
                                          std::string& failure,
                                          double& gauss_inf,
                                          double& gauss_sup,
                                          std::ofstream& fichier_out );
double mean_from_creatMat_isolde(std::vector<double> inputwei);
double intercept_from_creatMat_isolde(std::vector<double> inputwei);
double slope_from_creatMat_isolde(std::vector<double> inputwei);
std::vector<double> creatMat_isolde(std::vector<double> inputwei); 
int Fstat(bool identitybool,int Indic);
int Fis_Div(bool identitybool,bool warnOnExit);
int basic_info();
int BartonS86();
int Diffexact();
int LDexact();
int HWfileMenu();
#endif
