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
#include <limits>
#include <map>
#include <sstream>
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
#endif
#include "GenepopS.h"
#include "genepop.h"
#include "myutils.h"
#include "F_est.h"
#include "settings.h"
#include "tools.h"

long int mantelSeed=67144630;
bool GeometryInSettingsBool=false,IsolBDstatInSettingsBool=false;
bool mantelRankBool=false,meanDiffBool=false,singleGeneDiv=false,indivBool,explicitPerf=false; long int mantelPerms=-1;//,mindistInSettingsBool=false;long int mantelPerms=-1;
char enumMCindic=0;

using namespace std;

vector<int>poptypes;  // 'global' from the settings, for Mantel test
string typeSelection="all";
string Mode="Default";
int typeindex1,typeindex2;
double testPointslope=numeric_limits<double>::quiet_NaN();

// As the Rcpp-Exported functions use explicit assignments such as agv[10] = getOptionMantelPermutations(mantelPermutations)
// Then <keyword>=FALSE sets the boolean to true in cases where the <value> was ignored with default = false
// Then we nned the following function:

int evaluateBool(bool &boolean, string buf) { // safe assignment of value `buf' to 'boolean'
  stringstream strstr(buf);
  string locstring;
  strstr>>locstring;
  if(cmp_nocase(locstring,"")==0 || cmp_nocase(locstring,"T")==0 ||
     cmp_nocase(locstring,"True")==0 || cmp_nocase(locstring,"Yes")==0 ||
     cmp_nocase(locstring,"Y")==0)
    boolean=true;
  else if(cmp_nocase(locstring,"F")==0 || cmp_nocase(locstring,"False")==0 ||
          cmp_nocase(locstring,"No")==0 || cmp_nocase(locstring,"N")==0)
    boolean=false;
  else {
#ifdef COMPATIBILITYRCPP
#else
    cout<<"(!) Suspicious specification for a boolean: "<< buf <<endl;
    cout<<"(!) Only \"\", \"T\", \"True\", \"Yes\", \"Y\", \"F\", \"False\", \"No\", and \"N\" are allowed"<<endl;
    cout<<"I exit."<< endl;
    if (cinGetOnError) cin.get();
#endif
    genepop_exit(-1, "Suspicious specification for a boolean.");
  }
  return(0);
}


int seeks_settings_file_name(const string& cmdlinefilename,string& ref_settingsfilename) {
    string buf,var;
    size_t pos;
    ifstream settings(cmdlinefilename.c_str(),ios::in);
    if(!settings.is_open()) {
        noR_cout << "Unable to open file "<<cmdlinefilename<< endl;
        //cerr << "Unable to open file "<<cmdlinefilename<< endl; // modif temp jimmy
    } else do {
		getline(settings,buf);
		if(buf.length()==0) break;
		while((buf[0]==' ')||(buf[0]=='\t')) {buf.erase(0,1);}//vire les blancs initiaux
		pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
		var=buf.substr(0,pos).c_str();
		if(cmp_nocase(var,"SettingsFile")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>ref_settingsfilename;
			settings.close();
			return 0;
		}
		if(cmp_nocase(var,"help")==0) {
          noR_cout<<"** List of settings: \n";
          noR_cout<<" * General options: \n     SettingsFile [Genepop]InputFile MenuOptions  Mode CIcoverage \n";
          noR_cout<<" * Markov chain control:\n     Dememorization BatchLength BatchNumber RandomSeed\n";
          noR_cout<<" * Data manipulation:\n     EstimationPloidy PopTypes popTypeSelection AllelicDistance AlleleSizes\n";
          noR_cout<<" * HW tests:  HWtests DifferentiationTest \n";
          noR_cout<<"   -> Ad hoc file:  HWFile HWfileOptions\n";
          noR_cout<<" * LD tests:  gameticDiseqTest LDTest \n";
          noR_cout<<" * Differentiation tests: DifferentiationTest \n";
          noR_cout<<"   -> Ad hoc file:  strucFile\n";
          noR_cout<<" * Isolation by distance:\n     IsolBDstatistic GeographicScale testPoint\n";
          noR_cout<<"     MinimalDistance MantelPermutations MantelSeed PhylipMatrix\n";
          noR_cout<<"   -> Ad hoc files:  IsolationFile MultiMigFile GeoDistFile\n";
          noR_cout<<" * Null alleles:  NullAlleleMethod\n";
          noR_cout<<" * Performance:  Performance GenepopRootFileName JobMin JobMax\n";
          noR_cout<<" * Misc. information:  Maxima\n";
			genepop_exit(0, "Help");
		}
    } while(!settings.eof());
    settings.close();
return 0;
}

int read_settings_file(const char filename[]) {
string buf,var;
string::size_type pos; //can be long on x86_64 ??
int bidon = 0; // return type of .peek()

    ifstream settings(filename,ios::in);
    if(!settings.is_open()) {
        noR_cout << "Unable to open file "<<filename<< endl;
    	//cerr << "Unable to open file "<<filename<< endl; // modif temp jimmy
    	}
    else {
      noR_cout<<"\nReading settings file "<<filename<<"...\n";

        do {

                getline(settings,buf);
        		if(buf.length()==0) goto nextline; //tres dgrx � changer; parfois il lit une ligne vide la ou on en voit pas.
        		while((buf[0]==' ')||(buf[0]=='\t')) {buf.erase(0,1);}//vire les blancs initiaux
        // va ignorer toutes lignes sans =
        		if ((pos=buf.find('='))==string::npos) goto nextline; //pas de = dans la ligne
        		/*else*/ pos=std::min(pos,std::min(buf.find('\t'),buf.length())); // not sure about the need for this one
                if ((buf[pos])=='=') while (buf[pos-1]==' ') {buf.erase(pos-1,1); pos--;}// vire les blancs avant le =
        		var=buf.substr(0,pos).c_str();
                if (inputCheckBool) {
                    #ifdef COMPATIBILITYRCPP
                        Rcpp::Rcerr<<"|"<<var<<"|\n";
                    #else
                        cerr<<"|"<<var<<"|\n";cin.get();
                    #endif
                }
        //cout<<buf.find('=')<<" "<<pos<<" "<<var<<endl;getchar();
        		if(var.length()==0) goto nextline; // passe les lignes sans =
        		if(cmp_nocase(var,"SettingsFile")==0) {
                    //recognized keyword, but does nothing at this stage!
        			goto nextline;
        		}
        		if(cmp_nocase(var,"Dememorisation")==0 || cmp_nocase(var,"Dememorization")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>dem;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"Maxima")==0) {
               		vector<CLocusGP *>loci;
    		        vector<CPopulation *>pops;
    		        vector<CIndividual *>inds;
                  noR_cout<<"Maximum int value: "<<numeric_limits<int>::max()<<endl;
                  noR_cout<<"Maximum long int value: "<<numeric_limits<long int>::max()<<endl;
                  noR_cout<<"Maximum number of loci: "<<loci.max_size()<<endl;
                  noR_cout<<"Maximum number of populations: "<<pops.max_size()<<endl;
                  noR_cout<<"Maximum number of individuals per population: "<<inds.max_size()<<endl;
//                    cout<<"(Return) to continue."<<endl;
//                    getchar(); // pas de Mode a ce stade
//        			goto nextline;
                    genepop_exit(0, "Error Maxima value");
        		}
        		if(cmp_nocase(var,"InputCheck")==0) {
                    inputCheckBool=true; // may be coming too late though.
        			goto nextline;
        		}
        		if(cmp_nocase(var,"AllelicDistance")==0) { //Switch Type/ Taille
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if(cmp_nocase(locstring,"AlleleSize")==0||cmp_nocase(locstring,"Size")==0) identitySettingsBool=false;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"AlleleSizes")==0) { // indique si digits <> taille
        			vector<int>locvec(0);
        			int value;
        			size_t tailleSize=taille.size(); // 0 � la premi�re ligne AlleleSize... etc.
                    taille.resize(tailleSize+1);
        			string reste=buf.substr(pos+1);
        			stringstream strstr(reste);
        			while (!strstr.eof()) {
                        value=numeric_limits<int>::min();
        				strstr>>value;
        // il stocke les valeurs trouvees dans le vecteur continue � lire jusqu'a eostrstr/chiffre/virgule
        				if (value>numeric_limits<int>::min()) locvec.push_back(value);
                        while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!=',') strstr.get();
        // s'il est � un virgule il doit terminer un vecteur et en commencer un autre
                        if (bidon==',') {
                            if (locvec.size()!=2) {
                                #ifdef COMPATIBILITYRCPP
                                    // Rcpp::Rcerr<<"\nSyntax error in AlleleSizes around numbers ";
                                    //ostream_vector(locvec,Rcpp::Rcerr);
                                    // Rcpp::Rcerr<<"\n I will exit.";
                                #else
                                    cerr<<"\nSyntax error in AlleleSizes around numbers ";
                                    ostream_vector(locvec,cerr);
                                    cerr<<"\n I will exit.";if (cinGetOnError) cin.get();
                                #endif
                               genepop_exit(-1, "Syntax error in AlleleSizes around numbers");
                            } //ELSE
                            if (locvec[0]==0) {
                                #ifdef COMPATIBILITYRCPP
                                    // Rcpp::Rcerr<<"\nTrying to assign a size to allele `0' by the AlleleSizes setting.";
                                    // Rcpp::Rcerr<<"\nThis is not allowed: `0' means no information. I will exit.";if (cinGetOnError) cin.get();
                                #else
                                    cerr<<"\nTrying to assign a size to allele `0' by the AlleleSizes setting.";
                                    cerr<<"\nThis is not allowed: `0' means no information. I will exit.";if (cinGetOnError) cin.get();
                                #endif
                               genepop_exit(-1, "Trying to assign a size to allele `0' by the AlleleSizes setting.");
                            }
        // taille[tailleSize] est une map, locvec[0] va devenir un first et locvec[1] un second dans cette map
                            taille[tailleSize][locvec[0]]=locvec[1];
                            locvec.resize(0);
                            strstr.get(); // et maintenant il faut passer la virgule
                        } //sinon il est soit � la fin de la strstr (=>sortie de la boucle) ou il � un nouvelle value (=>strstr>>value;)
        			}
        // la il est arriv� � la fin du strstr; il a encore un vecteur � empiler, sauf si Toto a termin� sa ligne par une virgule (anticipons)
                    if (locvec.size()==2) taille[tailleSize][locvec[0]]=locvec[1];
        			strstr.clear(); // ??
        			goto nextline;
        		}
        		if(cmp_nocase(var,"BatchLength")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>batchlgth;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"BatchNumber")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>batchnbr;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"DifferentiationTest")==0) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if(cmp_nocase(locstring,"Proba")==0) genicProbaTestBool=true;
        			if(cmp_nocase(locstring,"AlleleNbr")==0) alleleNbrTestBool=true; //private option for EI
        			if(cmp_nocase(locstring,"geneDiv")==0) geneDivTestBool=true; //private option for EI
        			goto nextline;
        		}
                if(cmp_nocase(var,"GeneDivRanks")==0) {
                    sequenceGeneDivRanks.resize(0);
                    int value;
                    string reste=buf.substr(pos+1);
                    stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                        while (!strstr.eof() && !(isdigit(bidon=strstr.peek()))) strstr.get();
                        sequenceGeneDivRanks.push_back(value);
                    }
                    strstr.clear(); // c'est le truc essentiel pour le r�util... snif
                    goto nextline;
                }
        		if(cmp_nocase(var,"EstimationPloidy")==0) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if(cmp_nocase(locstring,"Haploid")==0) {
        			    estimDiploidBool=false;
        			}
        			goto nextline;
        		}
        		if((cmp_nocase(var,"gameticDiseqTest")==0) ||(cmp_nocase(var,"LDTest")==0)) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if(cmp_nocase(locstring,"Proba")==0) LDprobaTestBool=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"MenuOptions")==0) {
        			MenuOptions.resize(0); // safe but probably useless
        			vector<int>locvec(0);
        			int value;
        			string reste=buf.substr(pos+1);
        			stringstream strstr(reste);
        			while (!strstr.eof()) {
        				strstr>>value;
        // il garde la valeur en m�moire et continue � lire jusqu'a eostrstr/chiffre/virgule
        				locvec.push_back(value);
                        while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!=',') strstr.get();
        // s'il est � un virgule il doit terminer un vecteur et en commencer un autre
                        if (bidon==',') {
                            MenuOptions.push_back(locvec);
                            locvec.resize(0);
                            strstr.get(); // et maintenant il faut passer la virgule
                        } //sinon il est soit � la fin de la strstr (=>sortie de la boucle) ou il � un nouvelle value (=>strstr>>value;)
        			}
        // la il est arriv� � la fin du strstr; il a encore un vecteur � empiler, sauf si Toto a termin� sa ligne par une virgule (anticipons)
                    if (locvec.size()>0) MenuOptions.push_back(locvec);
        			goto nextline;
        		}
        		if(cmp_nocase(var,"HWfileOptions")==0) {
        			HWfileOptions.resize(0); // safe but probably useless
        			int value;
        			string reste=buf.substr(pos+1);
        			stringstream strstr(reste);
        			while (!strstr.eof()) {
        				strstr>>value;
        				HWfileOptions.push_back(value);
                        while (!strstr.eof() && (bidon=strstr.peek())!=',') strstr.get();
                        if (bidon==',') strstr.get(); // il faut passer la virgule
                        //sinon il est soit � la fin de la strstr (=>sortie de la boucle) ou il � un nouvelle value (=>strstr>>value;)
        			}
        			goto nextline;
        		}
        		if(cmp_nocase(var,"GenepopInputFile")==0 || cmp_nocase(var,"InputFile")==0) {
//"InputFile" is a bit ambiguous wrt other ad hoc file formats...
                    string locstring=gp_file;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>gp_file;
                    if (gp_fileInSettingsBool || HWfileBool || strucFileBool || isoldeFileBool) {
                      noR_cout<<"Note: new input file "<<gp_file;
                      noR_cout<<"\n replaces a previously declared input file (";
                      noR_cout<<locstring<<hw_file<<struc_file<<isolde_file<<")";
                       if (pauseGP) {
                         noR_cout<<"(Return) to continue"<<endl; getchar();
                       }
                       hw_file=struc_file=isolde_file="";
                       HWfileBool=strucFileBool=isoldeFileBool=false;
                    }
        			gp_fileInSettingsBool=true;
                    perf=false;
        			explicitPerf=false;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"PrivateOptions")==0) {
        			buf=buf.substr(pos+1);
        			pos=0;
                    while (buf.length()>0) {
        			   while (buf.length()>pos && !isalnum(buf[pos])) pos++;
        			   if (buf.length()==pos) goto nextline;
        			   buf=buf.substr(pos);
        			   pos=0;
        			   while (buf.length()>pos && isalnum(buf[pos])) pos++;
        			   //if(cmp_nocase(buf.substr(0,pos),"DiversitiesForDenom")==0) {HaploidIBD=true;} //obsolete !
        			   //if(cmp_nocase(buf.substr(0,pos),"HaploidIBD")==0) {HaploidIBD=true;}
        		       if (inputCheckBool) {
                         noR_cout<<"Parsed: "<<buf.substr(0,pos)<<endl;
                       }
        			   if (buf.length()>pos) {
            				buf=buf.substr(pos+1);
            			   	pos=0;
        			   } else goto nextline;
        			}
        			goto nextline;
        		}
        		if(cmp_nocase(var,"HWFile")==0) {
                    string locstring=hw_file;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>hw_file;
                    if (gp_fileInSettingsBool || HWfileBool || strucFileBool || isoldeFileBool) {
                      noR_cout<<"Note: new input file "<<hw_file;
                      noR_cout<<"\n replaces a previously declared input file (";
                      noR_cout<<gp_file<<locstring<<struc_file<<isolde_file<<")";
                       if (pauseGP) {
                         noR_cout<<"(Return) to continue"<<endl; getchar();
                       }
                       struc_file=gp_file=isolde_file="";
                       strucFileBool=gp_fileInSettingsBool=isoldeFileBool=false;
                    }
        			HWfileBool=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"HWtests")==0) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if(cmp_nocase(locstring,"Enumeration")==0) enumMCindic=1;
        			if(cmp_nocase(locstring,"MCMC")==0) enumMCindic=2;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"NullAlleleMethod")==0) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if((cmp_nocase(locstring,"B96")==0) || (cmp_nocase(locstring,"Brookfield96")==0)) {Brookfield96Bool=true;}
        			if((cmp_nocase(locstring,"C92")==0)) {nullIgnoredBool=true;} //undocum option!
        			if((cmp_nocase(locstring,"ApparentNulls")==0)) {NonNullfailuresBool=true;}
        			goto nextline;
        		}
        		if(cmp_nocase(var,"StrucFile")==0) {
                    string locstring=struc_file;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>struc_file;
                    if (gp_fileInSettingsBool || HWfileBool || strucFileBool || isoldeFileBool) {
                      noR_cout<<"Note: new input file "<<struc_file;
                      noR_cout<<"\n replaces a previously declared input file (";
                      noR_cout<<gp_file<<hw_file<<locstring<<isolde_file<<")";

                       if (pauseGP) {
                         noR_cout<<"(Return) to continue"<<endl; getchar();
                       }
                       hw_file=gp_file=isolde_file="";
                       HWfileBool=gp_fileInSettingsBool=isoldeFileBool=false;
                    }
        			strucFileBool=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"IsolationFile")==0) {
                    string locstring=isolde_file;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>isolde_file;
                    if (gp_fileInSettingsBool || HWfileBool || strucFileBool || isoldeFileBool || multiMigFileBool) {
                      noR_cout<<"Note: new input file "<<isolde_file;
                      noR_cout<<"\n replaces a previously declared input file (";
                      noR_cout<<gp_file<<hw_file<<struc_file<<locstring<<")";
                       if (pauseGP) {
                         noR_cout<<"(Return) to continue"<<endl; getchar();
                       }
                       hw_file=gp_file=struc_file="";
                       HWfileBool=gp_fileInSettingsBool=strucFileBool=multiMigFileBool=false;
                    }
        			isoldeFileBool=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"multiMigFile")==0) {
                    string locstring=isolde_file;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>isolde_file;
                    if (gp_fileInSettingsBool || HWfileBool || strucFileBool || isoldeFileBool || multiMigFileBool) {
                      noR_cout<<"Note: new input file "<<isolde_file;
                      noR_cout<<"\n replaces a previously declared input file (";
                      noR_cout<<gp_file<<hw_file<<struc_file<<locstring<<")";
                       if (pauseGP) {
                         noR_cout<<"(Return) to continue"<<endl; getchar();
                       }
                       hw_file=gp_file=struc_file="";
                       HWfileBool=gp_fileInSettingsBool=strucFileBool=isoldeFileBool=false;
                    }
        			multiMigFileBool=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"geoDistFile")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>geoDistFile;
        			geoDistFromGeoFile=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"MinimalDistance")==0) {
//                    mindistInSettingsBool=true;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>mindist;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"MaximalDistance")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>maxdist;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"PopTypes")==0) { //typically ecotypes to be distinguished in Mantel tests
                    poptypes.resize(0); // discards default values
                    int value;
                    string reste=buf.substr(pos+1);
                    stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
/*strstr.get(); to skip any one-character separator which is necessarily here
is not sufficient is there are whitespaces after the last value. If so, eof is not reached
but no new value is read -> the last value is duplicated */
//(((bidon=strstr.peek())<'0') || (bidon>'9'))
                        while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        poptypes.push_back(value);
                    }
                    strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                    goto nextline;
        		}
        		if(cmp_nocase(var,"PopTypeSelection")==0) {
        		    string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>typeSelection;
               		if(cmp_nocase(typeSelection,"only")==0) {
        			   strstr>>typeindex1;
               		} else if(cmp_nocase(typeSelection,"inter")==0) {
               		    strstr>>typeindex1>>typeindex2;
               		} else if(cmp_nocase(typeSelection,"inter_all_types")==0) {
               		  typeSelection="inter_all_types";
               		} else if(cmp_nocase(typeSelection,"intra_all_types")==0) {
               		  typeSelection="intra_all_types";
               		} else typeSelection="all";
        			goto nextline;
        		}
        		if(cmp_nocase(var,"MantelPermutations")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>mantelPerms;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"MantelRankTest")==0) {
        		  evaluateBool(mantelRankBool,buf.substr(pos+1));
        			goto nextline;
        		}
                if(cmp_nocase(var,"MeanDifferentiationTest")==0) {
                  evaluateBool(meanDiffBool,buf.substr(pos+1));
                    goto nextline;
                }
                if(cmp_nocase(var,"MantelSeed")==0) { /*not documented in GenepopS documentation*/
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>mantelSeed;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"PhylipMatrix")==0) {
        		  evaluateBool(phylipBool,buf.substr(pos+1));
        			goto nextline;
        		}
        		if(cmp_nocase(var,"RandomSeed")==0) { //tout sauf Mantel
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>alea_seed;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"Mode")==0) {
        /* Pause determines only two correct contexts for getchar():
               cerr<< error message + if(cinGetOnError) getchar() + exit // only false in Batch mode
        and
               cout<< some info + if(pauseGP) getchar() + execution continues //false in Batch and Batchdebug=cinGetOnError mode
        */
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			// Mode was introduced on 10/12/2011 and is not (yet) systematically used
        			if((cmp_nocase(locstring,"Batch")==0)) {Mode="Batch";pauseGP=false;alwaysAskBool=false;cinGetOnError=false;}
        			if((cmp_nocase(locstring,"BatchDebug")==0) || (cmp_nocase(locstring,"PauseOnError")==0)) {
        			    Mode="BatchDebug";pauseGP=false;alwaysAskBool=false;cinGetOnError=true;
                    } //private
        			if((cmp_nocase(locstring,"Ask")==0)) {Mode="Ask";pauseGP=true;alwaysAskBool=true;cinGetOnError=true;}
        // else defaults are alwaysAskBool = false mais pauseGP =true
        			if((cmp_nocase(locstring,"Default")==0)) {Mode="Default";pauseGP=true;alwaysAskBool=false;cinGetOnError=true;}
        			goto nextline;
        		}
        		if(cmp_nocase(var,"IsolBDstatistic")==0 || cmp_nocase(var,"IsolationStatistic")==0) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if((cmp_nocase(locstring,"a")==0)) {_a_stat=true;_e_stat=false;singleGeneDiv=false;}
        			if((cmp_nocase(locstring,"e")==0)) {_a_stat=false;_e_stat=true;singleGeneDiv=false;}
        			if((cmp_nocase(locstring,"singleGeneDiv")==0)) {_e_stat=false;_e_stat=false;singleGeneDiv=true;}
        			IsolBDstatInSettingsBool=true;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"GeographicScale")==0 || cmp_nocase(var,"Geometry")==0) {
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if((cmp_nocase(locstring,"Linear")==0) || (cmp_nocase(locstring,"1D")==0)) _logdist="identity";
        			if((cmp_nocase(locstring,"Planar")==0) || (cmp_nocase(locstring,"2D")==0) || (cmp_nocase(locstring,"Log")==0)) _logdist="log";
        //			GeometryInSettingsBool=true;
           			goto nextline;
        		}
        		if(cmp_nocase(var,"CIcoverage")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>widthCI;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"Performance")==0) {
        			explicitPerf=false;
        			string locstring;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>locstring;
        			if ((cmp_nocase(locstring,"aLinear")==0) || (cmp_nocase(locstring,"a1D")==0))
                       {indivBool=true;_a_stat=true;_e_stat=false;_logdist="identity";explicitPerf=true;}
        			else if ((cmp_nocase(locstring,"eLinear")==0) || (cmp_nocase(locstring,"e1D")==0))
                       {indivBool=true;_a_stat=false;_e_stat=true;_logdist="identity";explicitPerf=true;}
        			else if ((cmp_nocase(locstring,"aPlanar")==0) || (cmp_nocase(locstring,"a2D")==0))
                       {indivBool=true;_a_stat=true;_e_stat=false;_logdist="log";explicitPerf=true;}
        			else if ((cmp_nocase(locstring,"ePlanar")==0) || (cmp_nocase(locstring,"e2D")==0))
                       {indivBool=true;_a_stat=false;_e_stat=true;_logdist="log";explicitPerf=true;}
        			else if ((cmp_nocase(locstring,"FLinear")==0) || (cmp_nocase(locstring,"F1D")==0))
                       {indivBool=false;_logdist="identity";explicitPerf=true;}
        			else if ((cmp_nocase(locstring,"FPlanar")==0) || (cmp_nocase(locstring,"F2D")==0))
                       {indivBool=false;_logdist="log";explicitPerf=true;}
        			else if ((cmp_nocase(locstring,"")==0))
                       {}
                    else {
                        #ifdef COMPATIBILITYRCPP
                            // Rcpp::Rcerr<<"\nUnknown performance option "<<locstring;
                            // Rcpp::Rcerr<<"\n I will exit.";
                        #else
                            cerr<<"\nUnknown performance option "<<locstring;
                            cerr<<"\n I will exit.";if (cinGetOnError) cin.get();
                        #endif

                       genepop_exit(-1, "Unknown performance option");
                    } //ELSE
        			perf=true;
                    Mode="Batch";
        			pauseGP=false; // by default in batch mode
                    alwaysAskBool=false;
                    // note that cinGetOnError may still be true (only modified by Mode setting)
        			goto nextline;
        		}
        		if(cmp_nocase(var,"GenepopRootFileName")==0 || cmp_nocase(var,"GenepopRootFile")==0) { //02/2011 sufficient to activate 'performance' mode
        			perf=true;
                    Mode="Batch";
        			pauseGP=false; // by default in batch mode
                    alwaysAskBool=false;
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>NS_GP_PERF::gp_fileRoot;
                    // note that cinGetOnError may still be true (only modified by Mode setting)
                    goto nextline;
        		}
        		if(cmp_nocase(var,"JobMin")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>NS_GP_PERF::JobMin;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"JobMax")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>NS_GP_PERF::JobMax;
        			goto nextline;
        		}
        		if(cmp_nocase(var,"testPoint")==0) {
        			stringstream strstr(buf.substr(pos+1));
        			strstr>>testPointslope;
        			goto nextline;
        		}
        // ligne avec = (cf test au debut), keyword inconnu et pas commented out
                if (! (var[0]=='%' || var[0]=='#' || var[0]=='/' ) && cmp_nocase(var,"cmdlinefilename")!=0) {
                  noR_cout<<"(!) Unknown keyword \""<<var<<"\"\n";
                   if (pauseGP) {
                     noR_cout<<"(Return) to continue"<<endl; getchar();
                   }
                }
    nextline: ;	// the pleasure of sin :-)
        } while(!settings.eof());
    }
    settings.close();
return 0;
}

void initializeSetting() {
  GeometryInSettingsBool=false;
  IsolBDstatInSettingsBool=false;
  mantelRankBool=false;
  singleGeneDiv=false;
  indivBool=false;
  explicitPerf=false;
  //mantelPerms=-1;
  enumMCindic=0;
  poptypes.clear();
  typeSelection="all";
  //Mode="Default";
  typeindex1 = 0;
  typeindex2 = 0;
  testPointslope=numeric_limits<double>::quiet_NaN();
}

void cleanSetting() {
  poptypes.clear();
}
