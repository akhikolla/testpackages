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

#include "RGenepop.h"
#include "GenepopS.h"
#include "tools.h"

#include <sstream>
#include <unistd.h>


const bool debugMode = false;
static long RrandomSeed = 67144630;
static long RmantelSeed = 67144630;

/*std::ofstream outPut;
std::streambuf *coutbuf;
FILE *stream;*/

/* Tools */

#ifdef COMPATIBILITYRCPP
// function not used !
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}
#endif

void initRGenepop() {
 /* outPut.open("out.txt");
  coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(outPut.rdbuf());
  stream = freopen("out2.txt", "w", stdout);*/
}

void cleanRGenepop() {
 /* outPut.close();
  std::ofstream ofs;
  ofs.open("out.txt", std::ofstream::out | std::ofstream::trunc);
  ofs.close();
  ofs.open("out2.txt", std::ofstream::out | std::ofstream::trunc);
  ofs.close();*/
}

// //' Reset Random Seed
// // [[Rcpp::export]]
// void resetSeed() {
//   RrandomSeed = 67144630;
//   RmantelSeed = 67144630;
// }

//'@name genepop-utils
//'@title Programming utilities 
//'@description \code{getVersion} returns the version number of the C++ code (the same number that identifies the C++ executable). \code{set_restriction(TRUE)} sets the maximum number of populations and of loci to 300.
//'@param set logical: whether to set restrictions on number of populations and of loci
// [[Rcpp::export]]
std::string getVersion() { 
  return(getSetting("version")); 
}


std::string getNameProg() {
  return "Genepop";
}

void printCmd(int size, std::string argv[]) {
  if(debugMode) {
    for(int i=0;i<size;i++) Ronly_cout << argv[i] << " ";
  }
}

int getNumberLineFile(std::string file) {
  std::ifstream infile(file.c_str());
  int numberLine = 0;
  std::string line;
  while (std::getline(infile, line)) {
      numberLine++;
  }
  infile.close();
  return numberLine;
}


long getRandomSeed() {
  return RrandomSeed;
}

long getMantelSeed() {
  return RmantelSeed;
}

//' Set random generator seed (except for Mantel test)
//' @usage setRandomSeed(seed)
//' @param seed integer: the new seed
// [[Rcpp::export]]
void setRandomSeed(long seed) {
  RrandomSeed = seed;
}

//' Set random generator seed for Mantel test
//' @usage setMantelSeed(seed)
//' @param seed integer: the new seed
// [[Rcpp::export]]
void setMantelSeed(long seed) {
  RmantelSeed = seed;
}

/* Option Genepop */

std::string getOptionInputFile(std::string inputFile) {
  return "GenepopInputFile="+inputFile;
}

std::string getOptionMenu(std::string menuOptions) {
  return "MenuOptions="+menuOptions;
}

std::string getOptionEnumeration(bool enumeration) {
  if(enumeration) { return "HWtests=enumeration";}
  else { return "HWtests=MCMC";}
}

std::string getOptionDememorisation(int dememorization) {
  std::ostringstream stream;
  stream.clear();
  stream << "Dememorisation=" << dememorization;
  return stream.str();
}

std::string getOptionBatchNumber(int batches) {
  std::ostringstream stream;
  stream.clear();
  stream << "BatchNumber=" << batches;
  return stream.str();
}

std::string getOptionBatchLength(int batches) {
  std::ostringstream stream;
  stream.clear();
  stream << "BatchLength=" << batches;
  return stream.str();
}

std::string getOptionNullAlleleMethod(std::string method) {
  return "NullAlleleMethod="+method;
}

std::string getOptionICoverage(double CIcoverage) {
  std::ostringstream stream;
  stream.clear();
  stream << "CIcoverage=" << CIcoverage;
  return stream.str();
}

std::string getOptionEstimationPloidy(std::string dataType) {
  return "EstimationPloidy="+dataType;
}

std::string getOptionStructFile(std::string inputFile) {
  return "StrucFile="+inputFile;
}

std::string getOptionIsolationStatistic(std::string statistic) {
  return "IsolationStatistic="+statistic;
}

std::string getOptionGeographicScale(std::string scale) {
  return "Geometry="+scale;
}


std::string getOptionTestPoint(double testPoint) {
  std::ostringstream stream;
  stream.clear();
  stream << "testPoint=" << testPoint;
  return stream.str();
}

std::string getOptionMinimalDistance(double minDist) {
  std::ostringstream stream;
  stream.clear();
  stream << "MinimalDistance=" << minDist;
  return stream.str();
}
std::string getOptionMaximalDistance(double maxDist) {
  std::ostringstream stream;
  stream.clear();
  stream << "MaximalDistance=" << maxDist;
  return stream.str();
}

std::string getOptionMantelPermutations(int permutations) {
  std::ostringstream stream;
  stream.clear();
  stream << "MantelPermutations=" << permutations;
  return stream.str();
}

std::string getOptionMantelRankTest(bool mantelRank) {
  std::ostringstream stream;
  stream.clear();
  if(mantelRank){
      stream << "MantelRankTest=TRUE";
  }
  else {
      stream << "MantelRankTest=FALSE";
  }

  return stream.str();
}

std::string getOptionIsolationFile(std::string inputFile) {
  return "IsolationFile="+inputFile;
}

std::string getOptionRandomSeed(long seed) {
  std::ostringstream stream;
  stream.clear();
  stream << "RandomSeed=" << seed;
  return stream.str();
}

std::string getOptionMantelSeed(long seed) {
  std::ostringstream stream;
  stream.clear();
  stream << "MantelSeed=" << seed;
  return stream.str();
}

std::string getOptionAllelicDistance(std::string distance) {
  return "AllelicDistance="+distance;
}

std::string getOptionAlleleSizes(std::string sizes) {
  return "AlleleSizes="+sizes;
}

std::string getOptionHWFile(std::string inputFile)  {
  return "HWFile="+inputFile;
}

std::string getOptionHWFileMenu(int option) {
  std::ostringstream stream;
  stream.clear();
  stream << "HWfileOptions=" << option;
  return stream.str();
}

std::string getOptionPopTypes(std::string popTypes) {
  return "PopTypes="+popTypes;
}

std::string getOptionPopTypeSelection(std::string popTypesSelection) {
  return "PopTypeSelection="+popTypesSelection;
}

std::string getOptionMultiMigFile(std::string inputFile) {
  return "MultiMigFile="+inputFile;
}

std::string getOptionGeoDistFile(std::string inputFile) {
  return "geoDistFile="+inputFile;
}

std::string getOptionSettingsFile(std::string settingsFile) {
  return "settingsFile="+settingsFile;
}

std::string getOptionModeBatch() {
  return "Mode=Batch";
}

/* Extension output file */

std::string getOutPutFileMenu_1_1(std::string inputFile) {
  return inputFile+".D";
}

std::string getOutPutFileMenu_1_2(std::string inputFile) {
  return inputFile+".E";
}

std::string getOutPutFileMenu_1_3(std::string inputFile) {
  return inputFile+".P";
}

std::string getOutPutFileMenu_1_4(std::string inputFile) {
  return inputFile+".DG";
}

std::string getOutPutFileMenu_1_5(std::string inputFile) {
  return inputFile+".EG";
}

std::string getOutPutFileMenu_2_1(std::string inputFile) {
  return inputFile+".DIS";
}

std::string getOutPutFileMenu_2_2(std::string inputFile) {
  return inputFile+".TAB";
}

std::string getOutPutFileMenu_3_1(std::string inputFile) {
  return inputFile+".GE";
}

std::string getOutPutFileMenu_3_2(std::string inputFile) {
  return inputFile+".GE2";
}

std::string getOutPutFileMenu_3_3(std::string inputFile) {
  return inputFile+".G";
}

std::string getOutPutFileMenu_3_4(std::string inputFile) {
  return inputFile+".2G2";
}

std::string getOutPutFileMenu_4_1(std::string inputFile) {
  return inputFile+".PRI";
}

std::string getOutPutFileMenu_5_1(std::string inputFile) {
  return inputFile+".INF";
}

std::string getOutPutFileMenu_5_2(std::string inputFile) {
  return inputFile+".DIV";
}

std::string getOutPutFileMenu_5_3(std::string inputFile) {
  return inputFile+".MSD";
}

std::string getOutPutFileMenu_6_1(std::string inputFile) {
  return inputFile+".FST";
}

std::string getOutPutFileMenu_6_2(std::string inputFile) {
  return inputFile+".ST2";
}

std::string getOutPutFileMenu_6_2_b(std::string inputFile) {
  return inputFile+".MIG";
}

std::string getOutPutFileMenu_6_3(std::string inputFile) {
  return inputFile+".RHO";
}

std::string getOutPutFileMenu_6_4(std::string inputFile) {
  return inputFile+".ST2";
}

std::string getOutPutFileMenu_6_4_b(std::string inputFile) {
  return inputFile+".MIG";
}

std::string getOutPutFileMenu_6_5(std::string inputFile) {
  return inputFile+".ISO";
}

std::string getOutPutFileMenu_6_5_b(std::string inputFile) {
  return inputFile+".GRA";
}

std::string getOutPutFileMenu_6_5_c(std::string inputFile) {
  return inputFile+".MIG";
}

std::string getOutPutFileMenu_6_6(std::string inputFile) {
  return inputFile+".ISO";
}

std::string getOutPutFileMenu_6_6_b(std::string inputFile) {
  return inputFile+".GRA";
}

std::string getOutPutFileMenu_6_6_c(std::string inputFile) {
  return inputFile+".MIG";
}

std::string getOutPutFileMenu_7_1(std::string inputFile) {
  return inputFile+".DAT";
}

std::string getOutPutFileMenu_7_2(std::string inputFile) {
  return inputFile+".BIO";
}

std::string getOutPutFileMenu_7_3(std::string inputFile) {
  return inputFile+".BIO";
}

std::string getOutPutFileMenu_7_4(std::string inputFile) {
  return inputFile+".LKD";
}

std::string getOutPutFileMenu_8_1(std::string inputFile) {
  return inputFile+".NUL";
}

std::string getOutPutFileMenu_8_2(std::string inputFile) {
  std::string pathFile = inputFile.substr(0, inputFile.find_last_of("/\\")+1); // add by jimmy
  std::string nameFile = inputFile.substr(inputFile.find_last_of("/\\")+1); // add by jimmy
  return pathFile+"D"+nameFile;
}

std::string getOutPutFileMenu_8_3(std::string inputFile) {
  std::string pathFile = inputFile.substr(0, inputFile.find_last_of("/\\")+1); // add by jimmy
  std::string nameFile = inputFile.substr(inputFile.find_last_of("/\\")+1); // add by jimmy
  return pathFile+"N"+nameFile;
}

std::string getOutPutFileMenu_8_4(std::string inputFile) {
  std::string pathFile = inputFile.substr(0, inputFile.find_last_of("/\\")+1); // add by jimmy
  std::string nameFile = inputFile.substr(inputFile.find_last_of("/\\")+1); // add by jimmy
  return pathFile+"I"+nameFile;
}

std::string getOutPutFileMenu_8_5(std::string inputFile) {
  std::string pathFile = inputFile.substr(0, inputFile.find_last_of("/\\")+1); // add by jimmy
  std::string nameFile = inputFile.substr(inputFile.find_last_of("/\\")+1); // add by jimmy
  return pathFile+"I"+nameFile;
}

std::string getOutPutFileMenu_8_6(std::string inputFile) {
  std::string pathFile = inputFile.substr(0, inputFile.find_last_of("/\\")+1); // add by jimmy
  std::string nameFile = inputFile.substr(inputFile.find_last_of("/\\")+1); // add by jimmy
  return pathFile+"H"+nameFile;
}


/* Exact Hardy-Weinberg tests */

// [[Rcpp::export]]
std::string RHWEachLocusEachPopulationHD(std::string inputFile, std::string outputFile, bool enumeration, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 9;
  string agv[9];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:1");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionEnumeration(enumeration);
  agv[5] = getOptionBatchNumber(batches);
  agv[6] = getOptionBatchLength(iterations);
  agv[7] = getOptionRandomSeed(getRandomSeed());
  agv[8] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if( ! outputFile.empty()) {
    rename(getOutPutFileMenu_1_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_1(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWEachLocusEachPopulationHE(std::string inputFile, std::string outputFile, bool enumeration, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 9;
  string agv[9];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:2");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionEnumeration(enumeration);
  agv[5] = getOptionBatchNumber(batches);
  agv[6] = getOptionBatchLength(iterations);
  agv[7] = getOptionRandomSeed(getRandomSeed());
  agv[8] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_2(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWEachLocusEachPopulationProbability(std::string inputFile, std::string outputFile, bool enumeration, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 9;
  string agv[9];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:3");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionEnumeration(enumeration);
  agv[5] = getOptionBatchNumber(batches);
  agv[6] = getOptionBatchLength(iterations);
  agv[7] = getOptionRandomSeed(getRandomSeed());
  agv[8] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_3(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RHWGlobalHD(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:4");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_4(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_4(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWGlobalHE(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:5");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_5(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_5(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWEachLocusEachPopulationHDWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:1");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_1(inputFile).c_str()); 
  }

}

// [[Rcpp::export]]
std::string RHWEachLocusEachPopulationHEWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:2");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if( ! outputFile.empty()) {
    rename(getOutPutFileMenu_1_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_2(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWEachLocusEachPopulationProbabilityWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:3");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_3(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWGlobalHDWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:4");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_4(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_4(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RHWGlobalHEWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("1:5");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_1_5(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_1_5(inputFile).c_str()); 
  }
}

/* Exact Genotypic Disequilibrium tests */

// [[Rcpp::export]]
std::string RGDEachPairLociEachPopulation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("2:1");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_2_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_2_1(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RGDGenotypicContingency(std::string inputFile, std::string outputFile ) {
  using namespace std;
  
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("2:2");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_2_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_2_2(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RGDEachPairLociEachPopulationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("2:1");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_2_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_2_1(inputFile).c_str()); 
  }
}

/* Exact Population Differentiation tests */

// [[Rcpp::export]]
std::string RPDGenicAllPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:1");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_1(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RPDGenicAllPairPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:2");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_2(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RPDGenotypicAllPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:3");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_3(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RPDGenotypicAllPairPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations ) {
  using namespace std;
  
  int agc = 8;
  string agv[8];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:4");
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionBatchNumber(batches);
  agv[5] = getOptionBatchLength(iterations);
  agv[6] = getOptionRandomSeed(getRandomSeed());
  agv[7] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_4(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_4(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RPDGenicAllPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:1");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_1(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RPDGenicAllPairPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:2");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_2(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RPDGenotypicAllPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:3");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_3(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RPDGenotypicAllPairPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("3:4");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_3_4(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_3_4(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RAnalyzingSingleContingencyTable(std::string inputFile, int dememorization, int batches, int iterations) {
  using namespace std;
  int agc = 7;
  string agv[7];

  agv[0] = getNameProg();
  agv[1] = getOptionStructFile(inputFile);
  agv[2] = getOptionDememorisation(dememorization);
  agv[3] = getOptionBatchNumber(batches);
  agv[4] = getOptionBatchLength(iterations);
  agv[5] = getOptionRandomSeed(getRandomSeed());
  agv[6] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  
  return(inputFile.c_str()); 
}

// [[Rcpp::export]]
std::string RAnalyzingSingleContingencyTableWithSettingsFile(std::string inputFile, std::string settingsFile) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+3;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 2;
  agv[0] = getNameProg();
  agv[1] = getOptionStructFile(inputFile);
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  return(inputFile.c_str()); 
}

/* Nm estimates */

// [[Rcpp::export]]
std::string RNmEstimates(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("4");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_4_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_4_1(inputFile).c_str()); 
  }
}

/* Descriptif */

// [[Rcpp::export]]
std::string RDescriptifAlleleAndGenotypeFrequenciesPerLocusPerSample(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("5:1");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_5_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_5_1(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RDescriptifGeneDiversitiesAndFisUsingAlleleIdentity(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("5:2");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_5_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_5_2(inputFile).c_str()); 
  }
}
// MenuOptions=5:2 EstimationPloidy=Diploid

// [[Rcpp::export]]
std::string RDescriptifGeneDiversitiesAndFisUsingAlleleSize(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("5:3");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_5_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_5_3(inputFile).c_str()); 
  }
  
}

/* FstIBD */

// [[Rcpp::export]]
std::string REstimatingSpatialStructureAlleleIdentyAllPopulations(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("6:1");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_6_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_6_1(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string REstimatingSpatialStructureAlleleIdentyAllPopulationsPairs(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("6:2");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);
  mainJimmy(agc, agv);

  if(!outputFile.empty()) {

    rename(getOutPutFileMenu_6_2(inputFile).c_str(), outputFile.c_str());
    rename(getOutPutFileMenu_6_2_b(inputFile).c_str(), getOutPutFileMenu_6_2_b(outputFile).c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_6_2(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string REstimatingSpatialStructureAlleleSizeAllPopulations(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("6:3");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_6_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_6_3(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string REstimatingSpatialStructureAlleleSizeAllPopulationsPairs(std::string inputFile, std::string outputFile, std::string dataType ) {
  using namespace std;
  int agc = 6;
  string agv[6];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("6:4");
  agv[3] = getOptionEstimationPloidy(dataType);
  agv[4] = getOptionRandomSeed(getRandomSeed());
  agv[5] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_6_4(inputFile).c_str(), outputFile.c_str());
    //Rcpp::Rcerr <<"LA"<<getOutPutFileMenu_6_4_b(inputFile).c_str()<<" "<<getOutPutFileMenu_6_4_b(outputFile).c_str() ;
    rename(getOutPutFileMenu_6_4_b(inputFile).c_str(), getOutPutFileMenu_6_4_b(outputFile).c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_6_4(inputFile).c_str()); 
  }
  
}


 // [[Rcpp::export]]
 std::string RIsolationByDistanceBetweenIndividuals (std::string inputFile, std::string outputFile, std::string dataType, 
                                                     std::string statistic, std::string geographicScale,
                                                     double CIcoverage, double testPoint, double minimalDistance, 
                                                     double maximalDistance, int mantelPermutations, bool mantelRankTest ) {
   using namespace std;
   
    int agc = 14;
    string agv[14];

    agv[0] = getNameProg();
    agv[1] = getOptionInputFile(inputFile);
    agv[2] = getOptionMenu("6:5");
    agv[3] = getOptionEstimationPloidy(dataType);
    agv[4] = getOptionIsolationStatistic(statistic);
    agv[5] = getOptionGeographicScale(geographicScale);
    agv[6] = getOptionICoverage(CIcoverage);
    agv[7] = getOptionTestPoint(testPoint);
    agv[8] = getOptionMinimalDistance(minimalDistance);
    agv[9] = getOptionMaximalDistance(maximalDistance);
    agv[10] = getOptionMantelPermutations(mantelPermutations);
    agv[11] = getOptionMantelRankTest(mantelRankTest);
    agv[12] = getOptionMantelSeed(getMantelSeed());
    agv[13] = getOptionModeBatch();

    printCmd(agc, agv);
    
    mainJimmy(agc, agv);

    //system("rm -f LOCUS*"); // return value ignored

    if(!outputFile.empty()) {
      rename(getOutPutFileMenu_6_5(inputFile).c_str(), outputFile.c_str());
      rename(getOutPutFileMenu_6_5_b(inputFile).c_str(), getOutPutFileMenu_6_5_b(outputFile).c_str());
      rename(getOutPutFileMenu_6_5_c(inputFile).c_str(), getOutPutFileMenu_6_5_c(outputFile).c_str());
      return(outputFile.c_str()); 
    } else {
      return(getOutPutFileMenu_6_5(inputFile).c_str()); 
    }

}

  // [[Rcpp::export]]
std::string RIsolationByDistanceBetweenGroups(std::string inputFile, std::string outputFile, 
                                              std::string dataType, std::string statistic, std::string geographicScale,
                                              double CIcoverage, double testPoint, double minimalDistance, 
                                              double maximalDistance, int mantelPermutations, bool mantelRankTest ) {
  using namespace std;
  int agc = 14;
    string agv[14];

    agv[0] = getNameProg();
    agv[1] = getOptionInputFile(inputFile);
    agv[2] = getOptionMenu("6:6");
    agv[3] = getOptionEstimationPloidy(dataType);
    agv[4] = getOptionIsolationStatistic(statistic);
    agv[5] = getOptionGeographicScale(geographicScale);
    agv[6] = getOptionICoverage(CIcoverage);
    agv[7] = getOptionTestPoint(testPoint);
    agv[8] = getOptionMinimalDistance(minimalDistance);
    agv[9] = getOptionMaximalDistance(maximalDistance);
    agv[10] = getOptionMantelPermutations(mantelPermutations);
    agv[11] = getOptionMantelRankTest(mantelRankTest);
    agv[12] = getOptionMantelSeed(getMantelSeed());
    agv[13] = getOptionModeBatch();

    printCmd(agc, agv);

    mainJimmy(agc, agv);

    //system("rm -f LOCUS*"); // return value ignored

    if(!outputFile.empty()) {
      rename(getOutPutFileMenu_6_6(inputFile).c_str(), outputFile.c_str());
      rename(getOutPutFileMenu_6_6_b(inputFile).c_str(), getOutPutFileMenu_6_6_b(outputFile).c_str());
      rename(getOutPutFileMenu_6_6_c(inputFile).c_str(), getOutPutFileMenu_6_6_c(outputFile).c_str());
      return(outputFile.c_str()); 
    } else {
      return(getOutPutFileMenu_6_6(inputFile).c_str()); 
    }
    
}

// [[Rcpp::export]]
std::string RIsolationByDistanceBetweenIndividualsWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;
  
  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("6:5");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  
  delete[] agv;

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_6_5(inputFile).c_str(), outputFile.c_str());
    rename(getOutPutFileMenu_6_5_b(inputFile).c_str(), getOutPutFileMenu_6_5_b(outputFile).c_str());
    rename(getOutPutFileMenu_6_5_c(inputFile).c_str(), getOutPutFileMenu_6_5_b(outputFile).c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_6_5(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string RIsolationByDistanceBetweenGroupsWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile ) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("6:6");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_6_6(inputFile).c_str(), outputFile.c_str());
    rename(getOutPutFileMenu_6_6_b(inputFile).c_str(), getOutPutFileMenu_6_6_b(outputFile).c_str());
    rename(getOutPutFileMenu_6_6_c(inputFile).c_str(), getOutPutFileMenu_6_6_b(outputFile).c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_6_6(inputFile).c_str()); 
  }
}

/* Ecumenicism */

// [[Rcpp::export]]
std::string REcumenicismFstat(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("7:1");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_7_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_7_1(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string REcumenicismBiosysLetter(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("7:2");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_7_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_7_2(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string REcumenicismBiosysNumber(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("7:3");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_7_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_7_3(inputFile).c_str()); 
  }
}

// [[Rcpp::export]]
std::string REcumenicismLinkdos(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("7:4");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_7_4(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_7_4(inputFile).c_str()); 
  }
  
}

/* Misc */

// [[Rcpp::export]]
std::string RNullAlleleEstimateAlleleFrequencies(std::string inputFile, std::string outputFile, std::string nullAlleleMethod, double CIcoverage ) {
  using namespace std;
  int agc = 7;
  string agv[7];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:1");
  agv[3] = getOptionNullAlleleMethod(nullAlleleMethod);
  agv[4] = getOptionICoverage(CIcoverage);
  agv[5] = getOptionRandomSeed(getRandomSeed());
  agv[6] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_1(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RDiploidisationHaploidData(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:2");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_2(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_2(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RRelabelingAlleles(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:3");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_3(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_3(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RConversionToIndividualDataWithPopulationNames(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:4");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_4(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_4(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RConversionToIndividualDataWithIndividualNames(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:5");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_5(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_5(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RRandomSamplingOfHaploidGenotypesFromDiploidOnes(std::string inputFile, std::string outputFile ) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:6");
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_6(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_6(inputFile).c_str()); 
  }
  
}

// [[Rcpp::export]]
std::string RNullAlleleEstimateAlleleFrequenciesWithSettingsFile(std::string inputFile, 
                                                                 std::string outputFile, std::string settingsFile) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionInputFile(inputFile);
  agv[2] = getOptionMenu("8:1");
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  if(!outputFile.empty()) {
    rename(getOutPutFileMenu_8_1(inputFile).c_str(), outputFile.c_str());
    return(outputFile.c_str()); 
  } else {
    return(getOutPutFileMenu_8_1(inputFile).c_str()); 
  }
  
}

/* MOAO */

// [[Rcpp::export]]
std::string RHWtableHD(std::string inputFile, bool enumeration, int dememorization, int batches, int iterations) {
  using namespace std;
  
  int agc = 9;
  string agv[9];

  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(1);
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionEnumeration(enumeration);
  agv[5] = getOptionBatchNumber(batches);
  agv[6] = getOptionBatchLength(iterations);
  agv[7] = getOptionRandomSeed(getRandomSeed());
  agv[8] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  return(inputFile.c_str()); 
}

// [[Rcpp::export]]
std::string RHWtableHE(std::string inputFile, bool enumeration, int dememorization, int batches, int iterations) {
  using namespace std;
  
  int agc = 9;
  string agv[9];

  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(2);
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionEnumeration(enumeration);
  agv[5] = getOptionBatchNumber(batches);
  agv[6] = getOptionBatchLength(iterations);
  agv[7] = getOptionRandomSeed(getRandomSeed());
  agv[8] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  return(inputFile.c_str()); 
}

// [[Rcpp::export]]
std::string RHWtableProbability(std::string inputFile, bool enumeration, int dememorization, int batches, int iterations) {
  using namespace std;
  
  int agc = 9;
  string agv[9];

  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(3);
  agv[3] = getOptionDememorisation(dememorization);
  agv[4] = getOptionEnumeration(enumeration);
  agv[5] = getOptionBatchNumber(batches);
  agv[6] = getOptionBatchLength(iterations);
  agv[7] = getOptionRandomSeed(getRandomSeed());
  agv[8] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  return(inputFile.c_str()); 
}

// [[Rcpp::export]]
std::string RHWtableAlleleFrequenciesExpectedGenotypesFis(std::string inputFile) {
  using namespace std;
  int agc = 5;
  string agv[5];

  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(4);
  agv[3] = getOptionRandomSeed(getRandomSeed());
  agv[4] = getOptionModeBatch();

  printCmd(agc, agv);

  mainJimmy(agc, agv);

  return(inputFile.c_str()); 
}

// [[Rcpp::export]]
std::string RHWtableHDWithSettingsFile(std::string inputFile, std::string settingsFile) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(1);
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  return(inputFile.c_str()); 
}

// [[Rcpp::export]]
std::string RHWtableHEWithSettingsFile(std::string inputFile, std::string settingsFile) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(2);
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  return(inputFile.c_str()); 
}


// [[Rcpp::export]]
std::string RHWtableProbabilityWithSettingsFile(std::string inputFile, std::string settingsFile) {
  using namespace std;
  
  int numberLine = getNumberLineFile(settingsFile);
  const int agc = numberLine+4;
  string *agv=new string[agc];
  
  std::ifstream infile(settingsFile.c_str());
  std::string line;

  int i = 3;
  agv[0] = getNameProg();
  agv[1] = getOptionHWFile(inputFile);
  agv[2] = getOptionHWFileMenu(3);
  while (std::getline(infile, line)) {
      agv[i] = line;
      i++;
  }
  agv[agc-1] = getOptionModeBatch();

  infile.close();

  printCmd(agc, agv);

  mainJimmy(agc, agv);
  delete[] agv;
  
  return(inputFile.c_str()); 
}
