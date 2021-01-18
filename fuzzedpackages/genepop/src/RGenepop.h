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

#ifndef RGENEPOP_H
#define RGENEPOP_H

#include <iostream>
#include <cstring>

#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
void set_seed(unsigned int seed);
#endif


/* Tools */
std::string getNameProg();
void printCmd(int size, std::string argv[]);
long getRandomSeed();
long getMantelSeed();
void setRandomSeed(long seed);
void setMantelSeed(long seed);
void resetSeed();

void setVerboseGenepop(bool verbose);

void initRGenepop();
void cleanRGenepop();

int getNumberLineFile(std::string file);

std::string getVersion();

/* Option Genepop */
std::string getOptionInputFile(std::string inputFile);
std::string getOptionMenu(std::string menuOptions);
std::string getOptionEnumeration(bool enumeration);
std::string getOptionDememorisation(int dememorization);
std::string getOptionBatchNumber(int batches);
std::string getOptionBatchLength(int iterations);
std::string getOptionNullAlleleMethod(std::string method);
std::string getOptionICoverage(double icoverage);
std::string getOptionEstimationPloidy(std::string dataType);
std::string getOptionStructFile(std::string inputFile);
std::string getOptionIsolationStatistic(std::string statistic);
std::string getOptionGeographicScale(std::string scale);
std::string getOptionTestPoint(double testPoint);
std::string getOptionMinimalDistance(double minDist);
std::string getOptionMaximalDistance(double maxDist);
std::string getOptionMantelPermutations(int permutations);
std::string getOptionMantelRankTest(bool mantelRank);
std::string getOptionIsolationFile(std::string inputFile);
std::string getOptionRandomSeed(long seed);
std::string getOptionMantelSeed(long seed);
std::string getOptionAllelicDistance(std::string distance);
std::string getOptionAlleleSizes(std::string sizes);
std::string getOptionHWFile(std::string inputFile);
std::string getOptionHWFileMenu(int option);
std::string getOptionPopTypes(std::string popTypes);
std::string getOptionPopTypeSelection(std::string popTypesSelection);
std::string getOptionMultiMigFile(std::string inputFile);
std::string getOptionGeoDistFile(std::string inputFile);
std::string getOptionSettingsFile(std::string settingsFile);
std::string getOptionModeBatch();

/* Extension output file */

std::string getOutPutFileMenu_1_1(std::string inputFile);
std::string getOutPutFileMenu_1_2(std::string inputFile);
std::string getOutPutFileMenu_1_3(std::string inputFile);
std::string getOutPutFileMenu_1_4(std::string inputFile);
std::string getOutPutFileMenu_1_5(std::string inputFile);
std::string getOutPutFileMenu_2_1(std::string inputFile);
std::string getOutPutFileMenu_2_2(std::string inputFile);
std::string getOutPutFileMenu_3_1(std::string inputFile);
std::string getOutPutFileMenu_3_2(std::string inputFile);
std::string getOutPutFileMenu_3_3(std::string inputFile);
std::string getOutPutFileMenu_3_4(std::string inputFile);
std::string getOutPutFileMenu_4_1(std::string inputFile);
std::string getOutPutFileMenu_5_1(std::string inputFile);
std::string getOutPutFileMenu_5_2(std::string inputFile);
std::string getOutPutFileMenu_5_3(std::string inputFile);
std::string getOutPutFileMenu_6_1(std::string inputFile);
std::string getOutPutFileMenu_6_2(std::string inputFile);
std::string getOutPutFileMenu_6_2_b(std::string inputFile);
std::string getOutPutFileMenu_6_3(std::string inputFile);
std::string getOutPutFileMenu_6_4(std::string inputFile);
std::string getOutPutFileMenu_6_4_b(std::string inputFile);
std::string getOutPutFileMenu_6_5(std::string inputFile);
std::string getOutPutFileMenu_6_5_b(std::string inputFile);
std::string getOutPutFileMenu_6_5_c(std::string inputFile);
std::string getOutPutFileMenu_6_6(std::string inputFile);
std::string getOutPutFileMenu_6_6_b(std::string inputFile);
std::string getOutPutFileMenu_6_6_c(std::string inputFile);
std::string getOutPutFileMenu_7_1(std::string inputFile);
std::string getOutPutFileMenu_7_2(std::string inputFile);
std::string getOutPutFileMenu_7_3(std::string inputFile);
std::string getOutPutFileMenu_7_4(std::string inputFile);
std::string getOutPutFileMenu_8_1(std::string inputFile);
std::string getOutPutFileMenu_8_2(std::string inputFile);
std::string getOutPutFileMenu_8_3(std::string inputFile);
std::string getOutPutFileMenu_8_4(std::string inputFile);
std::string getOutPutFileMenu_8_5(std::string inputFile);
std::string getOutPutFileMenu_8_6(std::string inputFile);

/* Functionnality */

/* 1 : Exact Hardy-Weinberg tests */
std::string RHWEachLocusEachPopulationHD(std::string inputFile, std::string outputFile, bool enumeration, int dememorizationn, int batches, int iterations );
std::string RHWEachLocusEachPopulationHE(std::string inputFile, std::string outputFile, bool enumeration, int dememorization, int batches, int iterations );
std::string RHWEachLocusEachPopulationProbability(std::string inputFile, std::string outputFile, bool enumeration, int dememorization, int batches, int iterations );
std::string RHWGlobalHD(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );
std::string RHWGlobalHE(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );

/* 1 : with setting file */
std::string RHWEachLocusEachPopulationHDWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RHWEachLocusEachPopulationHEWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RHWEachLocusEachPopulationProbabilityWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RHWGlobalHDWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RHWGlobalHEWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );

/* 2 : Exact Genotypic Disequilibrium tests */
std::string RGDEachPairLociEachPopulation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );
std::string RGDGenotypicContingency(std::string inputFile, std::string outputFile );

/* 2 : with setting file */
std::string RGDEachPairLociEachPopulationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );

/* 3 : Exact Population Differentiation tests */
std::string RPDGenicAllPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );
std::string RPDGenicAllPairPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );
std::string RPDGenotypicAllPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );
std::string RPDGenotypicAllPairPopulationDifferentiation(std::string inputFile, std::string outputFile, int dememorization, int batches, int iterations );

/* 3 : with setting file */
std::string RPDGenicAllPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RPDGenicAllPairPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RPDGenotypicAllPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RPDGenotypicAllPairPopulationDifferentiationWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );

//output in input
std::string RAnalyzingSingleContingencyTable(std::string inputFile, int dememorization, int batches, int iterations);
std::string RAnalyzingSingleContingencyTableWithSettingsFile(std::string inputFile, std::string settingsFile);

/* 4 :  Nm estimates */
std::string RNmEstimates(std::string inputFile, std::string outputFile, std::string dataType );

/* 5 : Descriptif */
std::string RDescriptifAlleleAndGenotypeFrequenciesPerLocusPerSample(std::string inputFile, std::string outputFile );
std::string RDescriptifGeneDiversitiesAndFisUsingAlleleIdentity(std::string inputFile, std::string outputFile, std::string dataType );
std::string RDescriptifGeneDiversitiesAndFisUsingAlleleSize(std::string inputFile, std::string outputFile, std::string dataType );

/* 6 : FstIBD */
std::string REstimatingSpatialStructureAlleleIdentyAllPopulations(std::string inputFile, std::string outputFile, std::string dataType );
std::string REstimatingSpatialStructureAlleleIdentyAllPopulationsPairs(std::string inputFile, std::string outputFile, std::string dataType );
std::string REstimatingSpatialStructureAlleleSizeAllPopulations(std::string inputFile, std::string outputFile, std::string dataType );
std::string REstimatingSpatialStructureAlleleSizeAllPopulationsPairs(std::string inputFile, std::string outputFile, std::string dataType );

std::string RIsolationByDistanceBetweenIndividuals (std::string inputFile, std::string outputFile, std::string dataType, std::string statistic, std::string geographicScale, double CIcoverage, double testPoint, double minimalDistance, double maximalDistance, int mantelPermutations, bool mantelRankTest );
std::string RIsolationByDistanceBetweenGroups (std::string inputFile, std::string outputFile, std::string dataType, std::string statistic, std::string geographicScale, double CIcoverage, double testPoint, double minimalDistance, double maximalDistance, int mantelPermutations, bool mantelRankTest );
/* 6 : with setting file */

std::string RIsolationByDistanceBetweenIndividualsWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );
std::string RIsolationByDistanceBetweenGroupsWithSettingsFile(std::string inputFile, std::string outputFile, std::string settingsFile );

/* 7 : Ecumenicism */
std::string REcumenicismFstat(std::string inputFile, std::string outputFile );
std::string REcumenicismBiosysLetter(std::string inputFile, std::string outputFile );
std::string REcumenicismBiosysNumber(std::string inputFile, std::string outputFile );
std::string REcumenicismLinkdos(std::string inputFile, std::string outputFile );

/* 8 : Misc */
std::string RNullAlleleEstimateAlleleFrequencies(std::string inputFile, std::string outputFile, std::string nullAlleleMethod, double CIcoverage );
std::string RDiploidisationHaploidData(std::string inputFile, std::string outputFile );
std::string RRelabelingAlleles(std::string inputFile, std::string outputFile );
std::string RConversionToIndividualDataWithPopulationNames(std::string inputFile, std::string outputFile );
std::string RConversionToIndividualDataWithIndividualNames(std::string inputFile, std::string outputFile );
std::string RRandomSamplingOfHaploidGenotypesFromDiploidOnes(std::string inputFile, std::string outputFile );

/* 8 : with setting file */
std::string RNullAlleleEstimateAlleleFrequenciesWithSettingsFile(std::string inputFile, 
                                                                 std::string outputFile, std::string settingsFile);

/* HWfile */

// output in input
std::string RHWtableHD(std::string inputFile, bool enumeration, int dememorization, int batches, int iterations);
std::string RHWtableHE(std::string inputFile, bool enumeration, int dememorization, int batches, int iterations);
std::string RHWtableProbability(std::string inputFile, bool enumeration, int dememorization, int batches, int iterations);
std::string RHWtableAlleleFrequenciesExpectedGenotypesFis(std::string inputFile);

/* HWfile : with setting file */
std::string RHWtableHDWithSettingsFile(std::string inputFile, std::string settingsFile);
std::string RHWtableHEWithSettingsFile(std::string inputFile, std::string settingsFile);
std::string RHWtableProbabilityWithSettingsFile(std::string inputFile, std::string settingsFile);



#endif
