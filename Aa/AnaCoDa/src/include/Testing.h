#ifndef Testing_H
#define Testing_H


#include "Utility.h"
#include "SequenceSummary.h"
#include "Gene.h"
#include "Genome.h"
#include "base/Parameter.h"
#include "CovarianceMatrix.h"
#include "MCMCAlgorithm.h"


int testUtility();
int testSequenceSummary();
int testGene();
void testGenomePAHelper(Genome* genome, bool simulated); //Only used internally for testGenome
bool testGenomeSimulatedPAEqualityHelper(Genome genome1, Genome genome2); //Only used internally for testGenome
int testGenome(std::string testFileDir);
int testParameter(std::string testFileDir);
//int testParameterWithFile(std::string filename); //TODO: Rework or remove
int testCovarianceMatrix();
//int testPAParameter(); //TODO: Rework or remove
int testMCMCAlgorithm();

//Blank header
#endif // Testing_H
