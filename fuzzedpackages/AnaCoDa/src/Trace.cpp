#include "include/base/Trace.h"
#include "include/SequenceSummary.h"
#include <assert.h>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
//----------- Constructors & Destructors -----------//
//--------------------------------------------------//


Trace::Trace()
{
	categories = 0;
	numCodonSpecificParamTypes = 2;
	codonSpecificParameterTrace.resize(numCodonSpecificParamTypes);
	// TODO: fill this
}


Trace::~Trace()
{
	//dtor
}


Trace::Trace(unsigned _numCodonSpecificParamTypes)
{
	categories = 0;
	numCodonSpecificParamTypes = _numCodonSpecificParamTypes;
	codonSpecificParameterTrace.resize(numCodonSpecificParamTypes);
}


//-----------------------------------------------------//
//----------Private Initialization Functions ----------//
//-----------------------------------------------------//


void Trace::initializeSharedTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories,
	unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, std::vector<double> init_phi,
	std::vector<unsigned> init_mix_assign, unsigned numObservedPhiSets, bool estimateSynthesisRate)
{
	my_print("maxGrouping: %\n", maxGrouping);

	//numSelectionCategories always == numSynthesisRateCategories, so only one is passed in for convenience

	initStdDevSynthesisRateTrace(numSelectionCategories, samples);
	initSynthesisRateAcceptanceRateTrace(num_genes, numSelectionCategories);
	codonSpecificAcceptanceRateTrace.resize(maxGrouping);
	initSynthesisRateTrace(samples, num_genes, numSelectionCategories,init_phi,estimateSynthesisRate);
	
	initMixtureAssignmentTrace(samples, num_genes,init_mix_assign);
	initMixtureProbabilitiesTrace(samples, numMixtures);
	initSynthesisOffsetTrace(samples, numObservedPhiSets);
	initObservedSynthesisNoiseTrace(samples, numObservedPhiSets);
	categories = &_categories;
}


void Trace::initStdDevSynthesisRateTrace(unsigned numSelectionCategories, unsigned samples)
{
	stdDevSynthesisRateTrace.resize(numSelectionCategories);
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> temp(samples, 0.0);
		stdDevSynthesisRateTrace[i] = temp;
	}
}


void Trace::initSynthesisRateAcceptanceRateTrace(unsigned num_genes, unsigned numSynthesisRateCategories)
{
	synthesisRateAcceptanceRateTrace.resize(numSynthesisRateCategories);
	for (unsigned category = 0; category < numSynthesisRateCategories; category++)
	{
		synthesisRateAcceptanceRateTrace[category].resize(num_genes);
	}
	//NOTE: this is not sized to samples because push_back takes care of the initialization
}


void Trace::initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numSynthesisRateCategories,std::vector<double> init_phi, bool estimateSynthesisRate)
{
	synthesisRateTrace.resize(numSynthesisRateCategories);
	for (unsigned category = 0; category < numSynthesisRateCategories; category++)
	{
		synthesisRateTrace[category].resize(num_genes);
		for (unsigned i = 0; i < num_genes; i++)
		{
			// if (estimateSynthesisRate)
			// {
				std::vector<float> tempExpr(samples, init_phi[i]);
				synthesisRateTrace[category][i] = tempExpr;
			// }
			// else
			// {
			// 	synthesisRateTrace[category][i].resize(1);
			// 	synthesisRateTrace[category][i][0] = init_phi[i];
			// }
		}
	}
}


void Trace::initMixtureAssignmentTrace(unsigned samples, unsigned num_genes, std::vector<unsigned> init_mix_assign)
{
	mixtureAssignmentTrace.resize(num_genes);
	assert(init_mix_assign.size() == num_genes);
	for (unsigned i = 0u; i < num_genes; i++)
	{
		mixtureAssignmentTrace[i].resize(samples,init_mix_assign[i]);
	}
}


void Trace::initMixtureProbabilitiesTrace(unsigned samples, unsigned numMixtures)
{
	mixtureProbabilitiesTrace.resize(numMixtures);
	for (unsigned i = 0u; i < numMixtures; i++)
	{
		mixtureProbabilitiesTrace[i].resize(samples, 0.0);

	}
}


void Trace::initCodonSpecificParameterTrace(unsigned samples, unsigned numCategories, unsigned numParam, unsigned paramType)
{
	std::vector <std::vector <std::vector <float>>> tmp;
	tmp.resize(numCategories);
	for (unsigned category = 0; category < numCategories; category++)
	{
		tmp[category].resize(numParam);
		for (unsigned i = 0; i < numParam; i++)
		{
			std::vector <float> temp(samples, 0.0);
			tmp[category][i] = temp;
		}
	}


	//TODO: R output for error message here
	codonSpecificParameterTrace[paramType] = tmp;
	/*
	switch (paramType) {
	case 0:
		codonSpecificParameterTraceOne = tmp;
		break;
	case 1:
		codonSpecificParameterTraceTwo = tmp;
		break;
	default:
		my_printError("ERROR: Invalid paramType given, codon specific parameter trace not initialized.\n");
		break;
	}
	*/
}



void Trace::initSynthesisOffsetTrace(unsigned samples, unsigned numPhiGroupings)
{
	synthesisOffsetTrace.resize(numPhiGroupings);
	for (unsigned i = 0; i < numPhiGroupings; i++)
	{
		synthesisOffsetTrace[i].resize(samples);
	}

	synthesisOffsetAcceptanceRateTrace.resize(numPhiGroupings);
}


void Trace::initObservedSynthesisNoiseTrace(unsigned samples, unsigned numPhiGroupings)
{
	observedSynthesisNoiseTrace.resize(numPhiGroupings);
	for (unsigned i = 0; i < numPhiGroupings; i++)
	{
		observedSynthesisNoiseTrace[i].resize(samples);
	}
}


//----------------------------------//
//--------- PANSE Specific ---------//
//----------------------------------//
void Trace::initPartitionFunctionTrace(unsigned samples, unsigned numPartitionFunctionsGroupings)
{
    partitionFunctionTrace.resize(numPartitionFunctionsGroupings);
    for (unsigned i = 0; i < numPartitionFunctionsGroupings; i++)
    {
        partitionFunctionTrace[i].resize(samples);
    }
}



//------------------------------------//
//----------- FONSE Specific ---------//
//------------------------------------//

void Trace::initInitiationCostTrace(unsigned samples)
{
	initiationCostTrace.resize(samples);	
}


//----------------------------------------------------//
//---------- Model Initialization Functions ----------//
//----------------------------------------------------//

/* Note 1) -- on these Model Initialization Functions' initCodonSpecificParameterTrace functions
 * The last argument specifies the codon specific parameter type (dM and dEta for ROC)
 * You can check Parameter.cpp to check what the values of dM and dEta are.
 * The constants aren't used here because they are not available to the Trace object.
*/

void Trace::initializePATrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
	unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numObservedPhiSets, std::vector<double> init_phi,
	std::vector<unsigned> init_mix_assign,bool estimateSynthesisRate)
{
	initializeSharedTraces(samples, num_genes, numLambdaPrimeCategories, numMixtures,
		_categories, maxGrouping, init_phi, init_mix_assign,numObservedPhiSets,estimateSynthesisRate);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numAlphaCategories,  numParam, 0u); // alp
	initCodonSpecificParameterTrace(samples, numLambdaPrimeCategories, numParam, 1u); // lmPri

}


void Trace::initializeROCTrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
	unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numObservedPhiSets,std::vector<double> init_phi,
	std::vector<unsigned> init_mix_assign,bool estimateSynthesisRate)
{
	initializeSharedTraces(samples, num_genes, numSelectionCategories, numMixtures, _categories, maxGrouping,init_phi,init_mix_assign,estimateSynthesisRate);


	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numMutationCategories, numParam, 0u); // dM
	initCodonSpecificParameterTrace(samples, numSelectionCategories, numParam, 1u); // dEta

	// initSynthesisOffsetTrace(samples, numObservedPhiSets);
	// initObservedSynthesisNoiseTrace(samples, numObservedPhiSets);
}


void Trace::initializeFONSETrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
	unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numObservedPhiSets,std::vector<double> init_phi,
	std::vector<unsigned> init_mix_assign,bool estimateSynthesisRate)
{
	initializeSharedTraces(samples, num_genes, numSelectionCategories, numMixtures,
		 _categories, maxGrouping,init_phi,init_mix_assign,numObservedPhiSets,estimateSynthesisRate);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numMutationCategories, numParam, 0u); // dM
	initCodonSpecificParameterTrace(samples, numSelectionCategories, numParam, 1u); // dOmega
	initInitiationCostTrace(samples);
}


void Trace::initializePANSETrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
	unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numObservedPhiSets,std::vector<double> init_phi,
	std::vector<unsigned> init_mix_assign, bool estimateSynthesisRate)
{
    numCodonSpecificParamTypes = 3;
    codonSpecificParameterTrace.resize(numCodonSpecificParamTypes);

	initializeSharedTraces(samples, num_genes, numLambdaPrimeCategories, numMixtures,
		_categories, maxGrouping,init_phi,init_mix_assign,numObservedPhiSets,estimateSynthesisRate);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numAlphaCategories,  numParam, 0u); //alp
	initCodonSpecificParameterTrace(samples, numLambdaPrimeCategories, numParam, 1u); //lmPri
    initCodonSpecificParameterTrace(samples, numAlphaCategories, numParam, 2u); //nseRate
    initPartitionFunctionTrace(samples, numMixtures);

    nseSpecificAcceptanceRateTrace.resize(maxGrouping);
}


//--------------------------------------//
//---------- Getter Functions ----------//
//--------------------------------------//


std::vector<double> Trace::getStdDevSynthesisRateTrace(unsigned selectionCategory)
{
	return stdDevSynthesisRateTrace[selectionCategory];
}


std::vector<double> Trace::getExpectedSynthesisRateTrace()
{
	unsigned numGenes = (unsigned)synthesisRateTrace[0].size(); //number of genes
	unsigned samples = (unsigned)synthesisRateTrace[0][0].size(); //number of samples
	std::vector<double> RV(samples, 0.0);
	for (unsigned sample = 0; sample < samples; sample++)
	{
		for (unsigned geneIndex = 0; geneIndex < numGenes; geneIndex++)
		{
			unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][sample];
			unsigned category = getSynthesisRateCategory(mixtureElement);
			RV[sample] += synthesisRateTrace[category][geneIndex][sample];
		}
		RV[sample] /= numGenes;
	}
	return RV;
}


std::vector<double> Trace::getStdDevSynthesisRateAcceptanceRateTrace()
{
	return stdDevSynthesisRateAcceptanceRateTrace;
}


std::vector<std::vector<std::vector<float>>> Trace::getSynthesisRateTrace()
{
	return synthesisRateTrace;
}


std::vector<double> Trace::getSynthesisRateAcceptanceRateTraceByMixtureElementForGene(unsigned mixtureElement,
	unsigned geneIndex)
{
	unsigned category = getSynthesisRateCategory(mixtureElement);
	return synthesisRateAcceptanceRateTrace[category][geneIndex];
}


std::vector<std::vector<std::vector<double>>> Trace::getSynthesisRateAcceptanceRateTrace()
{
	return synthesisRateAcceptanceRateTrace;
}


std::vector<double> Trace::getCodonSpecificAcceptanceRateTraceForAA(std::string aa)
{
	aa[0] = (char)std::toupper(aa[0]);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
	return codonSpecificAcceptanceRateTrace[aaIndex];
}

std::vector<double> Trace::getCodonSpecificAcceptanceRateTraceForCodon(std::string codon)
{
	codon[0] = (char)std::toupper(codon[0]);
	unsigned codonIndex = SequenceSummary::codonToIndexWithReference.find(codon)->second;
	return codonSpecificAcceptanceRateTrace[codonIndex];
}



std::vector<float> Trace::getSynthesisRateTraceForGene(unsigned geneIndex)
{
	unsigned traceLength = (unsigned)synthesisRateTrace[0][0].size();

	std::vector<float> returnVector(traceLength, 0.0);
	for (unsigned i = 0u; i < traceLength; i++)
	{
		unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][i];
		unsigned category = getSynthesisRateCategory(mixtureElement);
		returnVector[i] = synthesisRateTrace[category][geneIndex][i];
	}
	return returnVector;
}


std::vector<float> Trace::getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
{


	unsigned category = getSynthesisRateCategory(mixtureElement);
	return synthesisRateTrace[category][geneIndex];
}


std::vector<unsigned> Trace::getMixtureAssignmentTraceForGene(unsigned geneIndex)
{
	return mixtureAssignmentTrace[geneIndex];
}


std::vector<double> Trace::getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex)
{
	return mixtureProbabilitiesTrace[mixtureIndex];
}


std::vector<std::vector<unsigned>> Trace::getMixtureAssignmentTrace()
{
	return mixtureAssignmentTrace;
}


std::vector<std::vector<double>> Trace::getMixtureProbabilitiesTrace()
{
	return mixtureProbabilitiesTrace;
}


std::vector<std::vector<double>> Trace::getCodonSpecificAcceptanceRateTrace()
{
	return codonSpecificAcceptanceRateTrace;
}

std::vector<std::vector<double>> Trace::getNseRateSpecificAcceptanceRateTrace()
{
	return nseSpecificAcceptanceRateTrace;
}


unsigned Trace::getSynthesisRateCategory(unsigned mixtureElement)
{
	return categories->at(mixtureElement).delEta;
}


std::vector<std::vector<std::vector<std::vector<float>>>>* Trace::getCodonSpecificParameterTrace()
{
	return &codonSpecificParameterTrace;
}


//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//

std::vector<float> Trace::getCodonSpecificParameterTraceByGeneElementForCodon(unsigned geneIndex,
	std::string& codon, unsigned paramType, bool withoutReference)
{
	std::vector <float> rv;
	unsigned codonIndex = SequenceSummary::codonToIndex(codon, withoutReference);
	unsigned tracelength = (unsigned)codonSpecificParameterTrace[0][0][0].size();
	rv.resize(tracelength);
	for(unsigned i = 0; i < tracelength; i++)
	{
		unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][i];
		unsigned category = getCodonSpecificCategory(mixtureElement, paramType);
		rv[i] = codonSpecificParameterTrace[paramType][category][codonIndex][i];
	}

	return rv;
}


std::vector<float> Trace::getCodonSpecificParameterTraceByMixtureElementForCodon(unsigned mixtureElement,
	std::string& codon, unsigned paramType, bool withoutReference)
{
	std::vector <float> rv;
	unsigned codonIndex = SequenceSummary::codonToIndex(codon, withoutReference);
	unsigned category = getCodonSpecificCategory(mixtureElement, paramType);
	rv = codonSpecificParameterTrace[paramType][category][codonIndex];
	return rv;
}


std::vector<double> Trace::getSynthesisOffsetTrace(unsigned index)
{
	return synthesisOffsetTrace[index];
}


std::vector<double> Trace::getSynthesisOffsetAcceptanceRateTraceForIndex(unsigned index)
{
	return synthesisOffsetAcceptanceRateTrace[index];
}


std::vector<double> Trace::getObservedSynthesisNoiseTrace(unsigned index)
{
	return observedSynthesisNoiseTrace[index];
}


std::vector<std::vector<std::vector<float>>> Trace::getCodonSpecificParameterTraceByParamType(unsigned paramType)
{
	return codonSpecificParameterTrace[paramType];
}


std::vector<std::vector<double>> Trace::getSynthesisOffsetAcceptanceRateTrace()
{
	return synthesisOffsetAcceptanceRateTrace;
}

//TO DO: How does this work for PA/PANSE
unsigned Trace::getCodonSpecificCategory(unsigned mixtureElement, unsigned paramType)
{
	unsigned rv = 0;
	switch (paramType) {
	case 0:
		rv = categories->at(mixtureElement).delM;
		break;
	case 1:
		rv = categories->at(mixtureElement).delEta;
		break;
    case 2:
        rv = categories->at(mixtureElement).delM;
        break;
	default:
		my_printError("ERROR: Unknown parameter type in getCodonSpecificCategory\n");
		break;
	}
	return rv;
}


//--------------------------------------//
//----------- PANSE Specific -----------//
//--------------------------------------//
std::vector<double> Trace::getPartitionFunctionTrace(unsigned mixtureIndex)
{
    return partitionFunctionTrace[mixtureIndex];
}

std::vector<std::vector<double>> Trace::getPartitionFunctionTraces()
{
    return partitionFunctionTrace;
}


std::vector<double> Trace::getPartitionFunctionAcceptanceRateTrace()
{
    return partitionFunctionTraceAcceptanceRateTrace;
}


//--------------------------------------//
//----------- FONSE Specific -----------//
//--------------------------------------//
std::vector<double> Trace::getInitiationCostTrace()
{
    return initiationCostTrace;
}


std::vector<double> Trace::getInitiationCostAcceptanceRateTrace()
{
    return initiationCostAcceptanceRateTrace;
}



//--------------------------------------//
//---------- Update Functions ----------//
//--------------------------------------//


void Trace::updateStdDevSynthesisRateTrace(unsigned sample, double stdDevSynthesisRate, unsigned synthesisRateCategory)
{
	stdDevSynthesisRateTrace[synthesisRateCategory][sample] = stdDevSynthesisRate;
}


void Trace::updateStdDevSynthesisRateAcceptanceRateTrace(double acceptanceLevel)
{
	stdDevSynthesisRateAcceptanceRateTrace.push_back(acceptanceLevel);
}


void Trace::updateSynthesisRateAcceptanceRateTrace(unsigned category, unsigned geneIndex, double acceptanceLevel)
{
	synthesisRateAcceptanceRateTrace[category][geneIndex].push_back(acceptanceLevel);
}


void Trace::updateCodonSpecificAcceptanceRateTrace(unsigned codonIndex, double acceptanceLevel)
{
	codonSpecificAcceptanceRateTrace[codonIndex].push_back(acceptanceLevel);
}



void Trace::updateSynthesisRateTrace(unsigned sample, unsigned geneIndex,
	std::vector<std::vector <double>> &currentSynthesisRateLevel)
{
	for (unsigned category = 0; category < synthesisRateTrace.size(); category++)
	{
		synthesisRateTrace[category][geneIndex][sample] = currentSynthesisRateLevel[category][geneIndex];
	}
}


void Trace::updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value)
{
	mixtureAssignmentTrace[geneIndex][sample] = value;
}


void Trace::updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities)
{
	for (unsigned category = 0; category < mixtureProbabilitiesTrace.size(); category++)
	{
		mixtureProbabilitiesTrace[category][samples] = categoryProbabilities[category];
	}
}


//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//


void Trace::updateCodonSpecificParameterTraceForAA(unsigned sample, std::string aa,
	std::vector<std::vector<double>> &curParam, unsigned paramType)
{
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
	for (unsigned category = 0; category < codonSpecificParameterTrace[paramType].size(); category++)
	{
		for (unsigned i = aaStart; i < aaEnd; i++)
		{
			//my_print("% % % %\n",category,i,aa,curParam[category][i]);
			codonSpecificParameterTrace[paramType][category][i][sample] = curParam[category][i];
		}
	}

}


void Trace::updateSynthesisOffsetTrace(unsigned index, unsigned sample, double value)
{
	synthesisOffsetTrace[index][sample] = value;
}


void Trace::updateSynthesisOffsetAcceptanceRateTrace(unsigned index, double value)
{
	synthesisOffsetAcceptanceRateTrace[index].push_back(value);
}


void Trace::updateObservedSynthesisNoiseTrace(unsigned index, unsigned sample, double value)
{
	observedSynthesisNoiseTrace[index][sample] = value;
}


//-------------------------------------//
//-------- PA/PANSE Specific ----------//
//-------------------------------------//

void Trace::updateCodonSpecificParameterTraceForCodon(unsigned sample, std::string codon,
	std::vector<std::vector<double>> &curParam, unsigned paramType)
{
	unsigned i = SequenceSummary::codonToIndex(codon);
	for (unsigned category = 0; category < codonSpecificParameterTrace[paramType].size(); category++)
	{
//        if(std::isnan(curParam[category][i])){
//            my_printError("\n Trace::updateCodonSpecificParameterTraceForCodon: Current parameter set contains NaN. \n");
//        }
		codonSpecificParameterTrace[paramType][category][i][sample] = curParam[category][i];
	}
}

void Trace::updatePartitionFunctionTrace(unsigned index, unsigned sample, double value)
{
   partitionFunctionTrace[index][sample] = value;
}


void Trace:: updatePartitionFunctionAcceptanceRateTrace(double value)
{
    partitionFunctionTraceAcceptanceRateTrace.push_back(value);
}

void Trace::updateNseRateSpecificAcceptanceRateTrace(unsigned codonIndex, double acceptanceLevel)
{
	nseSpecificAcceptanceRateTrace[codonIndex].push_back(acceptanceLevel);
}

//------------------------------------//
//----------- FONSE Specific ---------//
//------------------------------------//

void Trace::updateInitiationCostTrace(unsigned sample, double value)
{
   initiationCostTrace[sample] = value;
}


void Trace:: updateInitiationCostAcceptanceRateTrace(double value)
{
    initiationCostAcceptanceRateTrace.push_back(value);
}








// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//

#ifndef STANDALONE

//--------------------------------------//
//---------- Getter Functions ----------//
//--------------------------------------//
std::vector<double> Trace::getSynthesisRateAcceptanceRateTraceByMixtureElementForGeneR(unsigned mixtureElement,
	unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateAcceptanceRateTrace.size());
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	if (checkGene && checkMixtureElement)
	{
		RV = getSynthesisRateAcceptanceRateTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<float> Trace::getSynthesisRateTraceForGeneR(unsigned geneIndex)
{
	std::vector<float> RV;
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateTrace[0].size());
	if (checkGene)
	{
		RV = getSynthesisRateTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<float> Trace::getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex)
{
	std::vector<float> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateTrace[0].size());
	if (checkMixtureElement && checkGene)
	{
		RV = getSynthesisRateTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<unsigned> Trace::getMixtureAssignmentTraceForGeneR(unsigned geneIndex)
{
	std::vector <unsigned> RV;
	bool checkGene = checkIndex(geneIndex, 1, mixtureAssignmentTrace.size());
	if (checkGene)
	{
		RV = getMixtureAssignmentTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex)
{
	std::vector<double> RV;
	bool check = checkIndex(mixtureIndex, 1, mixtureProbabilitiesTrace.size());
	if (check)
	{
		RV = getMixtureProbabilitiesTraceForMixture(mixtureIndex - 1);
	}
	return RV;
}


std::vector<std::vector<double>> Trace::getStdDevSynthesisRateTraces()
{
	return stdDevSynthesisRateTrace;
}


unsigned Trace::getNumberOfMixtures()
{
	return mixtureProbabilitiesTrace.size();
}

std::vector<double> Trace::getPartitionFunctionTraceR(unsigned mixtureIndex){
    std::vector<double> RV;
    bool check = checkIndex(mixtureIndex, 1, partitionFunctionTrace.size());
    if (check)
    {
        RV = getPartitionFunctionTrace(mixtureIndex - 1);
    }
    return RV;
}

//--------------------------------------//
//---------- Setter Functions ----------//
//--------------------------------------//
void Trace::setStdDevSynthesisRateTraces(std::vector<std::vector<double>> _stdDevSynthesisRateTrace)
{
	stdDevSynthesisRateTrace = _stdDevSynthesisRateTrace;
}


void Trace::setStdDevSynthesisRateAcceptanceRateTrace(std::vector<double> _stdDevSynthesisRateAcceptanceRateTrace)
{
	stdDevSynthesisRateAcceptanceRateTrace = _stdDevSynthesisRateAcceptanceRateTrace;
}


void Trace::setSynthesisRateTrace(std::vector<std::vector<std::vector<float>>> _synthesisRateTrace)
{
	synthesisRateTrace = _synthesisRateTrace;
}


void Trace::setSynthesisRateAcceptanceRateTrace(std::vector<std::vector<std::vector<double>>>_synthesisRateAcceptanceRateTrace)
{
	synthesisRateAcceptanceRateTrace = _synthesisRateAcceptanceRateTrace;
}


void Trace::setMixtureAssignmentTrace(std::vector<std::vector<unsigned>> _mixtureAssignmentTrace)
{
	mixtureAssignmentTrace = _mixtureAssignmentTrace;
}


void Trace::setMixtureProbabilitiesTrace(std::vector<std::vector<double>> _mixtureProbabilitiesTrace)
{
	mixtureProbabilitiesTrace = _mixtureProbabilitiesTrace;
}


void Trace::setCodonSpecificAcceptanceRateTrace(std::vector<std::vector<double>> _cspAcceptanceRateTrace)
{
	codonSpecificAcceptanceRateTrace = _cspAcceptanceRateTrace;
}

void Trace::setNseRateSpecificAcceptanceRateTrace(std::vector<std::vector<double>> _nseAcceptanceRateTrace)
{
	nseSpecificAcceptanceRateTrace = _nseAcceptanceRateTrace;
}


void Trace::setCategories(std::vector<mixtureDefinition> &_categories)
{
	categories = &_categories;
}

//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//
std::vector<float> Trace::getCodonSpecificParameterTraceByMixtureElementForCodonR(unsigned mixtureElement,
	std::string& codon, unsigned paramType, bool withoutReference)
{
	std::vector<float> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, getNumberOfMixtures());
	if (checkMixtureElement)
	{
		RV = getCodonSpecificParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon, paramType, withoutReference);
	}
	return RV;
}

std::vector<std::vector<double>> Trace::getSynthesisOffsetTraceR()
{
	return synthesisOffsetTrace;
}


std::vector<std::vector<double>> Trace::getObservedSynthesisNoiseTraceR()
{
	return observedSynthesisNoiseTrace;
}


void Trace::setSynthesisOffsetTrace(std::vector<std::vector <double> > _NoiseOffsetTrace)
{
	synthesisOffsetTrace = _NoiseOffsetTrace;
}


void Trace::setSynthesisOffsetAcceptanceRateTrace(std::vector<std::vector <double> > _NoiseOffsetAcceptanceRateTrace)
{
	synthesisOffsetAcceptanceRateTrace = _NoiseOffsetAcceptanceRateTrace;
}


void Trace::setObservedSynthesisNoiseTrace(std::vector<std::vector <double> > _ObservedSynthesisNoiseTrace)
{
	observedSynthesisNoiseTrace = _ObservedSynthesisNoiseTrace;
}


void Trace::setCodonSpecificParameterTrace(std::vector<std::vector<std::vector<float>>> _parameterTrace, unsigned paramType)
{
	codonSpecificParameterTrace[paramType] = _parameterTrace;
	/*
	switch (paramType) {
	case 0:
		codonSpecificParameterTraceOne = _parameterTrace;
		break;
	case 1:
		codonSpecificParameterTraceTwo = _parameterTrace;
		break;
	default:
		my_printError("ERROR: Unknown parameter type in setCodonSpecificParameterTrace.\n");
		break;
	}
	*/
}


//----------------------------------//
//--------- PANSE Specific ---------//
//----------------------------------//
void Trace::setPartitionFunctionTraces(std::vector<std::vector <double> > _PartitionFunctionTrace)
{
    partitionFunctionTrace = _PartitionFunctionTrace;
}


//----------------------------------//
//--------- FONSE Specific ---------//
//----------------------------------//
void Trace::setInitiationCostTrace(std::vector <double> _InitiationCostTrace)
{
    initiationCostTrace = _InitiationCostTrace;
}





bool Trace::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		my_printError("ERROR: Index: % is out of bounds. Index must be between % & %\n", index, lowerbound, upperbound);
	}

	return check;
}

#endif


//Function implemented to resize codonSpecificParameterTrace for reloading parameter object.
//TO DO: Think of more elegant solution. 
void Trace::resizeNumberCodonSpecificParameterTrace(unsigned _numCodonSpecificParamTypes)
{
	codonSpecificParameterTrace.resize(_numCodonSpecificParamTypes);
}
