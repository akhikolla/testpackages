#include "include/base/Parameter.h"
#include <cfloat>

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


//C++ runs only
#ifdef STANDALONE
std::default_random_engine Parameter::generator( (unsigned) std::time(NULL));
#endif


//Definition of constant variables
const std::string Parameter::allUnique = "allUnique";
const std::string Parameter::selectionShared = "selectionShared";
const std::string Parameter::mutationShared = "mutationShared";

const unsigned Parameter::dM = 0;
const unsigned Parameter::dEta = 1;
const unsigned Parameter::dOmega = 1;
const unsigned Parameter::alp = 0;
const unsigned Parameter::lmPri = 1;
const unsigned Parameter::nse = 2;



//-------------------------------------------------//
//---------- Constructors & Destructors -----------//
//-------------------------------------------------//


Parameter::Parameter()
{
	lastIteration = 0u;
	numParam = 0u;
	obsPhiSets = 0u;
	adaptiveStepPrev = 0;
	adaptiveStepCurr = 0;
	stdDevSynthesisRate.resize(1);
	stdDevSynthesisRate_proposed.resize(1);
	numAcceptForStdDevSynthesisRate = 0u;
	bias_stdDevSynthesisRate = 0.0;
	bias_phi = 0.0;
	numMutationCategories = 0u;
	numSelectionCategories = 0u;
	numMixtures = 0u;
	std_stdDevSynthesisRate = 0.1;
	maxGrouping = 22;
}


Parameter::Parameter(unsigned _maxGrouping)
{
	lastIteration = 0u;
	numParam = 0u;
	obsPhiSets = 0u;
	adaptiveStepPrev = 0;
	adaptiveStepCurr = 0;
	stdDevSynthesisRate.resize(1);
	stdDevSynthesisRate_proposed.resize(1);
	numAcceptForStdDevSynthesisRate = 0u;
	bias_stdDevSynthesisRate = 0.0;
	bias_phi = 0.0;
	numMutationCategories = 0u;
	numSelectionCategories = 0u;
	numMixtures = 0u;
	std_stdDevSynthesisRate = 0.1;
	maxGrouping = _maxGrouping;
	numAcceptForCodonSpecificParameters.resize(maxGrouping, 0u);
}


Parameter& Parameter::operator=(const Parameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	numParam = rhs.numParam;

	stdDevSynthesisRate = rhs.stdDevSynthesisRate;
	stdDevSynthesisRate_proposed = rhs.stdDevSynthesisRate_proposed;

	numAcceptForStdDevSynthesisRate = rhs.numAcceptForStdDevSynthesisRate;
	obsPhiSets = rhs.obsPhiSets;
	categories = rhs.categories;

	adaptiveStepPrev = rhs.adaptiveStepPrev;
	adaptiveStepCurr = rhs.adaptiveStepCurr;

  	// proposal bias and std for phi values
  	bias_stdDevSynthesisRate = rhs.bias_stdDevSynthesisRate;
  	std_stdDevSynthesisRate = rhs.std_stdDevSynthesisRate;

  	// proposal bias and std for phi values
  	bias_phi = rhs.bias_phi;
  	std_phi = rhs.std_phi;

  	currentSynthesisRateLevel = rhs.currentSynthesisRateLevel;
  	proposedSynthesisRateLevel = rhs.proposedSynthesisRateLevel;
  	numAcceptForSynthesisRate = rhs.numAcceptForSynthesisRate;

  	numMutationCategories = rhs.numMutationCategories;
  	numSelectionCategories = rhs.numSelectionCategories;

	proposedCodonSpecificParameter = rhs.proposedCodonSpecificParameter;
	currentCodonSpecificParameter = rhs.currentCodonSpecificParameter;

  	numMixtures = rhs.numMixtures;

  	mutationSelectionState = rhs.mutationSelectionState;
  	selectionIsInMixture = rhs.selectionIsInMixture;
	mutationIsInMixture = rhs.mutationIsInMixture;
	maxGrouping = rhs.maxGrouping;
	groupList = rhs.groupList;
	mixtureAssignment = rhs.mixtureAssignment;
	categoryProbabilities = rhs.categoryProbabilities;
	traces = rhs.traces;
	numAcceptForCodonSpecificParameters = rhs.numAcceptForCodonSpecificParameters;
	std_csp = rhs.std_csp;
	covarianceMatrix = rhs.covarianceMatrix;

	noiseOffset = rhs.noiseOffset;
	noiseOffset_proposed = rhs.noiseOffset_proposed;
	std_NoiseOffset = rhs.std_NoiseOffset;
	numAcceptForNoiseOffset = rhs.numAcceptForNoiseOffset;
	return *this;
}


Parameter::~Parameter()
{
	//dtor
}





//--------------------------------------------------------------//
//---------- Initialization and Restart Functions --------------//
//--------------------------------------------------------------//


void Parameter::initParameterSet(std::vector<double> _stdDevSynthesisRate, unsigned _numMixtures,
	std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> mixtureDefinitionMatrix, bool splitSer,
    std::string _mutationSelectionState)
{
	// assign genes to mixture element
	unsigned numGenes = (unsigned)geneAssignment.size();
	mixtureAssignment.resize(numGenes, 0);

	for (unsigned i = 0u; i < numGenes; i++)
	{
		// Note: This section of code is because vectors in R are 1-indexed (i.e. for mixtureAssignment)
		//TODO:need to check index are correct, consecutive, and don't exceed numMixtures
		//possibly just use a set?
#ifndef STANDALONE
        mixtureAssignment[i] = geneAssignment[i] - 1;
		//mixtureAssignment[i + 1] = geneAssignment[i];
#else
		mixtureAssignment[i] = geneAssignment[i];
#endif
	}

	for (unsigned i = 0; i < getNumObservedPhiSets(); i++)
	{
		noiseOffset[i] = 0.1;
		noiseOffset_proposed[i] = 0.1;
		std_NoiseOffset[i] = 0.1;
		observedSynthesisNoise[i] = 0.1;
		numAcceptForNoiseOffset[i] = 0;
	}

	mutationSelectionState = _mutationSelectionState;
	// Lose one parameter if ser is split into two groups.
	// This is because one of the 2 codon set is now the reference codon
	numParam = ((splitSer) ? 40 : 41);
	numMixtures = _numMixtures;

	stdDevSynthesisRate = _stdDevSynthesisRate;
	stdDevSynthesisRate_proposed = _stdDevSynthesisRate;

	bias_stdDevSynthesisRate = 0;
	std_stdDevSynthesisRate = 0.1;

	numAcceptForStdDevSynthesisRate = 0u;
	std_csp.resize(numParam, 0.1);
	numAcceptForCodonSpecificParameters.resize(maxGrouping, 0u);
	// proposal bias and std for phi values
	bias_phi = 0;


	setNumMutationSelectionValues(_mutationSelectionState, mixtureDefinitionMatrix);
	mutationIsInMixture.resize(numMutationCategories);
	selectionIsInMixture.resize(numSelectionCategories);
	initCategoryDefinitions(_mutationSelectionState, mixtureDefinitionMatrix);

	categoryProbabilities.resize(numMixtures, 1.0/(double)numMixtures);


	//Set up vector of vectors:
	currentSynthesisRateLevel.resize(numSelectionCategories);
	proposedSynthesisRateLevel.resize(numSelectionCategories);
	numAcceptForSynthesisRate.resize(numSelectionCategories);

	std_phi.resize(numSelectionCategories);

	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> tempExpr(numGenes, 0.0);
		currentSynthesisRateLevel[i] = tempExpr;
		proposedSynthesisRateLevel[i] = tempExpr;

		std::vector<unsigned> tempAccExpr(numGenes, 0u);
		numAcceptForSynthesisRate[i] = tempAccExpr;

		std::vector<double> tempStdPhi(numGenes, 5);
		std_phi[i] = tempStdPhi;
	}
}


void Parameter::initBaseValuesFromFile(std::string filename)
{
	std::ifstream input;
	input.open(filename.c_str());
	if (input.fail())
		my_printError("Could not open file: % to initialize base values\n", filename.c_str());
	else
	{
		int cat = 0;
		std::vector<double> mat;
		std::string tmp, variableName;
		while (getline(input, tmp))
		{
			int flag;
			if (tmp[0] == '>') flag = 1;
			else if (input.eof()) flag = 2;
			else if (tmp[0] == '#') flag = 3;
			else flag = 4;

			if (flag == 1)
			{
				mat.clear();
				cat = 0;
				variableName = tmp.substr(1,tmp.size()-2);
			}
			else if (flag == 2)
			{
			}
			else if (flag == 3) //user comment, continue
			{
				continue;
			}
			else //store variable information
			{
				std::istringstream iss;
				if (variableName == "groupList")
				{
					std::string val;
					iss.str(tmp);
					while (iss >> val)
					{
						groupList.push_back(val);
					}
				}
				else if (variableName == "stdDevSynthesisRate")
				{
					stdDevSynthesisRate.resize(0);
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						stdDevSynthesisRate.push_back(val);
					}
				}
				else if (variableName == "numParam") {iss.str(tmp); iss >> numParam;}
				else if (variableName == "numMutationCategories") {iss.str(tmp); iss >> numMutationCategories;}
				else if (variableName == "numSelectionCategories") {iss.str(tmp); iss >> numSelectionCategories;}
				else if (variableName == "numMixtures") {iss.str(tmp); iss >> numMixtures;}
				else if (variableName == "mixtureAssignment")
				{
					unsigned val;
					iss.str(tmp);
					while (iss >> val)
					{
						mixtureAssignment.push_back(val);
					}
				}
				else if (variableName == "categories")
				{
					iss.str(tmp);
					mixtureDefinition K;
					iss >> K.delM;
					iss >> K.delEta;
					categories.push_back(K);
				}
				else if (variableName == "categoryProbabilities")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						categoryProbabilities.push_back(val);
					}
				}
				else if (variableName == "mutationIsInMixture")
				{
					if (tmp == "***")
					{
						mutationIsInMixture.resize(mutationIsInMixture.size() + 1);
						cat++;
					}
					else
					{
						unsigned val;
						iss.str(tmp);
						while (iss >> val)
						{
							mutationIsInMixture[cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "selectionIsInMixture")
				{
					if (tmp == "***")
					{
						selectionIsInMixture.resize(selectionIsInMixture.size() + 1);
						cat++;
					}
					else
					{
						unsigned val;
						iss.str(tmp);
						while (iss >> val)
						{
							selectionIsInMixture[cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "obsPhiSets")
				{
					iss.str(tmp);
					my_print("read\n");
					iss >> obsPhiSets;
					my_print("%",obsPhiSets);
				}
				else if (variableName == "currentSynthesisRateLevel")
				{
					if (tmp == "***")
					{
						currentSynthesisRateLevel.resize(currentSynthesisRateLevel.size() + 1);
						cat++;
					}
					else
					{
						double val;
						iss.str(tmp);
						while (iss >> val)
						{
							currentSynthesisRateLevel[cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "std_stdDevSynthesisRate")
				{
					iss.str(tmp);
					iss >> std_stdDevSynthesisRate;
				}
				else if (variableName == "std_phi")
				{
					if (tmp == "***")
					{
						std_phi.resize(std_phi.size() + 1);
						cat++;
					}
					iss.str(tmp);
					double val;
					while (iss >> val)
					{
						std_phi[cat - 1].push_back(val);
					}
				}
				else if (variableName == "noiseOffset")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						noiseOffset.push_back(val);
						noiseOffset_proposed.push_back(val);
					}
				}
				else if (variableName == "observedSynthesisNoise")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						observedSynthesisNoise.push_back(val);
					}
				}
				else if (variableName == "std_NoiseOffset")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						std_NoiseOffset.push_back(val);
					}
				}
			}
		}

		input.close();

		//initialize all the default Parameter values now.
		stdDevSynthesisRate_proposed = stdDevSynthesisRate;
		numAcceptForStdDevSynthesisRate = 0u;
		bias_stdDevSynthesisRate = 0;
		bias_phi = 0;
		//obsPhiSets = 0;
		numAcceptForNoiseOffset.resize(obsPhiSets, 0);
		numAcceptForSynthesisRate.resize(numSelectionCategories);
		proposedSynthesisRateLevel.resize(numSelectionCategories);
		for (unsigned i = 0; i < numSelectionCategories; i++)
		{
			proposedSynthesisRateLevel[i] = currentSynthesisRateLevel[i];
			std::vector <unsigned> tmp2(currentSynthesisRateLevel[i].size(), 0u);
			numAcceptForSynthesisRate[i] = tmp2;
		}
	}
}


void Parameter::writeBasicRestartFile(std::string filename)
{
	my_print("Begin writing restart file\n");

	std::ofstream out;
	std::string output = "";
	std::ostringstream oss;
	unsigned i, j;

	out.open(filename.c_str());
	if (out.fail())
		my_printError("Error: Could not open restart file % for writing\n", filename.c_str());
	else
	{
		oss << ">groupList:\n";
		for (i = 0; i < groupList.size(); i++)
		{
			oss << groupList[i];
			if ((i + 1) % 10 == 0) oss << "\n";
			else oss << " ";
		}
		if (i % 10 != 0) oss << "\n";
		oss << ">stdDevSynthesisRate:\n";
		for (i = 0; i < stdDevSynthesisRate.size(); i++)
		{
			oss << stdDevSynthesisRate[i];
			if ((i + 1) % 10 == 0) oss << "\n";
			else oss <<" ";
		}
		if (i % 10 != 0) oss << "\n";
		oss << ">numParam:\n" << numParam << "\n";
		oss << ">numMixtures:\n" << numMixtures << "\n";
		oss << ">std_stdDevSynthesisRate:\n" << std_stdDevSynthesisRate << "\n";
		//TODO: maybe clear the buffer
		oss << ">std_phi:\n";
		for (i = 0; i < std_phi.size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < std_phi[i].size(); j++)
			{
				oss << std_phi[i][j];
				if ((j + 1) % 10 == 0) oss << "\n";
				else oss <<" ";
			}
			if (j % 10 != 0) oss <<"\n";
		}
		oss << ">categories:\n";
		for (i = 0; i < categories.size(); i++)
		{
			oss << categories[i].delM << " " << categories[i].delEta << "\n";
		}

		oss << ">mixtureAssignment:\n";
		for (i = 0; i < mixtureAssignment.size(); i++)
		{
			oss << mixtureAssignment[i];
			if ((i + 1) % 50 == 0) oss <<"\n";
			else oss <<" ";
		}
		if (i % 50 != 0) oss <<"\n";
		oss << ">numMutationCategories:\n" << numMutationCategories << "\n";
		oss << ">numSelectionCategories:\n" << numSelectionCategories << "\n";

		oss << ">categoryProbabilities:\n";
		for (i = 0; i < categoryProbabilities.size(); i++)
		{
			oss << categoryProbabilities[i];
			if ((i + 1) % 10 == 0) oss << "\n";
			else oss <<" ";
		}
		if (i % 10 != 0) oss <<"\n";

		oss << ">selectionIsInMixture:\n";
		for (i = 0; i < selectionIsInMixture.size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < selectionIsInMixture[i].size(); j++)
			{
				oss << selectionIsInMixture[i][j] <<" ";
			}
			oss << "\n";
		}

		oss << ">mutationIsInMixture:\n";
		for (i = 0; i < mutationIsInMixture.size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < mutationIsInMixture[i].size(); j++)
			{
				oss << mutationIsInMixture[i][j] << " ";
			}
			oss << "\n";
		}
		oss << ">obsPhiSets:\n" << obsPhiSets << "\n";
		oss << ">currentSynthesisRateLevel:\n";
		for (i = 0; i < currentSynthesisRateLevel.size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentSynthesisRateLevel[i].size(); j++)
			{
				oss << currentSynthesisRateLevel[i][j];
				if ((j + 1) % 10 == 0) oss << "\n";
				else oss <<" ";
			}
			if (j % 10 != 0) oss << "\n";
		}
		oss << ">noiseOffset:\n";
		for (j = 0; j < noiseOffset.size(); j++)
		{
			oss << noiseOffset[j];
			if ((j + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (j % 10 != 0)
			oss << "\n";

		oss << ">observedSynthesisNoise:\n";
		for (j = 0; j < observedSynthesisNoise.size(); j++)
		{
			oss << observedSynthesisNoise[j];
			if ((j + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (j % 10 != 0)
			oss << "\n";
		oss << ">std_NoiseOffset:\n";
		for (j = 0; j < std_NoiseOffset.size(); j++)
		{
			oss << std_NoiseOffset[j];
			if ((j + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (j % 10 != 0)
			oss << "\n";
	}
	my_print("End writing restart file\n");

	output += oss.str();
	out << output;
	out.close();
}


void Parameter::initCategoryDefinitions(std::string _mutationSelectionState,
										std::vector<std::vector<unsigned>> mixtureDefinitionMatrix)
{
	std::set<unsigned> delMCounter;
	std::set<unsigned> delEtaCounter;

	for (unsigned i = 0u; i < numMixtures; i++)
	{
		categories.push_back(mixtureDefinition()); //push a blank mixtureDefinition on the vector, then alter.
		if (!mixtureDefinitionMatrix.empty())
		{
			categories[i].delM = mixtureDefinitionMatrix[i][0] - 1;
			categories[i].delEta = mixtureDefinitionMatrix[i][1] - 1; //need check for negative and consecutive checks
			mutationIsInMixture[mixtureDefinitionMatrix[i][0] - 1].push_back(i);
			selectionIsInMixture[mixtureDefinitionMatrix[i][1] - 1].push_back(i);
		}
		else if (_mutationSelectionState == selectionShared)
		{
			categories[i].delM = i;
			categories[i].delEta = 0;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[0].push_back(i);
		}
		else if (_mutationSelectionState == mutationShared)
		{
			categories[i].delM = 0;
			categories[i].delEta = i;
			mutationIsInMixture[0].push_back(i);
			selectionIsInMixture[i].push_back(i);
		}
		else //assuming the default of allUnique
		{
			categories[i].delM = i;
			categories[i].delEta = i;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[i].push_back(i);
		}
		delMCounter.insert(categories[i].delM);
		delEtaCounter.insert(categories[i].delEta);
	}

	//sets allow only the unique numbers to be added.
	//at the end, the size of the set is equal to the number
	//of unique categories.
}


/* InitializeSynthesisRate (using SCUO per genome) (RCPP EXPOSED VIA WRAPPER)
 * Arguments: //TODO
*/
void Parameter::InitializeSynthesisRate(Genome& genome, double sd_phi)
{
	unsigned genomeSize = genome.getGenomeSize();
	double* SCUOValues = new double[genomeSize]();
	double* expression = new double[genomeSize]();
	int* index = new int[genomeSize]();


	for (unsigned i = 0u; i < genomeSize; i++)
	{
		index[i] = i;
		//This used to be maxGrouping instead of 22, but PA model will not work that way
		SCUOValues[i] = calculateSCUO( genome.getGene(i));
		expression[i] = Parameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
	}

	quickSortPair(SCUOValues, index, 0, genomeSize);
	std::sort(expression, expression + genomeSize);

	for (unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for (unsigned j = 0u; j < genomeSize; j++)
		{
			currentSynthesisRateLevel[category][index[j]] = expression[j];
			std_phi[category][j] = 0.1;
			numAcceptForSynthesisRate[category][j] = 0u;
		}
	}

	delete [] SCUOValues;
	delete [] expression;
	delete [] index;
}


/* InitializeSynthesisRate (by random) (RCPP EXPOSED VIA WRAPPER)
 * Arguments: //TODO
*/
void Parameter::InitializeSynthesisRate(double sd_phi)
{	unsigned genomeSize = (unsigned)currentSynthesisRateLevel[0].size();
	for (unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for (unsigned i = 0u; i < genomeSize; i++)
		{
			currentSynthesisRateLevel[category][i] = Parameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
			std_phi[category][i] = 0.1;
			numAcceptForSynthesisRate[category][i] = 0u;
		}
	}
}


/* InitializeSynthesisRate (by list) (RCPP EXPOSED VIA WRAPPER)
 * Arguments: //TODO
*/
void Parameter::InitializeSynthesisRate(std::vector<double> expression)
{
	unsigned numGenes = (unsigned)currentSynthesisRateLevel[0].size();
	for (unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for (unsigned i = 0u; i < numGenes; i++)
		{
			currentSynthesisRateLevel[category][i] = expression[i];
			std_phi[category][i] = 0.1;
			numAcceptForSynthesisRate[category][i] = 0u;
		}
	}
}

//' @name readPhiValue
//' @title Read synthesis rate values from file. File should be two column file <gene_id,phi> and is expected to have a header row
//' @param filename name of file to be read
std::vector <double> Parameter::readPhiValues(std::string filename)
{
	std::size_t pos;
	std::ifstream currentFile;
	std::string tmpString;
	std::vector <double> RV;

	currentFile.open(filename);
	if (currentFile.fail())
		my_printError("Error opening file %\n", filename.c_str());
	else
	{
		currentFile >> tmpString; //trash the first line, no info given.
		while (currentFile >> tmpString)
		{
			pos = tmpString.find(',');
			if (pos != std::string::npos)
			{
				std::string val = tmpString.substr(pos + 1);
				//RV.push_back(std::stod(val));
				RV.push_back(std::atof(val.c_str()));
			}
		}
	}
	return RV;
}


//----------------------------------------------------------------------//
//-------------------------- Prior functions ---------------------------//
//----------------------------------------------------------------------//


double Parameter::getCodonSpecificPriorStdDev(unsigned paramType)
{
	return codonSpecificPrior[paramType];
}


//----------------------------------------------------------------------//
//---------- Mixture Definition Matrix and Category Functions ----------//
//----------------------------------------------------------------------//


void Parameter::setNumMutationSelectionValues(std::string _mutationSelectionState,
											  std::vector<std::vector<unsigned>> mixtureDefinitionMatrix)
{
	if (!mixtureDefinitionMatrix.empty())
	{
		//sets allow only the unique numbers to be added.
		//at the end, the size of the set is equal to the number
		//of unique categories.
		std::set<unsigned> delMCounter;
		std::set<unsigned> delEtaCounter;
		for (unsigned i = 0u; i < numMixtures; i++)
		{
			delMCounter.insert(mixtureDefinitionMatrix[i][0] - 1);
			delEtaCounter.insert(mixtureDefinitionMatrix[i][1] - 1);


		}
		numMutationCategories = (unsigned)delMCounter.size();
		numSelectionCategories = (unsigned)delEtaCounter.size();
	}
	else if (_mutationSelectionState == selectionShared)
	{
		numMutationCategories = numMixtures;
		numSelectionCategories = 1u;
	}
	else if (_mutationSelectionState == mutationShared)
	{
		numMutationCategories = 1u;
		numSelectionCategories = numMixtures;
	}
	else //assuming the default of allUnique
	{
		numMutationCategories = numMixtures;
		numSelectionCategories = numMixtures;
	}
}


void Parameter::printMixtureDefinitionMatrix()
{
	for (unsigned i = 0u; i < numMixtures; i++)
		my_print("%\t%\n", categories[i].delM, categories[i].delEta);
}


/* getCategoryProbability (NOT EXPOSED)
 * Arguments: A number representing a mixture element
 * Returns the category probability of the mixture element given.
*/
double Parameter::getCategoryProbability(unsigned mixtureElement)
{
	return categoryProbabilities[mixtureElement];
}


/* setCategoryProbability (NOT EXPOSED)
 * Arguments: A number representing a mixture element, a double representing a probability
 * Sets the probability for the category of the mixture element to the value given.
*/
void Parameter::setCategoryProbability(unsigned mixtureElement, double value)
{
	categoryProbabilities[mixtureElement] = value;
}


unsigned Parameter::getNumMutationCategories()
{
	return numMutationCategories;
}


unsigned Parameter::getNumSelectionCategories()
{
	return numSelectionCategories;
}


unsigned Parameter::getNumSynthesisRateCategories()
{
	return numSelectionCategories;
}

//Used to get alpha category in PA and PANSE
unsigned Parameter::getMutationCategory(unsigned mixtureElement)
{
	return categories[mixtureElement].delM;
}


//TODO: Considering renaming so that the two sides of the function match more than they currently do.
/* Note 1) -- on getSelectionCategory and getSynthesisRateCategory
 * These two functions are technically the same for readability.
 * Selection and synthesis rate are directly related even if they are not known
 * and thus are represented by the same variable. By splitting this
 * into selection and synthesis, we avoid confusion when otherwise
 * we may ask why one is used in the place of another.
*/

/* getSelectionCategory (RCPP EXPOSED VIA WRAPPER)
 * Arguments: A number representing a mixture element
 * Returns the selection category of the mixture element chosen.
 * See Note 1) above.
 * Wrapped by getSelectionCategoryForMixture on the R-side.
 */
//Used to get lambda category in PA and PANSE
unsigned Parameter::getSelectionCategory(unsigned mixtureElement)
{
	return categories[mixtureElement].delEta;
}


//TODO: Considering renaming so that the two sides of the function match more than they currently do.
/* getSynthesisRateCategory (RCPP EXPOSED VIA WRAPPER)
 * Arguments: A number representing a mixture element
 * Returns the synthesis rate category of the mixture element chosen.
 * See Note 1) above.
 * Wrapped by getSynthesisRateCategoryForMixture on the R-side.
 */
unsigned Parameter::getSynthesisRateCategory(unsigned mixtureElement)
{
	return categories[mixtureElement].delEta;
}


std::vector<unsigned> Parameter::getMixtureElementsOfMutationCategory(unsigned category)
{
	return mutationIsInMixture[category];
}


std::vector<unsigned> Parameter::getMixtureElementsOfSelectionCategory(unsigned category)
{
	return selectionIsInMixture[category];
}


std::string Parameter::getMutationSelectionState()
{
	return mutationSelectionState;
}


/* getNumAcceptForCspForIndex (NOT EXPOSED)
 * Arguments: index of numAcceptForCodonSpecificParameters to be returned
 * Returns the numAcceptForCodonSpecificParameters at the index given.
 * Note: Used in unit testing only.
*/
unsigned Parameter::getNumAcceptForCspForIndex(unsigned i)
{
	return numAcceptForCodonSpecificParameters[i];
}





// -------------------------------------------//
// ---------- Group List Functions -----------//
// -------------------------------------------//


//TODO: Hollis: Verify or implement the Group List such that the MCMC will now only run on the subset of codons set by
// this function.

/* setGroupList (RCPP EXPOSED)
 * Arguments: vector of strings representing a group list
 * Sets the group list to the argument after clearing the group list, adding elements only if they have no errors.
*/

//' @name setGroupList
//' @title Set amino acids (ROC, FONSE) or codons (PA, PANSE) for which parameters will be estimated. Note that non-default groupLists are still in beta testing and should be used with caution.
//' @param List of strings epresenting groups for parameters to be estimated. Should be one letter amino acid (ROC, FONSE) or list of sense codons (PA, PANSE). 
void Parameter::setGroupList(std::vector <std::string> gl)
{
	groupList.clear();
	for (unsigned i = 0; i < gl.size(); i++)
	{
		groupList.push_back(gl[i]);
	}
}


/* getGrouping (NOT EXPOSED)
 * Arguments: index of a group list element to be returned
 * Returns the group list element at the given index.
*/
std::string Parameter::getGrouping(unsigned index)
{
	return groupList[index];
}


/* getGroupList (RCPP EXPOSED)
 * Arguments: None
 * Returns the group list as a vector of strings.
*/

//' @name getGroupList
//' @title Get amino acids (ROC, FONSE) or codons (PA, PANSE) for which parameters will be estimated
std::vector<std::string> Parameter::getGroupList()
{
	return groupList;
}


/* getGroupListSize (NOT EXPOSED)
 * Arguments: None
 * Returns the size of the group list.
*/
unsigned Parameter::getGroupListSize()
{
	return (unsigned) groupList.size();
}





//----------------------------------------------------//
//---------- stdDevSynthesisRate Functions -----------//
//----------------------------------------------------//

void Parameter::fixStdDevSynthesis()
{
	fix_stdDevSynthesis = true;
}


double Parameter::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return (proposed ? stdDevSynthesisRate_proposed[selectionCategory] : stdDevSynthesisRate[selectionCategory]);
}


void Parameter::proposeStdDevSynthesisRate()
{
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{	
		if (!fix_stdDevSynthesis)
		{
			stdDevSynthesisRate_proposed[i] = std::exp(randNorm(std::log(stdDevSynthesisRate[i]), std_stdDevSynthesisRate));
		}
		else
		{
			stdDevSynthesisRate_proposed[i] = stdDevSynthesisRate[i];
		}
	}
}


void Parameter::setStdDevSynthesisRate(double _stdDevSynthesisRate, unsigned selectionCategory)
{
	stdDevSynthesisRate[selectionCategory] = _stdDevSynthesisRate;
}


double Parameter::getCurrentStdDevSynthesisRateProposalWidth()
{
	return std_stdDevSynthesisRate;
}


/* getNumAcceptForStdDevSynthesisRate (NOT EXPOSED)
 * Arguments: None
 * Returns the numAcceptForStdDevSynthesisRate.
 * Note: Used in unit testing only.
*/
unsigned Parameter::getNumAcceptForStdDevSynthesisRate()
{
	return numAcceptForStdDevSynthesisRate;
}


void Parameter::updateStdDevSynthesisRate()
{
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		stdDevSynthesisRate[i] = stdDevSynthesisRate_proposed[i];
	}
	numAcceptForStdDevSynthesisRate++;
}


/* getStdCspForIndex (NOT EXPOSED)
 * Arguments: index of standard deviation (proposal width) of the codon-specific parameter to be returned
 * Returns the standard deviation (proposal width) of the codon-specific parameter at the index given.
 * Note: Used in unit testing only.
*/
double Parameter::getStdCspForIndex(unsigned i)
{
	return std_csp[i];
}





//-----------------------------------------------//
//---------- Synthesis Rate Functions -----------//
//-----------------------------------------------//

double Parameter::getSynthesisRate(unsigned geneIndex, unsigned mixtureElement, bool proposed)
{
	unsigned category = getSelectionCategory(mixtureElement);
 	return  (proposed ? proposedSynthesisRateLevel[category][geneIndex] : currentSynthesisRateLevel[category][geneIndex]);
}


/* Note 2) -- on getCurrentSynthesisRateProposalWidth and getSynthesisRateProposalWidth
 * These two functions should perform the same action if properly used.
 * Similar to Note 1), these functions are based on how
 * synthesis rate category and selection category are directly related
 * but for readability two separated functions are created.
*/

/* getCurrentSynthesisRateProposalWidth (NOT EXPOSED)
 * Arguments: index of a gene in the genome, number representing the selected category
 * Returns the current synthesis rate proposal width of the category of the mixture element for the gene indexed.
 * See Note 2) above.
*/
double Parameter::getCurrentSynthesisRateProposalWidth(unsigned expressionCategory, unsigned geneIndex)
{
	return std_phi[expressionCategory][geneIndex];
}


/* getSynthesisRateProposalWidth (NOT EXPOSED)
 * Arguments: index of a gene in the genome, number representing a mixture element
 * Returns the synthesis rate proposal width of the category of the mixture element for the gene indexed.
 * See Note 2) above.
*/
double Parameter::getSynthesisRateProposalWidth(unsigned geneIndex, unsigned mixtureElement)
{
	unsigned category = getSelectionCategory(mixtureElement);
	return std_phi[category][geneIndex];
}


void Parameter::proposeSynthesisRateLevels()
{
	unsigned numSynthesisRateLevels = (unsigned) currentSynthesisRateLevel[0].size();
	for (unsigned category = 0; category < numSelectionCategories; category++)
	{
		for (unsigned i = 0u; i < numSynthesisRateLevels; i++)
		{
			// avoid adjusting probabilities for asymmetry of distribution
			proposedSynthesisRateLevel[category][i] = std::exp( randNorm( std::log(currentSynthesisRateLevel[category][i]),
																		  std_phi[category][i]) );
		}
	}
}


void Parameter::setSynthesisRate(double phi, unsigned geneIndex, unsigned mixtureElement)
{
	unsigned category = getSelectionCategory(mixtureElement);
	currentSynthesisRateLevel[category][geneIndex] = phi;
}


void Parameter::updateSynthesisRate(unsigned geneIndex)
{
	for (unsigned category = 0; category < numSelectionCategories; category++)
	{
		numAcceptForSynthesisRate[category][geneIndex]++;
		currentSynthesisRateLevel[category][geneIndex] = proposedSynthesisRateLevel[category][geneIndex];
	}
}


void Parameter::updateSynthesisRate(unsigned geneIndex, unsigned mixtureElement)
{
	unsigned category = getSelectionCategory(mixtureElement);
	numAcceptForSynthesisRate[category][geneIndex]++;
	currentSynthesisRateLevel[category][geneIndex] = proposedSynthesisRateLevel[category][geneIndex];
}


/* getNumAcceptForSynthesisRate (NOT EXPOSED)
 * Arguments: index of a gene in the genome, number representing the selected category
 * Returns the numAcceptForSynthesisRate of the category of the mixture element for the gene indexed.
 * Note: Used in unit testing only.
*/
unsigned Parameter::getNumAcceptForSynthesisRate(unsigned expressionCategory, unsigned geneIndex)
{
	return numAcceptForSynthesisRate[expressionCategory][geneIndex];
}





//------------------------------------------//
//---------- Iteration Functions -----------//
//------------------------------------------//


/* setLastIteration (RCPP EXPOSED)
 * Arguments: None
 * Returns the last iteration.
*/
unsigned Parameter::getLastIteration()
{
	return lastIteration;
}


/* setLastIteration (RCPP EXPOSED)
 * Arguments: unsigned value representing an iteration
 * Sets the last iteration to the argument given.
*/
void Parameter::setLastIteration(unsigned iteration)
{
	lastIteration = iteration;
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


unsigned Parameter::getNumParam()
{
	return numParam;
}


unsigned Parameter::getNumMixtureElements()
{
	return numMixtures;
}


/* getNumObservedPhiSets (NOT EXPOSED)
 * Arguments: None
 * Returns the observed number of phi sets.
*/
unsigned Parameter::getNumObservedPhiSets()
{
	return obsPhiSets;
}


/* setNumObservedPhiSets (NOT EXPOSED)
 * Arguments: unsigned value representing a new number of phi set groupings
 * Sets the observed number of phi sets to the argument given.
*/
// void Parameter::setNumObservedPhiSets(unsigned _phiGroupings)
// {
// 	obsPhiSets = _phiGroupings;
// }

void Parameter::setNumObservedPhiSets(unsigned _phiGroupings)
{
	obsPhiSets = _phiGroupings;
	noiseOffset.resize(obsPhiSets, 0.1);
	noiseOffset_proposed.resize(obsPhiSets, 0.1);
	std_NoiseOffset.resize(obsPhiSets, 0.1);
	numAcceptForNoiseOffset.resize(obsPhiSets, 0);
	observedSynthesisNoise.resize(obsPhiSets, 1.0);
}


void Parameter::setMixtureAssignment(unsigned gene, unsigned value)
{
	mixtureAssignment[gene] = value;
}


unsigned Parameter::getMixtureAssignment(unsigned gene)
{
	return mixtureAssignment[gene];
}


std::vector<std::vector<double>> Parameter::calculateSelectionCoefficients(unsigned sample)
{
	unsigned numGenes = (unsigned)mixtureAssignment.size();
	std::vector<std::vector<double>> selectionCoefficients;
	selectionCoefficients.resize(numGenes, std::vector<double> (61, 0));
	unsigned numGroupings = getGroupListSize();

	for (unsigned i = 0; i < numGenes; i++)
	{
		unsigned codon_index = 0;
		double phi = getSynthesisRatePosteriorMean(sample, i, false);

		for (unsigned j = 0; j < numGroupings; j++)
		{
			std::string aa = getGrouping(j);
			unsigned aaStart, aaEnd;
			SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
			std::vector<double> tmp;
			double minValue = 0.0;
			for (unsigned k = aaStart; k < aaEnd; k++)
			{
				std::string codon = SequenceSummary::codonArrayParameter[k];
				double x = getCodonSpecificPosteriorMean(i, sample, codon, 1, true, true);
				tmp.push_back(x);
				minValue = (x < minValue) ? x : minValue;
			}
			tmp.push_back(0.0);
			for (unsigned k = 0; k < tmp.size(); k++, codon_index++)
			{
				tmp[k] -= minValue;
				selectionCoefficients[i][codon_index] = -(phi * tmp[k]);
			}
		}
	}
	return selectionCoefficients;
}





//--------------------------------------//
//---------- Trace Functions -----------//
//--------------------------------------//


//' @name getTraceObject
//' @title Get Trace object stored by a Parameter object. Useful for plotting certain parameter traces.
Trace& Parameter::getTraceObject()
{
	return traces;
}


void Parameter::setTraceObject(Trace _trace)
{
	traces = _trace;
}

void Parameter::updateObservedSynthesisNoiseTraces(unsigned sample)
{
	for (unsigned i = 0; i < observedSynthesisNoise.size(); i++)
	{
		traces.updateObservedSynthesisNoiseTrace(i, sample, observedSynthesisNoise[i]);
	}
}


void Parameter::updateNoiseOffsetTraces(unsigned sample)
{
	for (unsigned i = 0; i < noiseOffset.size(); i++)
	{
		traces.updateSynthesisOffsetTrace(i, sample, noiseOffset[i]);
	}
}



void Parameter::updateStdDevSynthesisRateTrace(unsigned sample)
{
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		traces.updateStdDevSynthesisRateTrace(sample, stdDevSynthesisRate[i], i);
	}
}


void Parameter::updateSynthesisRateTrace(unsigned sample, unsigned geneIndex)
{
	traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel);
}


void Parameter::updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex)
{
	traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]);
}


void Parameter::updateMixtureProbabilitiesTrace(unsigned samples)
{
	traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities);
}


//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


void Parameter::adaptNoiseOffsetProposalWidth(unsigned adaptationWidth, bool adapt)
{
	for (unsigned i = 0; i < getNumObservedPhiSets(); i++)
	{
		double acceptanceLevel = numAcceptForNoiseOffset[i] / (double)adaptationWidth;
		traces.updateSynthesisOffsetAcceptanceRateTrace(i, acceptanceLevel);
		if (adapt)
		{
			if (acceptanceLevel < 0.2)
				std_NoiseOffset[i] *= 0.8;
			if (acceptanceLevel > 0.3)
				std_NoiseOffset[i] *= 1.2;

			numAcceptForNoiseOffset[i] = 0u;
		}
	}
}




//Adjust s_phi proposal distribution
void Parameter::adaptStdDevSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt)
{
	double acceptanceLevel = (double)numAcceptForStdDevSynthesisRate / (double)adaptationWidth;
	traces.updateStdDevSynthesisRateAcceptanceRateTrace(acceptanceLevel);
	if (adapt)
	{
		if (acceptanceLevel < 0.2)
			std_stdDevSynthesisRate *= 0.8;

		if (acceptanceLevel > 0.3)
			std_stdDevSynthesisRate *= 1.2;
	}
	numAcceptForStdDevSynthesisRate = 0u;
}


void Parameter::adaptSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt)
{
	unsigned acceptanceUnder = 0u;
	unsigned acceptanceOver = 0u;

	// mikeg: variables below should likely be made global or added to parameter so that the same criteria is used in all adaptive proposal width routines.
	// From below
	//Gelman BDA 3rd Edition suggests a target acceptance rate of 0.23
	// for high dimensional problems
	// We are dealing with single parameters here, however.

	double acceptanceTargetLow = 0.225; //below this value factor adjustment is applied
	double acceptanceTargetHigh = 0.325; //above this value factor adjustment is applied
	double factorCriteriaLow;
	double factorCriteriaHigh;
	double adjustFactorLow = 0.8; //factor by which to reduce proposal widths
	double adjustFactorHigh = 1.3; //factor by which to increase proposal widths
	double adjustFactor = 1.0; //variable assigned value of either adjustFactorLow or adjustFactorHigh depending on acceptance rate

	factorCriteriaLow = acceptanceTargetLow;
	factorCriteriaHigh = acceptanceTargetHigh;

	for (unsigned cat = 0u; cat < numSelectionCategories; cat++)
	{
		unsigned numGenes = (unsigned)numAcceptForSynthesisRate[cat].size();
		for (unsigned i = 0; i < numGenes; i++)
		{
			double acceptanceLevel = (double)numAcceptForSynthesisRate[cat][i] / (double)adaptationWidth;
			traces.updateSynthesisRateAcceptanceRateTrace(cat, i, acceptanceLevel);

			//Evaluate acceptance rates relative to target region
			if (acceptanceLevel < acceptanceTargetLow) acceptanceUnder++;
			else if (acceptanceLevel > acceptanceTargetHigh) acceptanceOver++;

			if (adapt)
			{
				if (acceptanceLevel < acceptanceTargetLow)
				{
					std_phi[cat][i] *= adjustFactorLow;

				}
				if (acceptanceLevel > acceptanceTargetHigh)
				{
					std_phi[cat][i] *= adjustFactorHigh;

				}
			}
			numAcceptForSynthesisRate[cat][i] = 0u;
		}
	}

	my_print("Acceptance rate for synthesis rate:\n");
	my_print("Target range: %-% \n", acceptanceTargetLow, acceptanceTargetHigh );
	my_print("Adjustment range: < % or > % \n", factorCriteriaLow, factorCriteriaHigh );
	my_print("\t acceptance rates below lower target of %: %\n", acceptanceTargetLow, acceptanceUnder);
	my_print("\t acceptance rate above upper target of %: %\n", acceptanceTargetHigh, acceptanceOver);
}


void Parameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt)
{
  //Gelman BDA 3rd Edition suggests a target acceptance rate of 0.23
  // for high dimensional problems
  // For CSP the combined selection and mutation dimensions range from 2 to 10
  //Adjust proposal variance to try and get within this range

  unsigned acceptanceUnder = 0u;
  unsigned acceptanceOver = 0u;

  double acceptanceTargetLow = 0.225; //below this value weighted sum adjustment is applied, was 0.2
  double acceptanceTargetHigh = 0.325;///above this value weighted sum adjustment is applied, was 0.3
  double diffFactorAdjust = 0.05; //sets when multiplication factor adjustment is applied, was 0.1 and 0.0, respectively
  double factorCriteriaLow;
  double factorCriteriaHigh;
  double adjustFactorLow = 0.8; //factor by which to reduce proposal widths
  double adjustFactorHigh = 1.2; //factor by which to increase proposal widths
  double adjustFactor = 1.0; //variable assigned value of either adjustFactorLow or adjustFactorHigh

  factorCriteriaLow = acceptanceTargetLow - diffFactorAdjust;  //below this value weighted sum and factor adjustments are applied
  factorCriteriaHigh = acceptanceTargetHigh + diffFactorAdjust;  //above this value weighted sum and factor adjustments are applied

  adaptiveStepPrev = adaptiveStepCurr;
  adaptiveStepCurr = lastIteration;
  unsigned samples = adaptiveStepCurr - adaptiveStepPrev;

  my_print("Acceptance rates for Codon Specific Parameters\n");
  my_print("Target range: %-% \n", factorCriteriaLow, factorCriteriaHigh );
  my_print("Adjustment range: < % or > % \n", acceptanceTargetLow, acceptanceTargetHigh );
  my_print("\tAA\tAcc.Rat\n"); //Prop.Width\n";

  for (unsigned i = 0; i < groupList.size(); i++) //cycle through all of the aa
  {
  	std::string aa = groupList[i];
    unsigned aaIndex = SequenceSummary::AAToAAIndex(aa);
    double acceptanceLevel = (double)numAcceptForCodonSpecificParameters[aaIndex] / (double)adaptationWidth;

    my_print("\t%:\t%\n", aa.c_str(), acceptanceLevel);

    traces.updateCodonSpecificAcceptanceRateTrace(aaIndex, acceptanceLevel);

    unsigned aaStart, aaEnd;
    SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);

      //Evaluate current acceptance ratio  performance
    if (acceptanceLevel < factorCriteriaLow) acceptanceUnder++;
    else if (acceptanceLevel > factorCriteriaHigh) acceptanceOver++;

    if (adapt)
	{
		if( (acceptanceLevel < acceptanceTargetLow) || (acceptanceLevel > acceptanceTargetHigh) )// adjust proposal width
	  	{
	  	  //Update cov matrix based on previous window to improve efficiency of sampling
	      CovarianceMatrix covcurr(covarianceMatrix[aaIndex].getNumVariates());
	      covcurr.calculateSampleCovariance(*traces.getCodonSpecificParameterTrace(), aa, samples, adaptiveStepCurr);
	      CovarianceMatrix covprev = covarianceMatrix[aaIndex];
	      covprev = (covprev*0.6);
	      covcurr = (covcurr*0.4);
	      covarianceMatrix[aaIndex] = covprev + covcurr;
	      //replace cov matrix based on previous window
	      //The is approach was commented out and above code uncommented to replace it in commit ec63bb21a1e9 (2016).  Should remove
	      //covarianceMatrix[aaIndex].calculateSampleCovariance(*traces.getCodonSpecificParameterTrace(), aa, samples, adaptiveStepCurr);

	      // define adjustFactor
	      if (acceptanceLevel < factorCriteriaLow)
	      {
		  	adjustFactor = adjustFactorLow;
		  }
	      else if(acceptanceLevel > factorCriteriaHigh)
	      {
		  	adjustFactor = adjustFactorHigh;
		  }
	      else //Don't adjust
	      {
		  	adjustFactor = 1.0;
		  }

		 if( adjustFactor != 1.0 )
		 {
		 
		    //Adjust proposal width for codon specific parameters
		   //	std_csp[k] *= adjustFactor;

		    //Adjust widths if using cov matrix
		   	covarianceMatrix[aaIndex] *= adjustFactor;
	  	   
		 }

		      //Decomposing of cov matrix to convert iid samples to covarying samples using matrix decomposition
		      //The decomposed matrix is used in the proposal of new samples
		 covarianceMatrix[aaIndex].choleskyDecomposition();
		}// end adjust loop
	} // end if(adapt)

	numAcceptForCodonSpecificParameters[aaIndex] = 0u;
   }
}


//------------------------------------------------------//
//---------- observedSynthesisNoise Functions ----------//
//------------------------------------------------------//


double Parameter::getObservedSynthesisNoise(unsigned index)
{
	return observedSynthesisNoise[index];
}


void Parameter::setObservedSynthesisNoise(unsigned index, double se)
{
	observedSynthesisNoise[index] = se;
}

void Parameter::setInitialValuesForSepsilon(std::vector<double> seps)
{
	if (seps.size() == observedSynthesisNoise.size())
	{
		for (unsigned i = 0; i < observedSynthesisNoise.size(); i++)
		{
			observedSynthesisNoise[i] = seps[i];
		}
	}
	else
	{
		my_printError("Parameter::setInitialValuesForSepsilon number of initial values (%) does not match number of expression sets (%)",
					  seps.size(), observedSynthesisNoise.size());
	}
}





//-------------------------------------------//
//---------- noiseOffset Functions ----------//
//-------------------------------------------//


double Parameter::getNoiseOffset(unsigned index, bool proposed)
{
	return (proposed ? noiseOffset_proposed[index] : noiseOffset[index]);
}


double Parameter::getCurrentNoiseOffsetProposalWidth(unsigned index)
{
	return std_NoiseOffset[index];
}


void Parameter::proposeNoiseOffset()
{
	for (unsigned i = 0; i < getNumObservedPhiSets(); i++)
	{
		noiseOffset_proposed[i] = randNorm(noiseOffset[i], std_NoiseOffset[i]);
	}
}


void Parameter::setNoiseOffset(unsigned index, double _noiseOffset)
{
	noiseOffset[index] = _noiseOffset;
}


void Parameter::updateNoiseOffset(unsigned index)
{
	noiseOffset[index] = noiseOffset_proposed[index];
	numAcceptForNoiseOffset[index]++;
}



//------------------------------------------------------------------//
//---------- Posterior, Variance, and Estimates Functions ----------//
//------------------------------------------------------------------//


double Parameter::getNoiseOffsetPosteriorMean(unsigned index, unsigned samples)
{
	double posteriorMean = 0.0;
	std::vector<double> NoiseOffsetTrace = traces.getSynthesisOffsetTrace(index);
	unsigned traceLength = lastIteration;

	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getNoiseOffsetPosteriorMean throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n", samples, traceLength);

		samples = traceLength;
	}
	unsigned start = traceLength - samples;

	for (unsigned i = start; i < traceLength; i++)
		posteriorMean += NoiseOffsetTrace[i];

	return posteriorMean / (double)samples;
}

//' @name getNoiseOffsetVariance
//' @title Calculate variance of noise offset parameter used when fitting model with empirical estimates of synthesis rates (ie. withPhi fits)
//' @param index mixture index to use. Should be number between 0 and n-1, where n is number of mixtures
//' @param samples number of samples over which to calculate variance
//' @param unbiased If TRUE, should calculate variance using unbiased (N-1). Otherwise, used biased (N) correction
//' @return returns variance for noise offset

double Parameter::getNoiseOffsetVariance(unsigned index, unsigned samples, bool unbiased)
{
	std::vector<double> NoiseOffsetTrace = traces.getSynthesisOffsetTrace(index);
	unsigned traceLength = lastIteration;
	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getNoiseOffsetVariance throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n", samples, traceLength);

		samples = traceLength;
	}
	double posteriorMean = getNoiseOffsetPosteriorMean(index, samples);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		double difference = NoiseOffsetTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
	return normalizationTerm * posteriorVariance;
}





/* getStdDevSynthesisRatePosteriorMean (RCPP EXPOSED)
 * Arguments: the number of samples from the end of the trace to examine, and the mixture element, both as unsigned values
 * Returns the standard deviation synthesis rate posterior mean of the mixture element.
 * This is calculated by simply gathering a number of traces from the end of the entire trace (up to all of it) to the end
 * of the standard deviation synthesis rate trace, and then getting the mean of these values.
*/

//' @name getStdDevSynthesisRatePosteriorMean
//' @title Calculate posterior mean of standard deviation parameter of lognormal describing distribution of synthesis rates
//' @param samples number of samples over which to calculate posterior mean
//' @param mixture mixture index to use. Should be number between 0 and n, where n is number of mixtures
//' @return returns posterior mean for standard deviation of lognormal distribution of synthesis rates
double Parameter::getStdDevSynthesisRatePosteriorMean(unsigned samples, unsigned mixture)
{
	double posteriorMean = 0.0;
	unsigned selectionCategory = getSelectionCategory(mixture);
	std::vector<double> stdDevSynthesisRateTrace = traces.getStdDevSynthesisRateTrace(selectionCategory);
	unsigned traceLength = lastIteration + 1;

	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getStdDevSynthesisRatePosteriorMean throws: Number of anticipated samples");
		my_printError("(%) is greater than the length of the available trace (%).", samples, traceLength);
		my_printError("Whole trace is used for posterior estimate!\n");

		samples = traceLength;
	}
	unsigned start = traceLength - samples;

	for (unsigned i = start; i < traceLength; i++)
		posteriorMean += stdDevSynthesisRateTrace[i];

	return posteriorMean / (double)samples;
}


//TODO: Considering renaming so that the two sides of the function match more than they currently do.
/* getSynthesisRatePosteriorMean (RCPP EXPOSED VIA WRAPPER)
 * Arguments: the number of samples from the end of the trace to examine, the gene to examine, and the mixture element.
 * Returns the posterior mean of the synthesis rate of the mixture element for a given gene.
 * This is calculated by simply gathering a number of traces from the end of the entire trace (up to all of it) to the end
 * of the gene's synthesis rate trace, and then getting the mean of these values.
 * Note: May return NaN if the gene was never in the category (expected resulted, OK).
 * Wrapped by getSynthesisRatePosteriorMeanByMixtureElementForGene on the R-side.
*/
double Parameter::getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, bool log_scale)
{
	float posteriorMean = 0.0;
	std::vector<float> synthesisRateTrace = traces.getSynthesisRateTraceForGene(geneIndex);
	if (synthesisRateTrace.size() == 1)
	{
		posteriorMean = synthesisRateTrace[0];
	}
	else
	{
		unsigned traceLength = lastIteration + 1;
		if (samples > lastIteration)
		{
			my_printError("Warning in Parameter::getSynthesisRatePosteriorMean throws: Number of anticipated samples");
			my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
						  samples, traceLength);

			samples = traceLength;
		}
		unsigned start = traceLength - samples;

		if(log_scale)
		{
			for (unsigned i = start; i < traceLength; i++)
			{
				synthesisRateTrace[i] = std::log10(synthesisRateTrace[i]);
			}
		}
		for (unsigned i = start; i < traceLength; i++)
		{
				posteriorMean += synthesisRateTrace[i];
		}
		posteriorMean = posteriorMean/float(samples);
	}
	return posteriorMean;
}





/* getCodonSpecificPosteriorMean (RCPP EXPOSED VIA WRAPPER)
 * Arguments: the mixture element, the number of samples from the end of the trace to examine, the codon to examine,
 * 			  the parameter type to examine, and whether or not it is with reference.
 * Returns the posterior mean of a codon-specific parameter of the mixture element.
 * This is calculated by simply gathering a number of traces from the end of the entire trace (up to all of it) to the end
 * of the codon-specific parameter's trace, and then getting the mean of these values.
 * Wrapped by getCodonSpecificPosteriorMeanForCodon on the R-side.
*/
double Parameter::getCodonSpecificPosteriorMean(unsigned element, unsigned samples, std::string &codon,
	unsigned paramType, bool withoutReference, bool byGene, bool log_scale)
{
	double posteriorMean = 0.0;
	std::vector<float> parameterTrace;
	if(byGene)
	{
		parameterTrace = traces.getCodonSpecificParameterTraceByGeneElementForCodon(
			element, codon, paramType, withoutReference);
	}else{
		parameterTrace = traces.getCodonSpecificParameterTraceByMixtureElementForCodon(
			element, codon, paramType, withoutReference);
	}

	unsigned traceLength = lastIteration + 1;


	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getCodonSpecificPosteriorMean throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
					  samples, traceLength);

		samples = traceLength;
	}
	unsigned start = traceLength - samples;

	for (unsigned i = start; i < traceLength; i++)
	{	
		if (log_scale)
		{
			posteriorMean += std::log10(parameterTrace[i]);
		}
		else
		{
			posteriorMean += parameterTrace[i];
		}
	}

	return posteriorMean / (double)samples;
}

//' @name getStdDevSynthesisRateVariance
//' @title Calculate variance of standard deviation parameter of lognormal describing distribution of synthesis rates
//' @param samples number of samples over which to calculate variance
//' @param mixture mixture index to use. Should be number between 0 and n, where n is number of mixtures
//' @param unbiased If TRUE, should calculate variance using unbiased (N-1). Otherwise, used biased (N) correction
//' @return returns variance for standard deviation of lognormal distribution of synthesis rates
double Parameter::getStdDevSynthesisRateVariance(unsigned samples, unsigned mixture, bool unbiased)
{
	unsigned selectionCategory = getSelectionCategory(mixture);
	std::vector<double> StdDevSynthesisRateTrace = traces.getStdDevSynthesisRateTrace(selectionCategory);
	unsigned traceLength = (unsigned)StdDevSynthesisRateTrace.size();
	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getSynthesisRateVariance throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
					  samples, traceLength);

		samples = traceLength;
	}
	double posteriorMean = getStdDevSynthesisRatePosteriorMean(samples, mixture);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		double difference = StdDevSynthesisRateTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1.0 / ((double)samples - 1.0)) : (1.0 / (double)samples);
	return normalizationTerm * posteriorVariance;
}


double Parameter::getSynthesisRateVariance(unsigned samples, unsigned geneIndex, bool unbiased, bool log_scale)
{
	double variance = 0.0;
	std::vector<float> synthesisRateTrace = traces.getSynthesisRateTraceForGene(geneIndex);
	if (synthesisRateTrace.size() != 1)
	{
		unsigned traceLength = lastIteration + 1;
		if (samples > traceLength)
		{
			my_printError("Warning in Parameter::getSynthesisRateVariance throws: Number of anticipated samples ");
			my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
						  samples, traceLength);

			samples = traceLength;
		}
		unsigned start = traceLength - samples;
		if(log_scale)
		{
			for (unsigned i = start; i < traceLength; i++)
			{
				synthesisRateTrace[i] = std::log10(synthesisRateTrace[i]);
			}
		}

		// NOTE: The loss of precision here is acceptable for storage purposes.
		float posteriorMean = (float)getSynthesisRatePosteriorMean(samples, geneIndex, log_scale);

		float posteriorVariance = 0.0;
		if (!std::isnan(posteriorMean))
		{
			double difference;
			for (unsigned i = start; i < traceLength; i++)
			{
				difference = synthesisRateTrace[i] - posteriorMean;
				posteriorVariance += difference * difference;
			}
		}
		float normalizationTerm = unbiased ? (1.0 / ((float)samples - 1.0)) : (1.0 / (float)samples);
		variance = normalizationTerm * posteriorVariance;
	}
	return variance;
}


double Parameter::getCodonSpecificVariance(unsigned mixtureElement, unsigned samples, std::string &codon,
	unsigned paramType, bool unbiased, bool withoutReference, bool log_scale)
{
	if (unbiased && samples == 1)
	{
		my_printError("Warning in Parameter::getCodonSpecificVariance throws: sample size is too small ");
		my_printError("to be considered unbiased (samples == 1). Setting as biased variance!\n");
		unbiased = false;
	}

	std::vector<float> parameterTrace = traces.getCodonSpecificParameterTraceByMixtureElementForCodon(
		mixtureElement, codon, paramType, withoutReference);
	unsigned traceLength = lastIteration + 1;
	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getCodonSpecificVariance throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
					  samples, traceLength);

		samples = traceLength;
	}

	double posteriorMean = getCodonSpecificPosteriorMean(mixtureElement, samples, codon, paramType, withoutReference, false,log_scale);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference;
	for (unsigned i = start; i < traceLength; i++)
	{
		if (log_scale)
		{
			difference = std::log10(parameterTrace[i]) - posteriorMean;
		}
		else
		{
			difference = parameterTrace[i] - posteriorMean;
		}
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1.0 / ((double)samples - 1.0)) : (1.0 / (double)samples);
	return normalizationTerm * posteriorVariance;
}



std::vector<double> Parameter::calculateQuantile(std::vector<float> &parameterTrace, unsigned samples, std::vector<double> probs, bool log_scale)
{
  unsigned traceLength = lastIteration + 1u;
    //unsigned traceEnd = parameterTrace.size() - (parameterTrace.size() - lastIteration); //currently unused
	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::calculateQuantile throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
					  samples, traceLength);

		samples = traceLength;
	}

    std::vector<double> samplesTrace(parameterTrace.begin() + (lastIteration - samples) + 1, (parameterTrace.begin() + lastIteration + 1));
    std::sort(samplesTrace.begin(), samplesTrace.end());

	if(log_scale)
	{
		for (unsigned i = 0u; i < samplesTrace.size(); i++)
		{
			samplesTrace[i] = std::log10(samplesTrace[i]);
		}
	}

    std::vector<double> retVec(probs.size());
    double N = samplesTrace.size();
    for (unsigned i = 0u; i < probs.size(); i++)
    {
		if( probs[i] < (2.0/3.0)/(N+(1.0/3.0)) )
		{
			retVec[i] = samplesTrace[0]; // first element
		}
		else if( probs[i] >= (N-(1.0/3.0))/(N+(1.0/3.0)) )
		{
			retVec[i] = samplesTrace[N - 1]; // last element
		}
		else
		{
			double h = (N*probs[i]) + (probs[i] + 1.0)/3.0;
			int low = std::floor(h);
			retVec[i] = samplesTrace[low] + (h - low)*(samplesTrace[low+1] - samplesTrace[low]);
		}
    }
    return retVec;
}

std::vector<double> Parameter::getExpressionQuantile(unsigned samples, unsigned geneIndex, std::vector<double> probs, bool log_scale)
{
	std::vector<double> quantile(probs.size());
	std::vector<float> parameterTrace = traces.getSynthesisRateTraceForGene(geneIndex);
	if (parameterTrace.size() == 1)
	{
		for (int i = 0; i < probs.size();i++)
		{
			quantile[i] = parameterTrace[0];
		}
	}
	else
	{
		quantile = calculateQuantile(parameterTrace, samples, probs, log_scale);
	}
	return quantile;
}

std::vector<double> Parameter::getCodonSpecificQuantile(unsigned mixtureElement, unsigned samples, std::string &codon,
	unsigned paramType, std::vector<double> probs, bool withoutReference, bool log_scale)
{
 	std::vector<float> parameterTrace = traces.getCodonSpecificParameterTraceByMixtureElementForCodon(
		mixtureElement, codon, paramType, withoutReference);

	return calculateQuantile(parameterTrace, samples, probs, log_scale);
}

//' @name getEstimatedMixtureAssignment
//' @title Get estimated mixture assignment for gene
//' @param samples number of samples over which to calculate mixture assignment
//' @param geneIndex corresponding index of gene in genome. Should be a number between 0 and length(genome) - 1. 
//' @return mixture returns value between 0 and n, where n is number of mixtures
unsigned Parameter::getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex)
{
	unsigned rv = 0u;
	double value = -1.0;
	std::vector <double> probabilities;
	probabilities = getEstimatedMixtureAssignmentProbabilities(samples, geneIndex);

	for (unsigned i = 0; i < probabilities.size(); i++)
	{
		if (value < probabilities[i])
		{
			value = probabilities[i];
			rv = i;
		}
	}
	return rv;
}


std::vector<double> Parameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
{
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	std::vector<double> probabilities(numMixtures, 0.0);
	unsigned traceLength = lastIteration + 1;

	if (samples > traceLength)
	{
		my_printError("Warning in Parameter::getEstimatedMixtureAssignmentProbabilities throws: Number of anticipated samples ");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n",
					  samples, traceLength);

		samples = traceLength;
	}

	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		unsigned value = mixtureAssignmentTrace[i];
		probabilities[value]++;
	}

	for (unsigned i = 0; i < numMixtures; i++)
		probabilities[i] /= (double)samples;

	return probabilities;
}


//--------------------------------------------------//
//---------- STATICS - Sorting Functions -----------//
//--------------------------------------------------//


/* sort array interval from first (included) to last (excluded)!!
// quick sort, sorting arrays a and b by a.
// Elements in b correspond to a, a will be sorted and it will be assured that b will be sorted by a */
void Parameter::quickSortPair(double a[], int b[], int first, int last)
{
	int pivotElement;

	if (first < last)
	{
		pivotElement = pivotPair(a, b, first, last);
		quickSortPair(a, b, first, pivotElement);
		quickSortPair(a, b, pivotElement + 1, last);
	}
}


int Parameter::pivotPair(double a[], int b[], int first, int last)
{
	int p = first;
	double pivotElement = a[first];

	for (int i = (first + 1) ; i < last ; i++)
	{
		/* If you want to sort the list in the other order, change "<=" to ">" */
		if (a[i] <= pivotElement)
		{
			p++;
			std::swap(a[i], a[p]);
			std::swap(b[i], b[p]);
		}
	}
	std::swap(a[p], a[first]);
	std::swap(b[p], b[first]);

	return p;
}


/* calculate SCUO values according to
// Wan et al. CodonO: a new informatics method for measuring synonymous codon usage bias within and across genomes
// International Journal of General Systems, Vol. 35, No. 1, February 2006, 109–125
// http://www.tandfonline.com/doi/pdf/10.1080/03081070500502967 */
double Parameter::calculateSCUO(Gene& gene)
{
	SequenceSummary *sequenceSummary = gene.getSequenceSummary();

	double totalDegenerateAACount = 0.0;
	unsigned maxAA = (SequenceSummary::AminoAcidArray).size();
	for (unsigned i = 0u; i < maxAA; i++)
	{
		std::string curAA = SequenceSummary::AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if (curAA == "X" || curAA == "M" || curAA == "W") continue;
		totalDegenerateAACount += (double)sequenceSummary->getAACountForAA(i);
	}

	double scuoValue = 0.0;
	for (unsigned i = 0u; i < maxAA; i++)
	{
		std::string curAA = SequenceSummary::AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if (curAA == "X" || curAA == "M" || curAA == "W") continue;
		double numDegenerateCodons = SequenceSummary::GetNumCodonsForAA(curAA);

		double aaCount = (double)sequenceSummary->getAACountForAA(i);
		if (aaCount == 0) continue;

		unsigned start, end;
		SequenceSummary::AAIndexToCodonRange(i, start, end, false);

		// calculate -sum(pij log(pij))
		double aaEntropy = 0.0;
		for (unsigned k = start; k < end; k++)
		{
			int currCodonCount = sequenceSummary->getCodonCountForCodon(k);
			if (currCodonCount == 0) continue;
			double codonProportion = (double)currCodonCount / aaCount;
			aaEntropy += codonProportion*std::log(codonProportion);
		}
		aaEntropy = -aaEntropy;
		// calculate max entropy -log(1/n_i)
		double maxEntropyForAA = -std::log(1.0 / numDegenerateCodons);
		// get normalized difference in entropy O_i
		double normalizedEntropyDiff = (maxEntropyForAA - aaEntropy) / maxEntropyForAA;

		// calculate the composition ratio F_i
		double compositionRatio = aaCount / totalDegenerateAACount;
		// SCUO is the sum(F_i * O_i) over all aa
		scuoValue += compositionRatio * normalizedEntropyDiff;
	}
	return scuoValue;
}


void Parameter::drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b),
									double* randomNumbers)
{
	for (unsigned i = 0u; i < draws; i++)
		randomNumbers[i] = (*proposal)(mean, sd);
}


void Parameter::drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumbers)
{
	for (unsigned i = 0u; i < draws; i++)
		randomNumbers[i] = (*proposal)(r);
}


double Parameter::randNorm(double mean, double sd)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rnorm(1, mean, sd);
	rv = xx[0];
#else
	std::normal_distribution<double> distribution(mean, sd);
	rv = distribution(generator);
#endif
	return rv;
}


double Parameter::randLogNorm(double m, double s)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rlnorm(1, m, s);
	rv = xx[0];
#else
	std::lognormal_distribution<double> distribution(m, s);
	rv = distribution(generator);
#endif
	return rv;
}


double Parameter::randExp(double r)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rexp(1, r);
	rv = xx[0];
#else
	std::exponential_distribution<double> distribution(r);
	rv = distribution(generator);
#endif
	return rv;
}


/* The R version and C++ differ because C++ uses the
// shape and scale parameter version while R uses the
// shape and rate. */
double Parameter::randGamma(double shape, double rate)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rgamma(1, shape, 1.0 / rate);
	rv = xx[0];
#else
	// Looking at the definition of std::gamma_distribution at
	// http://www.cplusplus.com/reference/random/gamma_distribution/
	// to verify this is correct.
	// It does appear correct and that std::gamma_distribution(shape, scale)
	// This is despite the documentation using alpha and beta to represent these
	// parameters.
	std::gamma_distribution<double> distribution(shape, 1.0 / rate);
	rv = distribution(generator);
#endif
	return rv;
}


// TODO: CHANGE THIS BACK TO DOUBLE*
void Parameter::randDirichlet(std::vector <double> &input, unsigned numElements, std::vector <double> &output)
{
	// draw y_i from Gamma(a_i, 1)
	// normalize y_i such that x_i = y_i / sum(y_i)

	double sumTotal = 0.0;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	for (unsigned i = 0; i < numElements; i++)
	{
		xx = rgamma(1, input[i], 1);
		output[i] = xx[0];
		sumTotal += xx[0];
	}
#else
	for (unsigned i = 0; i < numElements; i++)
	{
		std::gamma_distribution<double> distribution(input[i], 1);
		output[i] = distribution(generator);
		sumTotal += output[i];
	}
#endif
	for (unsigned i = 0; i < numElements; i++)
	{
		output[i] = output[i] / sumTotal;
	}
}


double Parameter::randUnif(double minVal, double maxVal)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, minVal, maxVal);
	rv = xx[0];
#else
	std::uniform_real_distribution<double> distribution(minVal, maxVal);
	rv = distribution(generator);
#endif
	return rv;
}


unsigned Parameter::randMultinom(std::vector <double> &probabilities, unsigned mixtureElements)
{
	// calculate cumulative sum to determine group boundaries
	double* cumulativeSum = new double[mixtureElements]();
	//std::vector<double> cumulativeSum(groups);
	cumulativeSum[0] = probabilities[0];

	for (unsigned i = 1u; i < mixtureElements; i++)
	{
		cumulativeSum[i] = cumulativeSum[i-1u] + probabilities[i];
	}
	// draw random number from U(0,1)
	double referenceValue;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, 0, 1);
	referenceValue = xx[0];
#else
	std::uniform_real_distribution<double> distribution(0, 1);
	referenceValue = distribution(generator);
#endif
	// check in which category the element falls
	unsigned returnValue = 0u;
	for (unsigned i = 0u; i < mixtureElements; i++)
	{
		if (referenceValue <= cumulativeSum[i])
		{
			returnValue = i;
			break;
		}
	}
	delete [] cumulativeSum;
	return returnValue;
}

unsigned Parameter::randMultinom(double *probabilities, unsigned mixtureElements)
{
	// calculate cumulative sum to determine group boundaries
	double* cumulativeSum = new double[mixtureElements]();
	//std::vector<double> cumulativeSum(groups);
	cumulativeSum[0] = probabilities[0];

	for (unsigned i = 1u; i < mixtureElements; i++)
	{
		cumulativeSum[i] = cumulativeSum[i - 1u] + probabilities[i];
	}
	// draw random number from U(0,1)
	double referenceValue;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, 0, 1);
	referenceValue = xx[0];
#else
	std::uniform_real_distribution<double> distribution(0, 1);
	referenceValue = distribution(generator);
#endif
	// check in which category the element falls
	unsigned returnValue = 0u;
	for (unsigned i = 0u; i < mixtureElements; i++)
	{
		if (referenceValue <= cumulativeSum[i])
		{
			returnValue = i;
			break;
		}
	}
	delete[] cumulativeSum;
	return returnValue;
}


double Parameter::densityNorm(double x, double mean, double sd, bool log)
{
	const double inv_sqrt_2pi = 0.3989422804014327;
	const double log_sqrt_2pi = 0.9189385332046727;
	double a = (x - mean) / sd;

	return log ? (-log_sqrt_2pi - std::log(sd) - (0.5 * a * a)) : ((inv_sqrt_2pi / sd) * std::exp(-0.5 * a * a));
}


double Parameter::densityLogNorm(double x, double mean, double sd, bool log)
{
	double returnValue = log ? -DBL_MAX : 0.0; // if log scale, instead of returning -Inf, -maximum possible double value
	// logN is only defined for x > 0 => all values less or equal to zero have probability 0
	if (x > 0.0)
	{
		const double inv_sqrt_2pi = 0.3989422804014327;
		const double log_sqrt_2pi = 0.9189385332046727;
		double a = (std::log(x) - mean) / sd;
		returnValue = log ? (-std::log(x * sd) - log_sqrt_2pi - (0.5 * a * a)) : ((inv_sqrt_2pi / (x * sd)) * std::exp(-0.5 * a * a));
	}
	return returnValue;
}





//-----------------------------------------------------------------------------------------------------//
//---------------------------------------- R SECTION --------------------------------------------------//
//-----------------------------------------------------------------------------------------------------//

#ifndef STANDALONE

//--------------------------------------------------------------//
//---------- Initialization and Restart Functions --------------//
//--------------------------------------------------------------//

//' @name initializeSynthesisRateByGenome
//' @title Initialize synthesis rates using SCUO values calcuated from the genome
//' @param genome a Genome object

void Parameter::initializeSynthesisRateByGenome(Genome& genome,double sd_phi)
{
	InitializeSynthesisRate(genome,sd_phi);
}

//' @name initializeSynthesisRateByRandom
//' @title Initialize synthesis rates by drawing a from a lognormal distribution with mean = -(sd_phi)^2/2 and sd = sd_phi
//' @param sd_phi a positive value which will be the standard deviation of the lognormal distribution

void Parameter::initializeSynthesisRateByRandom(double sd_phi)
{
	InitializeSynthesisRate(sd_phi);
}

//' @name initializeSynthesisRateByList
//' @title Initialize synthesis rates with values passed in as a list
//' @param expression a list of values to use as initial synthesis rate values. Should be same size as number of genes in genome.
void Parameter::initializeSynthesisRateByList(std::vector<double> expression)
{
	InitializeSynthesisRate(expression);
}


bool Parameter::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		my_printError("Error: Index % is out of bounds. Index must be between % & %\n", index, lowerbound, upperbound);
	}

	return check;
}





//----------------------------------------------------------------------//
//---------- Mixture Definition Matrix and Category Functions ----------//
//----------------------------------------------------------------------//


unsigned Parameter::getMutationCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delM + 1 : 0;
}


/* getSelectionCategoryForMixture
 * Arguments: A number representing a mixture element
 * Returns the selection category of the mixture element chosen.
 * See Note 1) above.
 * To implement the R version of this function, the index is also checked.
 * This is the R-wrapper for the C-side function "getSelectionCategory".
 */
unsigned Parameter::getSelectionCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delEta + 1 : 0;
}


/* getSynthesisRateCategoryForMixture
 * Arguments: A number representing a mixture element
 * Returns the synthesis rate category of the mixture element chosen.
 * See Note 1) above.
 * To implement the R version of this function, the index is also checked.
 * This is the R-wrapper for the C-side function "getSynthesisRateCategory".
 */
unsigned Parameter::getSynthesisRateCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delEta + 1 : 0;
}


std::vector<std::vector<unsigned>> Parameter::getCategories()
{
	unsigned size = (unsigned) categories.size();
	std::vector<std::vector<unsigned>> RV;
	for (unsigned i = 0; i < size; i++)
	{
		std::vector<unsigned> tmp;
		tmp.push_back(categories[i].delM);
		tmp.push_back(categories[i].delEta);
		RV.push_back(tmp);
	}

	return RV;
}


void Parameter::setCategories(std::vector<std::vector<unsigned>> _categories)
{
	for (unsigned i = 0; i < _categories.size(); i++)
	{
		categories.push_back(mixtureDefinition());
		categories[i].delM = _categories[i][0];
		categories[i].delEta = _categories[i][1];
	}
}


void Parameter::setCategoriesForTrace()
{
	traces.setCategories(categories);
}


void Parameter::setNumMutationCategories(unsigned _numMutationCategories)
{
	numMutationCategories = _numMutationCategories;
}


void Parameter::setNumSelectionCategories(unsigned _numSelectionCategories)
{
	numSelectionCategories = _numSelectionCategories;
}





//-----------------------------------------------//
//---------- Synthesis Rate Functions -----------//
//-----------------------------------------------//

//' @name getSynthesisRate
//' @title Get current synthesis rates for all genes and all mixtures 
//' @return 2 by 2 vector of numeric values
std::vector<std::vector<double>> Parameter::getSynthesisRateR()
{
	return currentSynthesisRateLevel;
}


std::vector<double> Parameter::getCurrentSynthesisRateForMixture(unsigned mixture)
{
	bool checkMixture = checkIndex(mixture, 1, numMixtures);
	unsigned exprCat = 0u;
	if (checkMixture)
	{
		exprCat = getSynthesisRateCategory(mixture - 1);
	}
	else
	{
		my_printError("WARNING: Mixture element % NOT found. Mixture element 1 is returned instead.\n", mixture);
	}
	return currentSynthesisRateLevel[exprCat];
}


//------------------------------------------------------------------//
//---------- Posterior, Variance, and Estimates Functions ----------//
//------------------------------------------------------------------//


/* getCodonSpecificPosteriorMeanForCodon
 * Arguments: the mixture element, the number of samples from the end of the trace to examine, the codon to examine,
 * 			  the parameter type to examine, and whether or not it is with reference.
 * Returns the posterior mean of a codon-specific parameter of the mixture element.
 * This is calculated by simply gathering a number of traces from the end of the entire trace (up to all of it) to the end
 * of the codon-specific parameter's trace, and then getting the mean of these values.
 * To implement the R version of this function, the index is also checked.
 * This is the R-wrapper for the C-side function "getCodonSpecificPosteriorMean".
*/


 //' @name getCodonSpecificPosteriorMeanForCodon
 //' @title Calculate codon-specific parameter (CSP) posterior mean
 //' @param mixtureElement mixture to calculate CSP posterior mean. Should be between 1 and n, where n is number of mixtures.
 //' @param samples number of samples to use for calculating posterior mean
 //' @param codon codon to calculate CSP
 //' @param paramType CSP to calculate posterior mean for. 0: Mutation (ROC,FONSE) or Alpha (PA, PANSE). 1: Selection (ROC,FONSE), Lambda (PANSE), Lambda^prime (PA). 2: NSERate (PANSE) 
 //' @param withoutReference If model uses reference codon, then ignore this codon (fixed at 0). Should be TRUE for ROC and FONSE. Should be FALSE for PA and PANSE.
 //' @param log_scale If true, calculate posterior mean on log scale. Should only be used for PA and PANSE.
double Parameter::getCodonSpecificPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon,
	unsigned paramType, bool withoutReference,bool log_scale)
{
	double rv = -1.0;
	codon[0] = (char)std::toupper(codon[0]);
	codon[1] = (char)std::toupper(codon[1]);
	codon[2] = (char)std::toupper(codon[2]);
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		rv = getCodonSpecificPosteriorMean(mixtureElement - 1, samples, codon, paramType, withoutReference, false, log_scale);
	}
	return rv;
}

 //' @name getCodonSpecificPosteriorVarianceForCodon
 //' @title Calculate codon-specific parameter (CSP) variance
 //' @param mixtureElement mixture to calculate CSP variance. Should be between 1 and n, where n is number of mixtures.
 //' @param samples number of samples to use for calculating variance
 //' @param codon codon to calculate CSP
 //' @param paramType CSP to calculate variance for. 0: Mutation (ROC,FONSE) or Alpha (PA, PANSE). 1: Selection (ROC,FONSE), Lambda (PANSE), Lambda^prime (PA). 2: NSERate (PANSE) 
 //' @param unbiased If TRUE, should calculate variance using unbiased (N-1). Otherwise, used biased (N) correction
 //' @param withoutReference If model uses reference codon, then ignore this codon (fixed at 0). Should be TRUE for ROC and FONSE. Should be FALSE for PA and PANSE.
 //' @param log_scale If true, calculate posterior mean on log scale. Should only be used for PA and PANSE.
double Parameter::getCodonSpecificVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon,
	unsigned paramType, bool unbiased, bool withoutReference, bool log_scale)
{
	double rv = -1.0;
	codon[0] = (char)std::toupper(codon[0]);
	codon[1] = (char)std::toupper(codon[1]);
	codon[2] = (char)std::toupper(codon[2]);
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		rv = getCodonSpecificVariance(mixtureElement - 1, samples, codon, paramType, unbiased, withoutReference, log_scale);
	}
	return rv;
}

 //' @name getCodonSpecificQuantilesForCodon
 //' @title Calculate quantiles of CSP traces
 //' @param mixtureElement mixture to calculate CSP variance. Should be between 1 and n, where n is number of mixtures.
 //' @param samples number of samples to use for calculating variance
 //' @param codon codon to calculate CSP
 //' @param paramType CSP to calculate variance for. 0: Mutation (ROC,FONSE) or Alpha (PA, PANSE). 1: Selection (ROC,FONSE), Lambda (PANSE), Lambda^prime (PA). 2: NSERate (PANSE) 
 //' @param probs vector of two doubles between 0 and 1 indicating range over which to calculate quantiles. <0.0275, 0.975> would give 95% quantiles.
 //' @param withoutReference If model uses reference codon, then ignore this codon (fixed at 0). Should be TRUE for ROC and FONSE. Should be FALSE for PA and PANSE.
 //' @param log_scale If true, calculate posterior mean on log scale. Should only be used for PA and PANSE.
//'  @return vector representing lower and upper bound of quantile
std::vector<double> Parameter::getCodonSpecificQuantileForCodon(unsigned mixtureElement, unsigned samples,
	std::string &codon, unsigned paramType, std::vector<double> probs, bool withoutReference, bool log_scale)
{
	std::vector<double> rv;
	codon[0] = (char)std::toupper(codon[0]);
	codon[1] = (char)std::toupper(codon[1]);
	codon[2] = (char)std::toupper(codon[2]);
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
        rv = getCodonSpecificQuantile(mixtureElement - 1, samples, codon, paramType, probs, withoutReference, log_scale);
    }
    return rv;
}

std::vector<double> Parameter::getExpressionQuantileForGene(unsigned samples,
	unsigned geneIndex, std::vector<double> probs, bool log_scale)
{
	std::vector<double> rv;
	bool checkGene = checkIndex(geneIndex, 1, (unsigned) mixtureAssignment.size());
	if (checkGene)
	{
        rv = getExpressionQuantile(samples, geneIndex - 1, probs, log_scale);
    }
    return rv;
}

/* getSynthesisRatePosteriorMeanForGene
 * Arguments: the number of samples from the end of the trace to examine, the gene to examine, and the mixture element.
 * Returns the posterior mean of the synthesis rate of the mixture element for a given gene.
 * This is calculated by simply gathering a number of traces from the end of the entire trace (up to all of it) to the end
 * of the gene's synthesis rate trace, and then getting the mean of these values.
 * To implement the R version of this function, the index is also checked.
 * This is the R-wrapper for the C-side function "getSynthesisRatePosteriorMean".
*/

//' @name getSynthesisRatePosteriorMeanForGene
//' @title Get posterior mean synthesis rate value for a gene
//' @param samples number of samples over which to calculate mean
//' @param geneIndex corresponding index of gene in genome for which posterior mean synthesis rate will be calculated. Should be a number between 1 and length(genome) 
//' @param log_scale Calculate posterior mean on log scale
//' @return posterior mean synthesis rate for gene
double Parameter::getSynthesisRatePosteriorMeanForGene(unsigned samples, unsigned geneIndex, bool log_scale)
{
	double rv = -1.0;
	bool checkGene = checkIndex(geneIndex, 1, (unsigned) mixtureAssignment.size());
	if (checkGene)
	{
		rv = getSynthesisRatePosteriorMean(samples, geneIndex - 1, log_scale);
	}
	return rv;
}

//' @name getSynthesisRatePosteriorVarianceForGene
//' @title Get synthesis rate variance for a gene
//' @param samples number of samples over which to calculate variance
//' @param geneIndex corresponding index of gene in genome for which synthesis rate variance will be calculated. Should be a number between 1 and length(genome) 
//' @param unbiased Should calculate variance using unbiased (N-1) or biased (N) correction
//' @param log_scale Calculate variance on log scale
//' @return posterior mean synthesis rate for gene
double Parameter::getSynthesisRateVarianceForGene(unsigned samples, unsigned geneIndex, bool unbiased, bool log_scale)
{
	double rv = -1.0;
	bool checkGene = checkIndex(geneIndex, 1, (unsigned) mixtureAssignment.size());
	if (checkGene)
	{
		rv = getSynthesisRateVariance(samples, geneIndex - 1, unbiased, log_scale);
	}
	return rv;
}


unsigned Parameter::getEstimatedMixtureAssignmentForGene(unsigned samples, unsigned geneIndex)
{
	bool check = checkIndex(geneIndex, 1, (unsigned) mixtureAssignment.size());
	return check ? getEstimatedMixtureAssignment(samples, geneIndex - 1) + 1 : 0;
}


std::vector<double> Parameter::getEstimatedMixtureAssignmentProbabilitiesForGene(unsigned samples, unsigned geneIndex)
{
	std::vector <double> probabilities;
	bool check = checkIndex(geneIndex, 1, (unsigned) mixtureAssignment.size());
	if (check)
	{
		probabilities = getEstimatedMixtureAssignmentProbabilities(samples, geneIndex - 1);
	}
	return probabilities;
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


SEXP Parameter::calculateSelectionCoefficientsR(unsigned sample)
{
	NumericMatrix RSelectionCoefficents(mixtureAssignment.size(), 61); //61 due to stop codons
	std::vector<std::vector<double>> selectionCoefficients = calculateSelectionCoefficients(sample);
	unsigned index = 0;
	for (unsigned i = 0; i < selectionCoefficients.size(); i++)
	{
		for (unsigned j = 0; j < selectionCoefficients[i].size(); j++, index++)
		{
			RSelectionCoefficents(i,j) = selectionCoefficients[i][j];
		}
	}
	return RSelectionCoefficents;
}


std::vector<unsigned> Parameter::getMixtureAssignmentR()
{
	return mixtureAssignment;
}


void Parameter::setMixtureAssignmentR(std::vector<unsigned> _mixtureAssignment)
{
	mixtureAssignment = _mixtureAssignment;
}


unsigned Parameter::getMixtureAssignmentForGeneR(unsigned geneIndex)
{
	unsigned rv = 0;
	bool check = checkIndex(geneIndex, 1, (unsigned)mixtureAssignment.size());
	if (check)
	{
		rv = getMixtureAssignment(geneIndex - 1) + 1;
	}
	return rv;
}


void Parameter::setMixtureAssignmentForGene(unsigned geneIndex, unsigned value)
{
	bool check = checkIndex(geneIndex, 1, (unsigned) mixtureAssignment.size());
	if (check)
	{
		mixtureAssignment[geneIndex - 1] = value;
	}
}


void Parameter::setNumMixtureElements(unsigned _numMixtures)
{
    numMixtures = _numMixtures;
}


#endif
