#include "include/ROC/ROCParameter.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


ROCParameter::ROCParameter() : Parameter()
{
	//CTOR
	bias_csp = 0;
	//mutation_prior_sd = 0.35;
	mutation_prior_mean.resize(40);
	mutation_prior_sd.resize(40);
	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);
}


ROCParameter::ROCParameter(std::string filename) : Parameter(22)
{
	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);
	initFromRestartFile(filename);
}


ROCParameter::ROCParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
		std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState) :
		Parameter(22)
{
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initROCParameterSet();
}


ROCParameter& ROCParameter::operator=(const ROCParameter& rhs)
{
	if (this == &rhs)
		return *this; // handle self assignment

	Parameter::operator=(rhs);

	covarianceMatrix = rhs.covarianceMatrix;
	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	noiseOffset = rhs.noiseOffset;
	noiseOffset_proposed = rhs.noiseOffset_proposed;
	std_NoiseOffset = rhs.std_NoiseOffset;
	numAcceptForNoiseOffset = rhs.numAcceptForNoiseOffset;

	return *this;
}


ROCParameter::~ROCParameter()
{
	//DTOR
}





//---------------------------------------------------------------//
//----------- Initialization, Restart, Index Checking -----------//
//---------------------------------------------------------------//


void ROCParameter::initROCParameterSet()
{
	//mutation_prior_sd = 0.35;
	mutation_prior_mean.resize(numMutationCategories);
	mutation_prior_sd.resize(numMutationCategories);
	for (int i=0; i < numMutationCategories; i++)
	{
	  // POTENTIAL ISSUE: Shouldn't we use (numParam) instead of (40)
		mutation_prior_mean[i].resize(40);
		mutation_prior_sd[i].resize(40);
		std::vector<double> tmp(40, 0.0);
		mutation_prior_mean[i] = tmp;
		mutation_prior_sd[i] = tmp;
	}
	// POTENTIAL ISSUE: This seems to assume splitSer = True
	groupList = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z"};
	// proposal bias and std for codon specific parameter
	bias_csp = 0;


	// for (unsigned i = 0; i < getNumObservedPhiSets(); i++)
	// {
	// 	noiseOffset[i] = 0.1;
	// 	noiseOffset_proposed[i] = 0.1;
	// 	std_NoiseOffset[i] = 0.1;
	// 	observedSynthesisNoise[i] = 0.1;
	// 	numAcceptForNoiseOffset[i] = 0;
	// }

	//may need getter fcts

	proposedCodonSpecificParameter.resize(2);
	currentCodonSpecificParameter.resize(2);

	currentCodonSpecificParameter[dM].resize(numMutationCategories);
	proposedCodonSpecificParameter[dM].resize(numMutationCategories);

	for (unsigned i = 0u; i < numMutationCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		currentCodonSpecificParameter[dM][i] = tmp;
		proposedCodonSpecificParameter[dM][i] = tmp;
	}

	currentCodonSpecificParameter[dEta].resize(numSelectionCategories);
	proposedCodonSpecificParameter[dEta].resize(numSelectionCategories);
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		proposedCodonSpecificParameter[dEta][i] = tmp;
		currentCodonSpecificParameter[dEta][i] = tmp;
	}

  for (unsigned i = 0; i < maxGrouping; i++)
  {
    std::string aa = SequenceSummary::AminoAcidArray[i];
    unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa, true);
    CovarianceMatrix m((numMutationCategories + numSelectionCategories) * numCodons);
    m.choleskyDecomposition();
    covarianceMatrix.push_back(m);
  }
}


void ROCParameter::initROCValuesFromFile(std::string filename)
{
	bool old_format = true;
	std::ifstream input;
	covarianceMatrix.resize(maxGrouping);
	std::vector <double> mat;
	input.open(filename.c_str());
	if (input.fail())
	{
		my_printError("Error opening file % to initialize from restart file.\n", filename.c_str());
		my_printError("please use absolute path");
		return;
	}
	else
	  {
		std::string tmp, variableName;
		unsigned cat = 0;
		while (getline(input, tmp))
		{
			int flag;
			if (tmp[0] == '>')
				flag = 1;
			else if (input.eof() || tmp == "\n")
				flag = 2;
			else if (tmp[0] == '#')
				flag = 3;
			else
				flag = 4;

			if (flag == 1)
			{
				mat.clear();
				cat = 0;
				variableName = tmp.substr(1,tmp.size()-2);
				if (variableName == "covarianceMatrix")
				{
					getline(input,tmp);
					//char aa = tmp[0];
					cat = SequenceSummary::AAToAAIndex(tmp); // ????
				}
			}
			else if (flag == 2)
				my_print("here\n");
			else if (flag == 3) //user comment, continue
				continue;
			else
			{
				std::istringstream iss;
				if (variableName == "currentMutationParameter")
				{
					if (tmp == "***")
					{
						currentCodonSpecificParameter[dM].resize(currentCodonSpecificParameter[dM].size() + 1);
						cat++;
					}
					else if (tmp == "\n")
						continue;
					else
					{
						double val;
						iss.str(tmp);
						while (iss >> val)
						{
							currentCodonSpecificParameter[dM][cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "currentSelectionParameter")
				{
					if (tmp == "***")
					{
						currentCodonSpecificParameter[dEta].resize(currentCodonSpecificParameter[dEta].size() + 1);
						cat++;
					}
					else if (tmp == "\n")
						continue;
					else
					{
						double val;
						iss.str(tmp);
						while (iss >> val)
						{
							currentCodonSpecificParameter[dEta][cat - 1].push_back(val);
						}
					}
				}

				else if (variableName == "std_csp")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						std_csp.push_back(val);
					}
				}
				// else if (variableName == "noiseOffset")
				// {
				// 	double val;
				// 	iss.str(tmp);
				// 	while (iss >> val)
				// 	{
				// 		noiseOffset.push_back(val);
				// 		noiseOffset_proposed.push_back(val);
				// 	}
				// }
				// else if (variableName == "observedSynthesisNoise")
				// {
				// 	double val;
				// 	iss.str(tmp);
				// 	while (iss >> val)
				// 	{
				// 		observedSynthesisNoise.push_back(val);
				// 	}
				// }
				// else if (variableName == "std_NoiseOffset")
				// {
				// 	double val;
				// 	iss.str(tmp);
				// 	while (iss >> val)
				// 	{
				// 		std_NoiseOffset.push_back(val);
				// 	}
				// }
				else if (variableName == "covarianceMatrix")
				{
					if (tmp == "***") //end of matrix
					{
						CovarianceMatrix CM(mat);
						CM.choleskyDecomposition();
						covarianceMatrix[cat] = CM;
					}
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						mat.push_back(val);
					}
				}
				else if (variableName == "mutation_prior_mean")
				{
					if (tmp == "***")
					{
						old_format = false;
						mutation_prior_mean.resize(mutation_prior_mean.size() + 1);
						cat++;
					}
					else if (tmp == "\n")
						continue;
					else
					{
						double val;
						iss.str(tmp);
						while (iss >> val)
						{
							mutation_prior_mean[cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "mutation_prior_sd")
				{
					if (tmp == "***")
					{
						old_format = false;
						mutation_prior_sd.resize(mutation_prior_sd.size() + 1);
						cat++;
					}
					else if (tmp == "\n")
						continue;
					else
					{
						double val = 0.0;
						iss.str(tmp);
						if (!old_format)
						{
							while (iss >> val)
							{
								mutation_prior_sd[cat - 1].push_back(val);
							}
						}
						else
						{
							mutation_prior_sd.resize(numMutationCategories);
							for (int i = 0; i < numMutationCategories;i++)
							{
								mutation_prior_sd[i].resize(40);
								for (int j = 0; j < 40; j++)
								{
									mutation_prior_sd[i][j] = val;
								}
							}
						}
					}
				}
			}
		}
	}
	input.close();
	if (!old_format)
	{
		mutation_prior_mean.resize(numMutationCategories);
		for (int i = 0; i < numMutationCategories;i++)
		{
			mutation_prior_mean[i].resize(40);
			for (int j = 0; j < 40; j++)
			{
				mutation_prior_mean[i][j] = 0.0;
			}
		}
	}
	//init other values
	//numAcceptForNoiseOffset.resize(obsPhiSets, 0);
	bias_csp = 0;
	proposedCodonSpecificParameter[dM].resize(numMutationCategories);
	proposedCodonSpecificParameter[dEta].resize(numSelectionCategories);
	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		proposedCodonSpecificParameter[dM][i] = currentCodonSpecificParameter[dM][i];
	}
	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		proposedCodonSpecificParameter[dEta][i] = currentCodonSpecificParameter[dEta][i];
	}
}


void ROCParameter::writeEntireRestartFile(std::string filename)
{
	writeBasicRestartFile(filename);
	writeROCRestartFile(filename);
}


void ROCParameter::writeROCRestartFile(std::string filename)
{
	std::ofstream out;
	out.open(filename.c_str(), std::ofstream::app);
	if (out.fail())
		my_printError("Error opening file % to write restart file.\n", filename.c_str());
		
	else
	{
		std::ostringstream oss;
		unsigned j;
		// oss << ">noiseOffset:\n";
		// for (j = 0; j < noiseOffset.size(); j++)
		// {
		// 	oss << noiseOffset[j];
		// 	if ((j + 1) % 10 == 0)
		// 		oss << "\n";
		// 	else
		// 		oss << " ";
		// }
		// if (j % 10 != 0)
		// 	oss << "\n";

		// oss << ">observedSynthesisNoise:\n";
		// for (j = 0; j < observedSynthesisNoise.size(); j++)
		// {
		// 	oss << observedSynthesisNoise[j];
		// 	if ((j + 1) % 10 == 0)
		// 		oss << "\n";
		// 	else
		// 		oss << " ";
		// }
		// if (j % 10 != 0)
		// 	oss << "\n";

		//oss << ">mutation_prior_sd:\n" << mutation_prior_sd << "\n";
		oss << ">mutation_prior_mean:\n";
		for (unsigned i = 0; i < mutation_prior_mean.size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < mutation_prior_mean[i].size(); j++)
			{
				oss << mutation_prior_mean[i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}
		oss << ">mutation_prior_sd:\n";
		for (unsigned i = 0; i < mutation_prior_sd.size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < mutation_prior_sd[i].size(); j++)
			{
				oss << mutation_prior_sd[i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}
		// oss << ">std_NoiseOffset:\n";
		// for (j = 0; j < std_NoiseOffset.size(); j++)
		// {
		// 	oss << std_NoiseOffset[j];
		// 	if ((j + 1) % 10 == 0)
		// 		oss << "\n";
		// 	else
		// 		oss << " ";
		// }
		// if (j % 10 != 0)
		// 	oss << "\n";
		oss << ">std_csp:\n";
		for (j = 0; j < std_csp.size(); j++)
		{
			oss << std_csp[j];
			if ((j + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (j % 10 != 0)
			oss << "\n";
		oss << ">currentMutationParameter:\n";
		for (unsigned i = 0; i < currentCodonSpecificParameter[dM].size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentCodonSpecificParameter[dM][i].size(); j++)
			{
				oss << currentCodonSpecificParameter[dM][i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}

		oss << ">currentSelectionParameter:\n";
		for (unsigned i = 0; i < currentCodonSpecificParameter[dEta].size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentCodonSpecificParameter[dEta][i].size(); j++)
			{
				oss << currentCodonSpecificParameter[dEta][i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}

		for (unsigned i = 0; i < groupList.size(); i++)
		{
			std::string aa = groupList[i];
			oss << ">covarianceMatrix:\n" << aa << "\n";
			CovarianceMatrix m = covarianceMatrix[SequenceSummary::AAToAAIndex(aa)];
			std::vector<double>* tmp = m.getCovMatrix();
			int size = m.getNumVariates();
			for (unsigned k = 0; k < size * size; k++)
			{
				if (k % size == 0 && k != 0) { oss << "\n"; }
				oss << tmp->at(k) << "\t";
			}
			oss << "\n***\n";
		}
		std::string output = oss.str();
		out << output;
	}
	out.close();
}


void ROCParameter::initFromRestartFile(std::string filename)
{
	initBaseValuesFromFile(filename);
	initROCValuesFromFile(filename);
}


void ROCParameter::initAllTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
	traces.initializeROCTrace(samples, num_genes, numMutationCategories, numSelectionCategories, numParam,
						 numMixtures, categories, maxGrouping, obsPhiSets,currentSynthesisRateLevel[0],mixtureAssignment,estimateSynthesisRate);
}


//' @name initMutationCategories
//' @title Initialize values for mutation CSP. File should be of comma-separated with header. Three columns should be of order Amino_acid,Codon,Value
//' @param files list of files containing starting values. Number of files should equal the number of categories.
//' @param numCategories number of mutation categories (should be less than or equal to number of mixtures)
//' @param fix Can use this parameter to fix mutation at current values (won't change over course of MCMC run)

void ROCParameter::initMutationCategories(std::vector<std::string> files, unsigned numCategories, bool fix)
{
	for (unsigned category = 0; category < numCategories; category++)
	{
		//Open the file for the category
		std::ifstream currentFile;
		currentFile.open(files[category].c_str());
		if (currentFile.fail())
		{
			my_printError("Error opening file % to initialize mutation values.\n", category);
			my_printError("please use absolute path");
			return;
		}

		else
		{
			std::string tmp;
			currentFile >> tmp; //The first line is a header (Amino Acid, Codon, Value, Std_deviation)
			fix_dM = fix;
			while (currentFile >> tmp)
			{
				//Get the Codon and Index
				std::size_t pos = tmp.find(',', 2); //Amino Acid and a comma will always be the first 2 characters
				std::string codon = tmp.substr(2, pos - 2);
				unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);

				//get the value to store
				std::size_t pos2 = tmp.find(',', pos + 1);
				//my_print("%\n", tmp.substr(pos + 1, pos2 - pos - 1 ));
				double value = std::atof(tmp.substr(pos + 1, pos2 - pos - 1).c_str());

				currentCodonSpecificParameter[dM][category][codonIndex] = value;
				proposedCodonSpecificParameter[dM][category][codonIndex] = value;
			}
		}

		currentFile.close();
	} //END OF A CATEGORY/FILE
}

//' @name initSelectionCategories
//' @title Initialize values for selection CSP. File should be of comma-separated with header. Three columns should be of order Amino_acid,Codon,Value
//' @param files list of files containing starting values. Number of files should equal the number of categories.
//' @param numCategories number of mutation categories (should be less than or equal to number of mixtures)
//' @param fix Can use this parameter to fix selection at current values (won't change over course of MCMC run)

void ROCParameter::initSelectionCategories(std::vector<std::string> files, unsigned numCategories, bool fix)
{
	for (unsigned category = 0; category < numCategories; category++)
	{
		//Open the file for the category
		std::ifstream currentFile;
		currentFile.open(files[category].c_str());
		if (currentFile.fail())
		{
			my_printError("Error opening file % to initialize selection values.\n", category);
			my_printError("please use absolute path");
			return;
		}

		else
		{
			std::string tmp;
			currentFile >> tmp; //The first line is a header (Amino Acid, Codon, Value, Std_deviation)
			fix_dEta = fix;
			while (currentFile >> tmp)
			{
				//Get the Codon and Index
				std::size_t pos = tmp.find(',', 2); //Amino Acid and a comma will always be the first 2 characters
				std::string codon = tmp.substr(2, pos - 2);
				unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);

				//get the value to store
				std::size_t pos2 = tmp.find(',', pos + 1);
				//	my_print("%\n", tmp.substr(pos + 1, pos2 - pos - 1 ));
				double value = std::atof(tmp.substr(pos + 1, pos2 - pos - 1).c_str());

				currentCodonSpecificParameter[dEta][category][codonIndex] = value;
				proposedCodonSpecificParameter[dEta][category][codonIndex] = value;
			}
		}

		currentFile.close();
	} //END OF A CATEGORY/FILE
}





//--------------------------------------//
//---------- Trace Functions -----------//
//--------------------------------------//

// void ROCParameter::updateObservedSynthesisNoiseTraces(unsigned sample)
// {
// 	for (unsigned i = 0; i < observedSynthesisNoise.size(); i++)
// 	{
// 		traces.updateObservedSynthesisNoiseTrace(i, sample, observedSynthesisNoise[i]);
// 	}
// }


// void ROCParameter::updateNoiseOffsetTraces(unsigned sample)
// {
// 	for (unsigned i = 0; i < noiseOffset.size(); i++)
// 	{
// 		traces.updateSynthesisOffsetTrace(i, sample, noiseOffset[i]);
// 	}
// }

void ROCParameter::updateCodonSpecificParameterTrace(unsigned sample, std::string grouping)
{
	traces.updateCodonSpecificParameterTraceForAA(sample, grouping, currentCodonSpecificParameter[dM], dM);
	traces.updateCodonSpecificParameterTraceForAA(sample, grouping, currentCodonSpecificParameter[dEta], dEta);
}


//' @name fixDM
//' @title Fix the value of mutation its current value

void ROCParameter::fixDM()
{
	fix_dM = true;
}

//' @name fixDEta
//' @title Fix the value of selection its current value

void ROCParameter::fixDEta()
{
	fix_dEta = true;
}

bool ROCParameter::isDMFixed()
{
	return fix_dM;
}

bool ROCParameter::isDEtaFixed()
{
	return fix_dEta;
}

//------------------------------------------//
//---------- Covariance Functions ----------//
//------------------------------------------//


CovarianceMatrix& ROCParameter::getCovarianceMatrixForAA(std::string aa)
{
	aa[0] = (char) std::toupper(aa[0]);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa) -> second;
	return covarianceMatrix[aaIndex];
}





//------------------------------------------------------//
//---------- observedSynthesisNoise Functions ----------//
//------------------------------------------------------//


// double ROCParameter::getObservedSynthesisNoise(unsigned index)
// {
// 	return observedSynthesisNoise[index];
// }


// void ROCParameter::setObservedSynthesisNoise(unsigned index, double se)
// {
// 	observedSynthesisNoise[index] = se;
// }





// //-------------------------------------------//
// //---------- noiseOffset Functions ----------//
// //-------------------------------------------//


// double ROCParameter::getNoiseOffset(unsigned index, bool proposed)
// {
// 	return (proposed ? noiseOffset_proposed[index] : noiseOffset[index]);
// }


// double ROCParameter::getCurrentNoiseOffsetProposalWidth(unsigned index)
// {
// 	return std_NoiseOffset[index];
// }


// void ROCParameter::proposeNoiseOffset()
// {
// 	for (unsigned i = 0; i < getNumObservedPhiSets(); i++)
// 	{
// 		noiseOffset_proposed[i] = randNorm(noiseOffset[i], std_NoiseOffset[i]);
// 	}
// }


// void ROCParameter::setNoiseOffset(unsigned index, double _noiseOffset)
// {
// 	noiseOffset[index] = _noiseOffset;
// }


// void ROCParameter::updateNoiseOffset(unsigned index)
// {
// 	noiseOffset[index] = noiseOffset_proposed[index];
// 	numAcceptForNoiseOffset[index]++;
// }


//-----------------------------------//
//---------- Noise Functions --------//
//-----------------------------------//

// void ROCParameter::setInitialValuesForSepsilon(std::vector<double> seps)
// {
// 	if (seps.size() == observedSynthesisNoise.size())
// 	{
// 		for (unsigned i = 0; i < observedSynthesisNoise.size(); i++)
// 		{
// 			observedSynthesisNoise[i] = seps[i];
// 		}
// 	}
// 	else
// 	{
// 		my_printError("ROCParameter::setInitialValuesForSepsilon number of initial values (%) does not match number of expression sets (%)",
// 					  seps.size(), observedSynthesisNoise.size());
// 	}
// }


//-----------------------------------//
//---------- CSP Functions ----------//
//-----------------------------------//


double ROCParameter::getCurrentCodonSpecificProposalWidth(unsigned aa)
{
	unsigned aaStart, aaEnd;
	SequenceSummary::AAIndexToCodonRange(aa, aaStart, aaEnd, true);
	return std_csp[aaStart];
}

/* Cedric: I decided to use a normal distribution to propose Sphi and phi instead of a lognormal because:
 * 1. It is a symmetric distribution and you therefore do not have to account for the unsymmetry in jump probabilities
 * 2. The one log and exp operation that have to be performed per parameter are cheaper than the operations necessary to draw from a lognormal
 * 3. phi has to be on a non log scale for the likelihood evaluation thus it does not help to keep phi on th elog scale all the time
 * 4. the adjustment of the likelihood by the jacobian that arises from this transformation is cheap and by grouping everything in one class it takes place more or less at the same place
*/
void ROCParameter::proposeCodonSpecificParameter()
{

	for (unsigned k = 0; k < getGroupListSize(); k++)
	{
		std::vector<double> iidProposed;
		std::string aa = getGrouping(k);
		unsigned aaStart, aaEnd;
		SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
		unsigned numCodons = aaEnd - aaStart;
		for (unsigned i = 0u; i < (numCodons * (numMutationCategories + numSelectionCategories)); i++)
		{
			iidProposed.push_back(randNorm(0.0, 1.0));
		}

		std::vector<double> covaryingNums;
		covaryingNums = covarianceMatrix[SequenceSummary::AAToAAIndex(aa)].transformIidNumbersIntoCovaryingNumbers(iidProposed);
		for (unsigned i = 0; i < numMutationCategories; i++)
		{
			for (unsigned j = i * numCodons, l = aaStart; j < (i * numCodons) + numCodons; j++, l++)
			{
				if (fix_dM)
				{
					proposedCodonSpecificParameter[dM][i][l] = currentCodonSpecificParameter[dM][i][l];
				}
				else if (propose_by_prior)
				{
					proposedCodonSpecificParameter[dM][i][l] = randNorm(mutation_prior_mean[i][l],mutation_prior_sd[i][l]);
				}
				else
				{
					proposedCodonSpecificParameter[dM][i][l] = currentCodonSpecificParameter[dM][i][l] + covaryingNums[j];
				}
			}
		}
		for (unsigned i = 0; i < numSelectionCategories; i++)
		{
			for (unsigned j = i * numCodons, l = aaStart; j < (i * numCodons) + numCodons; j++, l++)
			{
				if (fix_dEta)
				{
					proposedCodonSpecificParameter[dEta][i][l] = currentCodonSpecificParameter[dEta][i][l];
				}
				else
				{
					proposedCodonSpecificParameter[dEta][i][l] = currentCodonSpecificParameter[dEta][i][l]
												   + covaryingNums[(numMutationCategories * numCodons) + j];
				}
			}
		}
	}
}

void ROCParameter::setProposeByPrior(bool _propose_by_prior)
{
	propose_by_prior = _propose_by_prior;
}


void ROCParameter::completeUpdateCodonSpecificParameter()
{
	for (std::string grouping : CSPToUpdate)
	{
		unsigned aaStart, aaEnd;
		SequenceSummary::AAToCodonRange(grouping, aaStart, aaEnd, true);
		unsigned aaIndex = SequenceSummary::aaToIndex.find(grouping)->second;
		numAcceptForCodonSpecificParameters[aaIndex]++;
		if (!fix_dM)
		{
			for (unsigned k = 0u; k < numMutationCategories; k++)
			{
				for (unsigned i = aaStart; i < aaEnd; i++)
				{
					currentCodonSpecificParameter[dM][k][i] = proposedCodonSpecificParameter[dM][k][i];
				}
			}
		}
		if (!fix_dEta)
		{
			for (unsigned k = 0u; k < numSelectionCategories; k++)
			{
				for (unsigned i = aaStart; i < aaEnd; i++)
				{
					currentCodonSpecificParameter[dEta][k][i] = proposedCodonSpecificParameter[dEta][k][i];
				}
			}
		}
	}
	CSPToUpdate.clear();
}


void ROCParameter::updateCodonSpecificParameter(std::string grouping)
{
	//CSPToUpdate.push_back(grouping);
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(grouping, aaStart, aaEnd, true);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(grouping)->second;
	numAcceptForCodonSpecificParameters[aaIndex]++;
	if (!fix_dM)
	{
		for (unsigned k = 0u; k < numMutationCategories; k++)
		{
			for (unsigned i = aaStart; i < aaEnd; i++)
			{
				currentCodonSpecificParameter[dM][k][i] = proposedCodonSpecificParameter[dM][k][i];
			}
		}
	}
	if (!fix_dEta)
	{
		for (unsigned k = 0u; k < numSelectionCategories; k++)
		{
			for (unsigned i = aaStart; i < aaEnd; i++)
			{
				currentCodonSpecificParameter[dEta][k][i] = proposedCodonSpecificParameter[dEta][k][i];
			}
		}
	}
}

//-------------------------------------//
//---------- Prior Functions ----------//
//-------------------------------------//

std::vector<std::vector<double>> ROCParameter::getMutationPriorMean()
{
	return mutation_prior_mean;
}

std::vector<std::vector<double>> ROCParameter::getMutationPriorStandardDeviation()
{
	return mutation_prior_sd;
}


std::vector<double> ROCParameter::getMutationPriorMeanForCategory(unsigned category)
{
	return mutation_prior_mean[category];
}

std::vector<double> ROCParameter::getMutationPriorStandardDeviationForCategory(unsigned category)
{
	return mutation_prior_sd[category];
}

void ROCParameter::getMutationPriorMeanForCategoryForGroup(unsigned category, std::string aa, double *returnSet)
{
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
	std::vector<double> mutation_prior_mean_category = mutation_prior_mean[category];
	unsigned j = 0u;
	for (unsigned i = aaStart; i < aaEnd; i++, j++)
	{
		returnSet[j] = mutation_prior_mean_category[i];
	}
}


void ROCParameter::getMutationPriorStandardDeviationForCategoryForGroup(unsigned category, std::string aa, double *returnSet)
{
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
	std::vector<double> mutation_prior_sd_category = mutation_prior_sd[category];
	unsigned j = 0u;
	for (unsigned i = aaStart; i < aaEnd; i++, j++)
	{
		returnSet[j] = mutation_prior_sd_category[i];
	}
}


void ROCParameter::setMutationPriorMean(std::vector<std::vector<double>> _mutation_prior_mean)
{
	mutation_prior_mean = _mutation_prior_mean;
}


void ROCParameter::setMutationPriorStandardDeviation(std::vector<std::vector<double>> _mutation_prior_sd)
{
	mutation_prior_sd = _mutation_prior_sd;
}





//------------------------------------------------------------------//
//---------- Posterior, Variance, and Estimates Functions ----------//
//------------------------------------------------------------------//


// double ROCParameter::getNoiseOffsetPosteriorMean(unsigned index, unsigned samples)
// {
// 	double posteriorMean = 0.0;
// 	std::vector<double> NoiseOffsetTrace = traces.getSynthesisOffsetTrace(index);
// 	unsigned traceLength = lastIteration;

// 	if (samples > traceLength)
// 	{
// 		my_printError("Warning in ROCParameter::getNoiseOffsetPosteriorMean throws: Number of anticipated samples ");
// 		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n", samples, traceLength);

// 		samples = traceLength;
// 	}
// 	unsigned start = traceLength - samples;

// 	for (unsigned i = start; i < traceLength; i++)
// 		posteriorMean += NoiseOffsetTrace[i];

// 	return posteriorMean / (double)samples;
// }


// double ROCParameter::getNoiseOffsetVariance(unsigned index, unsigned samples, bool unbiased)
// {
// 	std::vector<double> NoiseOffsetTrace = traces.getSynthesisOffsetTrace(index);
// 	unsigned traceLength = lastIteration;
// 	if (samples > traceLength)
// 	{
// 		my_printError("Warning in ROCParameter::getNoiseOffsetVariance throws: Number of anticipated samples ");
// 		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n", samples, traceLength);

// 		samples = traceLength;
// 	}
// 	double posteriorMean = getNoiseOffsetPosteriorMean(index, samples);

// 	double posteriorVariance = 0.0;

// 	unsigned start = traceLength - samples;
// 	for (unsigned i = start; i < traceLength; i++)
// 	{
// 		double difference = NoiseOffsetTrace[i] - posteriorMean;
// 		posteriorVariance += difference * difference;
// 	}
// 	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
// 	return normalizationTerm * posteriorVariance;
// }





//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


// void ROCParameter::adaptNoiseOffsetProposalWidth(unsigned adaptationWidth, bool adapt)
// {
// 	for (unsigned i = 0; i < getNumObservedPhiSets(); i++)
// 	{
// 		double acceptanceLevel = numAcceptForNoiseOffset[i] / (double)adaptationWidth;
// 		traces.updateSynthesisOffsetAcceptanceRateTrace(i, acceptanceLevel);
// 		if (adapt)
// 		{
// 			if (acceptanceLevel < 0.2)
// 				std_NoiseOffset[i] *= 0.8;
// 			if (acceptanceLevel > 0.3)
// 				std_NoiseOffset[i] *= 1.2;

// 			numAcceptForNoiseOffset[i] = 0u;
// 		}
// 	}
// }





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


// void ROCParameter::setNumObservedPhiSets(unsigned _phiGroupings)
// {
// 	obsPhiSets = _phiGroupings;
// 	noiseOffset.resize(obsPhiSets, 0.1);
// 	noiseOffset_proposed.resize(obsPhiSets, 0.1);
// 	std_NoiseOffset.resize(obsPhiSets, 0.1);
// 	numAcceptForNoiseOffset.resize(obsPhiSets, 0);
// 	observedSynthesisNoise.resize(obsPhiSets, 1.0);
// }


void ROCParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal,
										   double *returnSet)
{
	std::vector<double> *tempSet;
	tempSet = (proposal ? &proposedCodonSpecificParameter[paramType][category] : &currentCodonSpecificParameter[paramType][category]);

	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);

	unsigned j = 0u;
	for (unsigned i = aaStart; i < aaEnd; i++, j++)
	{
		returnSet[j] = tempSet->at(i);
	}
}





//-----------------------------------------------------------------------------------------------------//
//---------------------------------------- R SECTION --------------------------------------------------//
//-----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE


//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


ROCParameter::ROCParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment,
						std::vector<unsigned> _matrix, bool splitSer) : Parameter(22)
{
	unsigned _numMixtures = _matrix.size() / 2;
	std::vector<std::vector<unsigned>> thetaKMatrix;
	thetaKMatrix.resize(_numMixtures, std::vector<unsigned> (2, 0));
	unsigned index = 0;
	for (unsigned j = 0; j < 2; j++)
	{
		for (unsigned i = 0; i < _numMixtures; i++,index++)
		{
			thetaKMatrix[i][j] = _matrix[index];
		}
	}
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, "");
	initROCParameterSet();

}

ROCParameter::ROCParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
							bool splitSer, std::string _mutationSelectionState) : Parameter(22)
{
	std::vector<std::vector<unsigned>> thetaKMatrix;
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initROCParameterSet();
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//


void ROCParameter::initCovarianceMatrix(SEXP _matrix, std::string aa)
{
	std::vector<double> tmp;
	NumericMatrix matrix(_matrix);

	for (unsigned i = 0u; i < aa.length(); i++)	aa[i] = (char)std::toupper(aa[i]);

	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa) -> second;
	unsigned numRows = matrix.nrow();
	std::vector<double> covMatrix(numRows * numRows);

	//NumericMatrix stores the matrix by column, not by row. The loop
	//below transposes the matrix when it stores it.
	unsigned index = 0;
	for (unsigned i = 0; i < numRows; i++)
	{
		for (unsigned j = i; j < numRows * numRows; j += numRows, index++)
		{
			covMatrix[index] = matrix[j];
		}
	}
	CovarianceMatrix m(covMatrix);
	m.choleskyDecomposition();
	covarianceMatrix[aaIndex] = m;
}


void ROCParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa)
{
	//TODO: separate out the R wrapper functionality and make the wrapper
	//currentMutationParameter
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;

        unsigned category = getMutationCategory(mixtureElement);
        aa[0] = (char) std::toupper(aa[0]);

        unsigned aaStart, aaEnd;
        SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        for (unsigned i = aaStart, j = 0; i < aaEnd; i++, j++)
        {
            currentCodonSpecificParameter[dM][category][i] = mutationValues[j];
        }
	}
}


void ROCParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa)
{
	//TODO: separate out the R wrapper functionality and make the wrapper
	//currentSelectionParameter
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
        mixtureElement--;
        int category = getSelectionCategory(mixtureElement);

        aa[0] = (char) std::toupper(aa[0]);

        unsigned aaStart, aaEnd;
        SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        for (unsigned i = aaStart, j = 0; i < aaEnd; i++, j++)
        {
            currentCodonSpecificParameter[dEta][category][i] = selectionValues[j];
        }
    }
}

void ROCParameter::setMutationPriorMeanR(std::vector<double> _mutation_prior_mean)
{
	unsigned index = 0;
	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		for (unsigned j = 0; j < 40; j++,index++)
		{
			mutation_prior_mean[i][j] = _mutation_prior_mean[index];
		}
	}
}



void ROCParameter::setMutationPriorStandardDeviationR(std::vector<double> _mutation_prior_sd)
{
	unsigned index = 0;
	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		for (unsigned j = 0; j < 40; j++,index++)
		{
			mutation_prior_sd[i][j] = _mutation_prior_sd[index];
		}
	}
}



// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


std::vector<std::vector<double>> ROCParameter::getProposedMutationParameter()
{
	return proposedCodonSpecificParameter[dM];
}


std::vector<std::vector<double>> ROCParameter::getCurrentMutationParameter()
{
	return currentCodonSpecificParameter[dM];
}


std::vector<std::vector<double>> ROCParameter::getProposedSelectionParameter()
{
	return proposedCodonSpecificParameter[dEta];
}


std::vector<std::vector<double>> ROCParameter::getCurrentSelectionParameter()
{
	return currentCodonSpecificParameter[dEta];
}


void ROCParameter::setProposedMutationParameter(std::vector<std::vector<double>> _proposedMutationParameter)
{
	proposedCodonSpecificParameter[dM] = _proposedMutationParameter;
}


void ROCParameter::setCurrentMutationParameter(std::vector<std::vector<double>> _currentMutationParameter)
{
	currentCodonSpecificParameter[dM] = _currentMutationParameter;
}


void ROCParameter::setProposedSelectionParameter(std::vector<std::vector<double>> _proposedSelectionParameter)
{
	proposedCodonSpecificParameter[dEta] = _proposedSelectionParameter;
}


void ROCParameter::setCurrentSelectionParameter(std::vector<std::vector<double>> _currentSelectionParameter)
{
	currentCodonSpecificParameter[dEta] = _currentSelectionParameter;
}




// ------------------------------------------------------------------//
// ---------- Posterior, Variance, and Estimates Functions ----------//
// ------------------------------------------------------------------//


#endif





std::vector<double> ROCParameter::propose(std::vector<double> currentParam, double (*proposal)(double a, double b),
		double A, std::vector<double> B)
{
	unsigned _numParam = (unsigned)currentParam.size();
	std::vector<double> proposedParam(_numParam, 0.0);
	for (unsigned i = 0; i < _numParam; i++)
	{
		proposedParam[i] = (*proposal)(A + currentParam[i], B[i]);
	}
	return proposedParam;
}
