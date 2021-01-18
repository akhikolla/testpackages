#include "include/PA/PAParameter.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//

/* PAParameter Constructor (RCPP EXPOSED)
 * Arguments: None
 * Initialize the object with the default values
*/
PAParameter::PAParameter() : Parameter()
{
	//ctor
	bias_csp = 0;
	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);
}


/* PAParameter Constructor (RCPP EXPOSED)
 * Arguments: filename
 * Initialize the object with values from a restart file.
*/
PAParameter::PAParameter(std::string filename) : Parameter(64)
{
	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);
	initFromRestartFile(filename);
	numParam = 61;
}


/* PAParameter Constructor (NOT EXPOSED)
 * Arguments: synthesis rate values (vector), number of mixtures, vector containing gene assignments, vector of vectors
 *   representation of a category matrix, boolean to tell if ser should be split, keyword for mutation/selection state.
 * Initializes the object from given values. If thetaK matrix is null or empty, the mutationSelectionState keyword
 * is used to generate the matrix.
*/
PAParameter::PAParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures,
		std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer,
		std::string _mutationSelectionState) : Parameter(64)
{
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initPAParameterSet();
}


/* PAParameter Assignment operator (NOT EXPOSED)
 * Arguments: PAParameter object
 * Assign the given PAParameter object to another.
*/
PAParameter& PAParameter::operator=(const PAParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	Parameter::operator=(rhs);

	lambdaValues = rhs.lambdaValues;

	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	return *this;
}


/* PAParameter Deconstructor (NOT EXPOSED)
 * Arguments: None
 * Standard deconstructor.
*/
PAParameter::~PAParameter()
{
	//dtor
	//TODO: Need to call Parameter's deconstructor?
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//

/* initPAParameterSet (NOT EXPOSED)
 * Arguments: None
 * Initializes the variables that are specific to the PA Parameter object. The group list is set to all codons from
 * table 1 minus the stop codons. This will be corrected in CodonTable.
*/
void PAParameter::initPAParameterSet()
{
	unsigned alphaCategories = getNumMutationCategories();
	unsigned lambdaPrimeCategories = getNumSelectionCategories();

	currentCodonSpecificParameter.resize(2);
	proposedCodonSpecificParameter.resize(2);

	currentCodonSpecificParameter[alp].resize(alphaCategories);
	proposedCodonSpecificParameter[alp].resize(alphaCategories);
	currentCodonSpecificParameter[lmPri].resize(lambdaPrimeCategories);
	proposedCodonSpecificParameter[lmPri].resize(lambdaPrimeCategories);
	lambdaValues.resize(lambdaPrimeCategories);
	numParam = 61;

	for (unsigned i = 0; i < alphaCategories; i++)
	{
		std::vector <double> tmp(numParam,1.0);
		currentCodonSpecificParameter[alp][i] = tmp;
		proposedCodonSpecificParameter[alp][i] = tmp;
	}
	for (unsigned i = 0; i < lambdaPrimeCategories; i++)
	{
		std::vector <double> tmp(numParam,0.1);
		currentCodonSpecificParameter[lmPri][i] = tmp;
		proposedCodonSpecificParameter[lmPri][i] = tmp;
		lambdaValues[i] = tmp; //Maybe we don't initialize this one? or we do it differently?
	}

	bias_csp = 0;
	std_csp.resize(numParam, 0.1);

	groupList = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
		"TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
		"ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
		"AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
		"CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
		"ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
		"AGT"};
}


/* initRFPValuesFromFile (NOT EXPOSED)
 * Arguments: filename
 * Opens a restart file to initialize RFP specific values.
 */
void PAParameter::initRFPValuesFromFile(std::string filename)
{
	std::ifstream input;
	input.open(filename.c_str());
	if (input.fail())
		my_printError("ERROR: Could not open file to initialize RFP values\n");
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
				cat = 0;
				variableName = tmp.substr(1, tmp.size() - 2);
			}
			else if (flag == 2)
			{
				my_print("here\n");
			}
			else if (flag == 3) //user comment, continue
			{
				continue;
			}
			else
			{
				std::istringstream iss;
				if (variableName == "currentAlphaParameter")
				{
					if (tmp == "***")
					{
						currentCodonSpecificParameter[alp].resize(currentCodonSpecificParameter[alp].size() + 1);
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
							currentCodonSpecificParameter[alp][cat - 1].push_back(val);
						}
					}
				}
				else if (variableName == "currentLambdaPrimeParameter")
				{
					if (tmp == "***")
					{
						currentCodonSpecificParameter[lmPri].resize(currentCodonSpecificParameter[lmPri].size() + 1);
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
							currentCodonSpecificParameter[lmPri][cat - 1].push_back(val);
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
			}
		}
	}
	input.close();

	bias_csp = 0;
	proposedCodonSpecificParameter[alp].resize(numMutationCategories);
	proposedCodonSpecificParameter[lmPri].resize(numSelectionCategories);
	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		proposedCodonSpecificParameter[alp][i] = currentCodonSpecificParameter[alp][i];
	}
	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		proposedCodonSpecificParameter[lmPri][i] = currentCodonSpecificParameter[lmPri][i];
	}
}


/* writeEntireRestartFile (NOT EXPOSED)
 * Arguments: filename
 * Takes a filename and passes it to the write functions for a restart file (basic and model specific functions).
 */
void PAParameter::writeEntireRestartFile(std::string filename)
{
	writeBasicRestartFile(filename);
	writePARestartFile(filename);
}


/* writePARestartFile (NOT EXPOSED)
 * Arguments: filename
 * Appends the RFP specific values to a restart file. writeBasicRestartFile should be called previous to this by calling
 * writeEntireRestartFile.
 */
void PAParameter::writePARestartFile(std::string filename)
{
	std::ofstream out;
	std::string output = "";
	std::ostringstream oss;
	unsigned i, j;
	out.open(filename.c_str(), std::ofstream::app);
	if (out.fail())
		my_printError("ERROR: Could not open restart file for writing\n");
	else
	{
		oss << ">currentAlphaParameter:\n";
		for (i = 0; i < currentCodonSpecificParameter[alp].size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentCodonSpecificParameter[alp][i].size(); j++)
			{
				oss << currentCodonSpecificParameter[alp][i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}

		oss << ">currentLambdaPrimeParameter:\n";
		for (i = 0; i < currentCodonSpecificParameter[lmPri].size(); i++)
		{
			oss << "***\n";
			for (j = 0; j < currentCodonSpecificParameter[lmPri][i].size(); j++)
			{
				oss << currentCodonSpecificParameter[lmPri][i][j];
				if ((j + 1) % 10 == 0)
					oss << "\n";
				else
					oss << " ";
			}
			if (j % 10 != 0)
				oss << "\n";
		}

		oss << ">std_csp:\n";
		my_print("%\n", std_csp.size());
		for (i = 0; i < std_csp.size(); i++)
		{
			oss << std_csp[i];
			if ((i + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (i % 10 != 0)
			oss << "\n";

		output = oss.str();
		out << output;
	}
	out.close();

}


/* initFromRestartFile (NOT EXPOSED)
 * Arguments: filename
 * Load Parameter values in from a restart file by calling initialization functions (basic and model specific).
 */
void PAParameter::initFromRestartFile(std::string filename)
{
	initBaseValuesFromFile(filename);
	initRFPValuesFromFile(filename);
}


/* initAllTraces (NOT EXPOSED)
 * Arguments: number of samples, number of genes
 * Initializes all traces, base traces and those specific to RFP.
 */
void PAParameter::initAllTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
	traces.initializePATrace(samples, num_genes, numMutationCategories, numSelectionCategories, numParam,
						 numMixtures, categories, (unsigned)groupList.size(),obsPhiSets,currentSynthesisRateLevel[0],mixtureAssignment, estimateSynthesisRate);
}


/* initAlpha (RCPP EXPOSED VIA WRAPPER)
 * Arguments: alpha value, mixture element, codon string (all caps)
 * Gets the category and index to index into the alpha vector by looking at the mixtureElement and codon respectively.
 * Puts the alphaValue into the indexed location.
 */
void PAParameter::initAlpha(double alphaValue, unsigned mixtureElement, std::string codon)
{
	unsigned category = getMutationCategory(mixtureElement);
	unsigned index = SequenceSummary::codonToIndex(codon);
	currentCodonSpecificParameter[alp][category][index] = alphaValue;
}


/* initLambdaPrime (RCPP EXPOSED VIA WRAPPER)
 * Arguments: lambda prime value, mixture element, codon string (all caps)
 * Gets the category and index to index into the alpha vector by looking at the mixtureElement and codon respectively.
 * Puts the lambdaPrimeValue into the indexed location.
 */
void PAParameter::initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon)
{
	unsigned category = getSelectionCategory(mixtureElement);
	unsigned index = SequenceSummary::codonToIndex(codon);
	currentCodonSpecificParameter[lmPri][category][index] = lambdaPrimeValue;
}


/* initMutationSelectionCategories (RCPP EXPOSED VIA WRAPPER)
 * Arguments: vector of file names, number of categories, parameter type to initialize
 * From a file, initialize the alpha or lambda prime values for all categories. The files vector length should
 * be the same number as numCategories.
 * TODO: Rename to initAlphaLambdaCategories. This is difficult since this function is derived from a base function
 * in Parameter.cpp.
*/
void PAParameter::initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType)
{
	std::ifstream currentFile;
	std::string tmpString;
	std::string type;

	if (paramType == PAParameter::alp)
		type = "alpha";
	else
		type = "lambda";

	//TODO: Might consider doing a size check before going through all of this.
	for (unsigned i = 0; i < numCategories; i++)
	{
		std::vector<double> temp(numParam, 0.0);

		//open the file, make sure it opens
		currentFile.open(files[i].c_str());
		if (currentFile.fail())
			my_printError("Error opening file % in the file vector.\n", i);
		else
		{
			currentFile >> tmpString; //trash the first line, no info given.

			//expecting Codon,paramType Value as the current format
			while (currentFile >> tmpString)
			{
				std::string codon = tmpString.substr(0, 3);
				std::size_t pos = tmpString.find(',', 3);
				std::string val = tmpString.substr(pos + 1, std::string::npos);
				unsigned index = SequenceSummary::codonToIndex(codon, false);
				temp[index] = std::atof(val.c_str());
			}
			unsigned altered = 0u;
			for (unsigned j = 0; j < categories.size(); j++)
			{
				if (paramType == PAParameter::alp && categories[j].delM == i)
				{
					currentCodonSpecificParameter[alp][j] = temp;
					proposedCodonSpecificParameter[alp][j] = temp;
					altered++;
				}
				else if (paramType == PAParameter::lmPri && categories[j].delEta == i)
				{
					currentCodonSpecificParameter[lmPri][j] = temp;
					proposedCodonSpecificParameter[lmPri][j] = temp;
					altered++;
				}
				if (altered == numCategories)
					break; //to not access indices out of bounds.
			}
		}
		currentFile.close();
	}
}





// --------------------------------------//
// ---------- Trace Functions -----------//
// --------------------------------------//

/* updateCodonSpecificParameterTrace (NOT EXPOSED)
 * Arguments: sample index to update, codon given as a string
 * Takes a sample as an index into the trace and will eventually convert the codon into
 * an index into the trace as well.
 */
void PAParameter::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
    traces.updateCodonSpecificParameterTraceForCodon(sample, codon, currentCodonSpecificParameter[alp], alp);
    traces.updateCodonSpecificParameterTraceForCodon(sample, codon, currentCodonSpecificParameter[lmPri], lmPri);
}





// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


/* getCurrentCodonSpecificProposalWidth (NOT EXPOSED)
 * Arguments: index into the vector
 * Returns the current codon specific proposal width for the given index.
*/
double PAParameter::getCurrentCodonSpecificProposalWidth(unsigned index)
{
	return std_csp[index];
}


/* proposeCodonSpecificParameter (NOT EXPOSED)
 * Arguments: None
 * Proposes a new alpha and lambda prime value for every category and codon.
*/
void PAParameter::proposeCodonSpecificParameter()
{
unsigned numAlpha = (unsigned)currentCodonSpecificParameter[alp][0].size();
unsigned numLambdaPrime = (unsigned)currentCodonSpecificParameter[lmPri][0].size();

	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		for (unsigned j = 0; j < numAlpha; j++)
		{
			proposedCodonSpecificParameter[alp][i][j] = std::exp( randNorm( std::log(currentCodonSpecificParameter[alp][i][j]) , std_csp[j]) );
		}
	}

	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		for (unsigned j = 0; j < numLambdaPrime; j++)
		{
			proposedCodonSpecificParameter[lmPri][i][j] = std::exp( randNorm( std::log(currentCodonSpecificParameter[lmPri][i][j]) , std_csp[j]) );
		}
	}/*
    if (std::isnan(l) || std::isnan(a)){
        div_flag = TRUE;
        bool isAlpha = isnan(a);
        my_print("First divergence is alpha %\n The Current state is:
        \n", isAlpha);
    }*/
}

/* updateCodonSpecificParameter (NOT EXPOSED)
 * Arguments: string representation of a grouping (amino acid, codon...)
 * Updates the count of accepted values for codon specific parameters and updates
 * the current value to the accepted proposed value for all codon specific parameters.
*/
void PAParameter::completeUpdateCodonSpecificParameter()
{
	for (std::string codon : CSPToUpdate)
	{
    	unsigned i = SequenceSummary::codonToIndex(codon);
		numAcceptForCodonSpecificParameters[i]++;
		for(unsigned j = 0; j < getNumMixtureElements(); j++)
		{
		    currentCodonSpecificParameter[alp][j][i] = proposedCodonSpecificParameter[alp][j][i];
		    currentCodonSpecificParameter[lmPri][j][i] = proposedCodonSpecificParameter[lmPri][j][i];
		}
	}
	CSPToUpdate.clear();
}

/* updateCodonSpecificParameter (NOT EXPOSED)
 * Arguments: string representation of a grouping (amino acid, codon...)
 * Updates the count of accepted values for codon specific parameters and updates
 * the current value to the accepted proposed value for all codon specific parameters.
*/
void PAParameter::updateCodonSpecificParameter(std::string grouping)
{
	unsigned i = SequenceSummary::codonToIndex(grouping);
	numAcceptForCodonSpecificParameters[i]++;
	for(unsigned j = 0; j < getNumMixtureElements(); j++)
	{
	    currentCodonSpecificParameter[alp][j][i] = proposedCodonSpecificParameter[alp][j][i];
	    currentCodonSpecificParameter[lmPri][j][i] = proposedCodonSpecificParameter[lmPri][j][i];
	}
	
}



// ----------------------------------------------//
// ---------- Adaptive Width Functions ----------//
// ----------------------------------------------//


/* adaptCodonSpecificParameterProposalWidth (NOT EXPOSED)
 * Arguments: adaptionWidth, last iteration (NOT USED), adapt (bool)
 * Calculates the acceptance level for each codon in the group list and updates the ratio trace. If adapt is turned on,
 * meaning true, then if the acceptance level is in a certain range we change the width.
 * NOTE: This function extends Parameter's adaptCodonSpecificParameterProposalWidth function!
 */
void PAParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt)
{
	my_print("Acceptance rate for Codon Specific Parameter\n");
	my_print("\tCodon\tAcc.Rat\n"); //Prop.Width\n";
	for (unsigned i = 0; i < groupList.size(); i++)
	{
        unsigned codonIndex = SequenceSummary::codonToIndex(groupList[i]);
		double acceptanceLevel = (double)numAcceptForCodonSpecificParameters[codonIndex] / (double)adaptationWidth;
		my_print("\t%:\t%\n", groupList[i].c_str(), acceptanceLevel);
		traces.updateCodonSpecificAcceptanceRateTrace(codonIndex, acceptanceLevel);
		if (adapt)
		{
			if (acceptanceLevel < 0.2)
				std_csp[i] *= 0.8;
			if (acceptanceLevel > 0.3)
				std_csp[i] *= 1.2;
		}
		numAcceptForCodonSpecificParameters[codonIndex] = 0u;
	}
}





// -------------------------------------//
// ---------- Other Functions ----------//
// -------------------------------------//


/* getParameterForCategory (RCPP EXPOSED VIA WRAPPER)
 * Arguments: category, parameter type, codon (as a string), where or not proposed or current
 * Gets the value for a given codon specific parameter type and codon based off of if the value needed is the
 * proposed or current one.
*/
double PAParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal)
{
	double rv;
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	rv = (proposal ? proposedCodonSpecificParameter[paramType][category][codonIndex] : currentCodonSpecificParameter[paramType][category][codonIndex]);

	return rv;
}





// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//


#ifndef STANDALONE


//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


PAParameter::PAParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer) : Parameter(64)
{
  unsigned _numMixtures = _matrix.size() / 2;
  std::vector<std::vector<unsigned>> thetaKMatrix;
  thetaKMatrix.resize(_numMixtures);

	for (unsigned i = 0; i < _numMixtures; i++)
	{
		std::vector<unsigned> temp(2, 0);
		thetaKMatrix[i] = temp;
	}


  unsigned index = 0;
  for (unsigned j = 0; j < 2; j++)
  {
	for (unsigned i = 0; i < _numMixtures; i++,index++)
	{
			thetaKMatrix[i][j] = _matrix[index];
	}
  }
  initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, "");
  initPAParameterSet();

}


PAParameter::PAParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer, std::string _mutationSelectionState) :
Parameter(64)
{
  std::vector<std::vector<unsigned>> thetaKMatrix;
  initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
  initPAParameterSet();
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//


void PAParameter::initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		codon[0] = (char)std::toupper(codon[0]);
		codon[1] = (char)std::toupper(codon[1]);
		codon[2] = (char)std::toupper(codon[2]);

		initAlpha(alphaValue, mixtureElement, codon);
	}
}


void PAParameter::initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		codon[0] = (char)std::toupper(codon[0]);
		codon[1] = (char)std::toupper(codon[1]);
		codon[2] = (char)std::toupper(codon[2]);

		initLambdaPrime(lambdaPrimeValue, mixtureElement, codon);
	}
}


void PAParameter::initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories,
													std::string paramType)
{
	unsigned value = 0;
	bool check = true;
	if (paramType == "Alpha")
	{
		value = PAParameter::alp;
	}
	else if (paramType == "LambdaPrime")
	{
		value = PAParameter::lmPri;
	}
	else
	{
		my_printError("Bad paramType given. Expected \"Alpha\" or \"LambdaPrime\".\nFunction not being executed!\n");
		check = false;
	}
	if (files.size() != numCategories) //we have different sizes and need to stop
	{
		my_printError("The number of files given and the number of categories given differ. Function will not be executed!\n");
		check = false;
	}

	if (check)
	{
		initMutationSelectionCategories(files, numCategories, value);
	}
}

// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


std::vector<std::vector<double>> PAParameter::getProposedAlphaParameter()
{
	return proposedCodonSpecificParameter[alp];
}


std::vector<std::vector<double>> PAParameter::getProposedLambdaPrimeParameter()
{
	return proposedCodonSpecificParameter[lmPri];
}


std::vector<std::vector<double>> PAParameter::getCurrentAlphaParameter()
{
	return currentCodonSpecificParameter[alp];
}


std::vector<std::vector<double>> PAParameter::getCurrentLambdaPrimeParameter()
{
	return currentCodonSpecificParameter[lmPri];
}


void PAParameter::setProposedAlphaParameter(std::vector<std::vector<double>> alpha)
{
	proposedCodonSpecificParameter[alp] = alpha;
}


void PAParameter::setProposedLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime)
{
	proposedCodonSpecificParameter[lmPri] = lambdaPrime;
}


void PAParameter::setCurrentAlphaParameter(std::vector<std::vector<double>> alpha)
{
	currentCodonSpecificParameter[alp] = alpha;
}


void PAParameter::setCurrentLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime)
{
	currentCodonSpecificParameter[lmPri] = lambdaPrime;
}


// -------------------------------------//
// ---------- Other Functions ----------//
// -------------------------------------//


double PAParameter::getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal)
{
	double rv = 0.0;
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		unsigned category = 0;
		codon[0] = (char)std::toupper(codon[0]);
		codon[1] = (char)std::toupper(codon[1]);
		codon[2] = (char)std::toupper(codon[2]);
		if (paramType == PAParameter::alp)
		{
			//TODO THIS NEEDS TO CHANGE, NAMING!!!!
			category = getMutationCategory(mixtureElement); //really alpha here
		}
		else if (paramType == PAParameter::lmPri)
		{
			category = getSelectionCategory(mixtureElement);
		}
		rv = getParameterForCategory(category, paramType, codon, proposal);
	}
	return rv;
}

#endif
