#include "include/PANSE/PANSEParameter.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//

/* PANSEParameter Constructor (RCPP EXPOSED)
 * Arguments: None
 * Initialize the object with the default values
*/
PANSEParameter::PANSEParameter() : Parameter()
{
	//ctor
	bias_csp = 0;
	currentCodonSpecificParameter.resize(3);
	proposedCodonSpecificParameter.resize(3);
}


/* PANSEParameter Constructor (RCPP EXPOSED)
 * Arguments: filename
 * Initialize the object with values from a restart file.
*/
PANSEParameter::PANSEParameter(std::string filename) : Parameter(61)
{
	currentCodonSpecificParameter.resize(3);
	proposedCodonSpecificParameter.resize(3);
	initFromRestartFile(filename);
	numParam = 61;
}


/* PANSEParameter Constructor (NOT EXPOSED)
 * Arguments: synthesis rate values (vector), number of mixtures, vector containing gene assignments, vector of vector
 * representation of a category matrix, boolean to tell if ser should be split, keyword for mutation/selection state.
 * Initializes the object from given values. If thetaK matrix is null or empty, the mutationSelectionState keyword
 * is used to generate the matrix.
*/
PANSEParameter::PANSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures,
		std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer,
		std::string _mutationSelectionState) : Parameter(61)
{
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initPANSEParameterSet();
}


/* PANSEParameter Assignment operator (NOT EXPOSED)
 * Arguments: PANSEParameter object
 * Assign the given PANSEParameter object to another.
*/
PANSEParameter& PANSEParameter::operator=(const PANSEParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	Parameter::operator=(rhs);

	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;
	covarianceMatrix = rhs.covarianceMatrix;
	std_partitionFunction = rhs.std_partitionFunction;
	partitionFunction_proposed = rhs.partitionFunction_proposed;
	partitionFunction = rhs.partitionFunction;


	return *this;
}


/* PANSEParameter Deconstructor (NOT EXPOSED)
 * Arguments: None
 * Standard deconstructor.
*/
PANSEParameter::~PANSEParameter()
{
	//dtor
	//TODO: Need to call Parameter's deconstructor?
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//

/* initPANSEParameterSet (NOT EXPOSED)
 * Arguments: None
 * Initializes the variables that are specific to the PANSE Parameter object. The group list is set to all codons from
 * table 1 minus the stop codons. This will be corrected in CodonTable.
*/
void PANSEParameter::initPANSEParameterSet()
{
	unsigned alphaCategories = getNumMutationCategories();
	unsigned lambdaPrimeCategories = getNumSelectionCategories();
	unsigned nonsenseErrorCategories = getNumMutationCategories();
	unsigned partitionFunctionCategories = numMixtures;

	numAcceptForPartitionFunction = 0u;

	currentCodonSpecificParameter.resize(3);
	proposedCodonSpecificParameter.resize(3);

	currentCodonSpecificParameter[alp].resize(alphaCategories);
	proposedCodonSpecificParameter[alp].resize(alphaCategories);
	currentCodonSpecificParameter[lmPri].resize(lambdaPrimeCategories);
	proposedCodonSpecificParameter[lmPri].resize(lambdaPrimeCategories);
	currentCodonSpecificParameter[nse].resize(nonsenseErrorCategories);
	proposedCodonSpecificParameter[nse].resize(nonsenseErrorCategories);
	
	partitionFunction_proposed.resize(partitionFunctionCategories, 1.0);
	partitionFunction.resize(partitionFunctionCategories, 1.0);

	numParam = 61;

	numAcceptForNSERates.resize(numParam,0u);


	for (unsigned i = 0; i < alphaCategories; i++)
	{
		std::vector <double> tmp(numParam,1.0);
		currentCodonSpecificParameter[alp][i] = tmp;
		proposedCodonSpecificParameter[alp][i] = tmp;
	}
	for (unsigned i = 0; i < lambdaPrimeCategories; i++)
	{
		std::vector <double> tmp(numParam,1.0);
		currentCodonSpecificParameter[lmPri][i] = tmp;
		proposedCodonSpecificParameter[lmPri][i] = tmp;
	}

    for (unsigned i = 0; i < nonsenseErrorCategories; i++)
    {
        std::vector <double> tmp(numParam,0.00005);
        currentCodonSpecificParameter[nse][i] = tmp;
        proposedCodonSpecificParameter[nse][i] = tmp;
    }
    for (unsigned i = 0; i < numParam; i++)
  	{
    	//CovarianceMatrix m((numMutationCategories+numSelectionCategories+numMutationCategories));
    	CovarianceMatrix m((numMutationCategories+numSelectionCategories));
    	m.choleskyDecomposition();
    	covarianceMatrix.push_back(m);
  	}

 
	bias_csp = 0;
	std_csp.resize(numParam,0.1);
	std_nse.resize(numParam,0.1);
	std_partitionFunction = 0.1;

	groupList = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
		"TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
		"ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
		"AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
		"CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
		"ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
		"AGT"};

}


/* initPANSEValuesFromFile (NOT EXPOSED)
 * Arguments: filename
 * Opens a restart file to initialize PANSE specific values.
 */
void PANSEParameter::initPANSEValuesFromFile(std::string filename)
{
	std::ifstream input;
	input.open(filename.c_str());
	std::vector <double> mat;
	if (input.fail())
		my_printError("ERROR: Could not open file to initialize PANSE values\n");
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
				else if (variableName == "currentNSERateParameter")
                {
                    if (tmp == "***")
                    {
                        currentCodonSpecificParameter[nse].resize(currentCodonSpecificParameter[nse].size() + 1);
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
                            currentCodonSpecificParameter[nse][cat - 1].push_back(val);
                        }
                    }
                }
          
                else if (variableName == "partitionFunction")
				{
					partitionFunction.resize(0);
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						partitionFunction.push_back(val);
					}
				}
                else if (variableName == "std_partitionFunction")
				{
					iss.str(tmp);
					iss >> std_partitionFunction ;
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
				else if (variableName == "std_nse")
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						std_nse.push_back(val);
					}
				}
				else if (variableName == "total_y")
				{
					iss.str(tmp);
					iss >> Y;
				}
			}
		}
	}
	input.close();

	bias_csp = 0;
	proposedCodonSpecificParameter[alp].resize(numMutationCategories);
	proposedCodonSpecificParameter[lmPri].resize(numSelectionCategories);
    proposedCodonSpecificParameter[nse].resize(numMutationCategories);
  
    partitionFunction_proposed.resize(partitionFunction.size());
    for (unsigned i = 0; i < partitionFunction.size(); i++)
    {
    	partitionFunction_proposed[i] = partitionFunction[i];
    }
    numAcceptForPartitionFunction = 0u;

	groupList = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
		"TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
		"ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
		"AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
		"CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
		"ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
		"AGT"};
	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		proposedCodonSpecificParameter[alp][i] = currentCodonSpecificParameter[alp][i];
	}
	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		proposedCodonSpecificParameter[lmPri][i] = currentCodonSpecificParameter[lmPri][i];
	}

    for (unsigned i = 0; i < numMutationCategories; i++)
    {
        proposedCodonSpecificParameter[nse][i] = currentCodonSpecificParameter[nse][i];
    }
    numAcceptForNSERates.resize(groupList.size(),0u);

}


/* writeEntireRestartFile (NOT EXPOSED)
 * Arguments: filename
 * Takes a filename and passes it to the write functions for a restart file (basic and model specific functions).
 */
void PANSEParameter::writeEntireRestartFile(std::string filename)
{
	writeBasicRestartFile(filename);
	writePANSERestartFile(filename);
}


/* writePANSERestartFile (NOT EXPOSED)
 * Arguments: filename
 * Appends the PANSE specific values to a restart file. writeBasicRestartFile should be called previous to this by calling
 * writeEntireRestartFile.
 */
void PANSEParameter::writePANSERestartFile(std::string filename)
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
        oss << ">currentNSERateParameter:\n";
        for (i = 0; i < currentCodonSpecificParameter[nse].size(); i++)
        {
            oss << "***\n";
            for (j = 0; j < currentCodonSpecificParameter[nse][i].size(); j++)
            {
                oss << currentCodonSpecificParameter[nse][i][j];
                if ((j + 1) % 10 == 0)
                    oss << "\n";
                else
                    oss << " ";
            }
            if (j % 10 != 0)
                oss << "\n";
        }
        oss << ">partitionFunction:\n";
		for (i = 0; i < partitionFunction.size(); i++)
		{
			oss << partitionFunction[i];
			if ((i + 1) % 10 == 0) oss << "\n";
			else oss <<" ";
		}
		if (i % 10 != 0) oss << "\n";
		oss << ">total_y:\n" << Y << "\n";
		oss << ">std_partitionFunction:\n" << std_partitionFunction << "\n";
		oss << ">std_csp:\n";
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
		oss << ">std_nse:\n";
		for (i = 0; i < std_nse.size(); i++)
		{
			oss << std_nse[i];
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
void PANSEParameter::initFromRestartFile(std::string filename)
{
	initBaseValuesFromFile(filename);
	initPANSEValuesFromFile(filename);
}


/* initAllTraces (NOT EXPOSED)
 * Arguments: number of samples, number of genes
 * Initializes all traces, base traces and those specific to PANSE.
 */
void PANSEParameter::initAllTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
	traces.initializePANSETrace(samples, num_genes, numMutationCategories, numSelectionCategories, numParam,
						 numMixtures, categories, (unsigned)groupList.size(),obsPhiSets,currentSynthesisRateLevel[0],mixtureAssignment,estimateSynthesisRate);
}


/* initAlpha (RCPP EXPOSED VIA WRAPPER)
 * Arguments: alpha value, mixture element, codon string (all caps)
 * Gets the category and index to index into the alpha vector by looking at the mixtureElement and codon respectively.
 * Puts the alphaValue into the indexed location.
 */
void PANSEParameter::initAlpha(double alphaValue, unsigned mixtureElement, std::string codon)
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
void PANSEParameter::initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon)
{
	unsigned category = getSelectionCategory(mixtureElement);
	unsigned index = SequenceSummary::codonToIndex(codon);
	currentCodonSpecificParameter[lmPri][category][index] = lambdaPrimeValue;
}


/* initNonsenseErrorRate(RCPP EXPOSED VIA WRAPPER)
 * Arguments: Nonsense Error Rate value, mixture element, codon string (all caps)
 * Gets the category and index to index into the alpha vector by looking at the mixtureElement and codon respectively.
 * Puts the nonsenseErrorRate into the indexed location.
 */
void PANSEParameter::initNonsenseErrorRate(double nonsenseErrorRateValue, unsigned mixtureElement, std::string codon)
{
    unsigned category = getMutationCategory(mixtureElement);
    unsigned index = SequenceSummary::codonToIndex(codon);
    currentCodonSpecificParameter[nse][category][index] = nonsenseErrorRateValue;
    //nse_rates[index] = nonsenseErrorRateValue;
}

/* initMutationSelectionCategories (RCPP EXPOSED VIA WRAPPER)
 * Arguments: vector of file names, number of categories, parameter type to initialize
 * From a file, initialize the alpha or lambda prime values for all categories. The files vector length should
 * be the same number as numCategories.
*/
void PANSEParameter::initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType)
{
	std::ifstream currentFile;
	std::string tmpString;
	std::string type;


	if (paramType == PANSEParameter::alp)
		type = "alpha";
	else if (paramType == PANSEParameter::nse)
	    type = "nse";
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
				if (paramType == PANSEParameter::alp && categories[j].delM == i)
				{
					currentCodonSpecificParameter[alp][j] = temp;
					proposedCodonSpecificParameter[alp][j] = temp;
					altered++;
				}
				else if (paramType == PANSEParameter::lmPri && categories[j].delEta == i)
				{
					currentCodonSpecificParameter[lmPri][j] = temp;
					proposedCodonSpecificParameter[lmPri][j] = temp;
					altered++;
				}
                else if (paramType == PANSEParameter::nse && categories[j].delM == i)
                {
                    currentCodonSpecificParameter[nse][j] = temp;
                    proposedCodonSpecificParameter[nse][j] = temp;
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
void PANSEParameter::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
	traces.updateCodonSpecificParameterTraceForCodon(sample, codon, currentCodonSpecificParameter[alp], alp);
	traces.updateCodonSpecificParameterTraceForCodon(sample, codon, currentCodonSpecificParameter[lmPri], lmPri);
    traces.updateCodonSpecificParameterTraceForCodon(sample, codon, currentCodonSpecificParameter[nse], nse);
}

void PANSEParameter::updatePartitionFunctionTrace(unsigned sample)
{
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        traces.updatePartitionFunctionTrace(i, sample, partitionFunction[i]);
    }
}


// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


/* getCurrentCodonSpecificProposalWidth (NOT EXPOSED)
 * Arguments: index into the vector
 * Returns the current codon specific proposal width for the given index.
*/
double PANSEParameter::getCurrentCodonSpecificProposalWidth(unsigned index)
{
	return std_csp[index];
}


/* proposeCodonSpecificParameter (NOT EXPOSED)
 * Arguments: None
 * Proposes a new alpha and lambda prime value for every category and codon.
*/
void PANSEParameter::proposeCodonSpecificParameter()
{
	unsigned numAlpha = (unsigned)currentCodonSpecificParameter[alp][0].size();
	unsigned numLambdaPrime = (unsigned)currentCodonSpecificParameter[lmPri][0].size();
    unsigned numNSE = (unsigned)currentCodonSpecificParameter[nse][0].size();


	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		for (unsigned j = 0; j < numAlpha; j++)
		{
			if (fix_alpha)
			{
				proposedCodonSpecificParameter[alp][i][j] = currentCodonSpecificParameter[alp][i][j];
			}
			else
			{
				proposedCodonSpecificParameter[alp][i][j] = std::exp( randNorm( std::log(currentCodonSpecificParameter[alp][i][j]) , std_csp[j]) );
			}
		}
	}
	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		for (unsigned j = 0; j < numLambdaPrime; j++)
		{
			if (fix_lp)
			{
				proposedCodonSpecificParameter[lmPri][i][j] = currentCodonSpecificParameter[lmPri][i][j];

			}
			else
			{
				if (j == 29) // M/ATG has index 29, see SequenceSummary. This will need to be updated to be more general.
        		{
        			//This guarantees \alpha/\lambda will be 1. 
        			proposedCodonSpecificParameter[lmPri][i][j] = proposedCodonSpecificParameter[alp][i][j];
        		}
        		else
        		{
					proposedCodonSpecificParameter[lmPri][i][j] = std::exp( randNorm( std::log(currentCodonSpecificParameter[lmPri][i][j]) , std_csp[j]) );
				}
			}
		}
	}

    for (unsigned i = 0; i < numMutationCategories; i++)
    {
        for (unsigned j = 0; j < numNSE; j++)
        {
        	if (fix_nse)
        	{
        		proposedCodonSpecificParameter[nse][i][j] = currentCodonSpecificParameter[nse][i][j];
        	}
        	else
        	{
        		
        		proposedCodonSpecificParameter[nse][i][j] = std::exp( randNorm( std::log(currentCodonSpecificParameter[nse][i][j]) , std_nse[j]) );
        		
        	}
    	}
    }
	
}




void PANSEParameter::updateCodonSpecificParameter(std::string grouping)
{
	my_print("ERROR: PANSEParameter::updateCodonSpecificParameterstd::string grouping) should not be called for PANSE.\n");
}




/* updateCodonSpecificParameter (NOT EXPOSED)
 * Arguments: string representation of a grouping (amino acid, codon...)
 * Updates the count of accepted values for codon specific parameters and updates
 * the current value to the accepted proposed value for all codon specific parameters.
*/
void PANSEParameter::updateCodonSpecificParameter(std::string grouping,std::string param)
{
	unsigned i = SequenceSummary::codonToIndex(grouping);
	if (param == "Elongation")
	{
		numAcceptForCodonSpecificParameters[i]++;
		//my_print("% % % %\n",grouping,i,param,numAcceptForCodonSpecificParameters[i]);
	    for (unsigned k = 0u; k < numMutationCategories; k++)
	    {
	        currentCodonSpecificParameter[alp][k][i] = proposedCodonSpecificParameter[alp][k][i];
	    }
	    for (unsigned k = 0u; k < numSelectionCategories; k++)
	    {
	        currentCodonSpecificParameter[lmPri][k][i] = proposedCodonSpecificParameter[lmPri][k][i];
	    }
	}
	else
	{
		numAcceptForNSERates[i]++;
		//my_print("% % % %\n",grouping,i,param,numAcceptForNSERates[i]);
	    for (unsigned k = 0u; k < numMutationCategories; k++)
	    {
	        currentCodonSpecificParameter[nse][k][i] = proposedCodonSpecificParameter[nse][k][i];
	    }
	}
}


/* updateCodonSpecificParameter (NOT EXPOSED)
 * Arguments: string representation of a grouping (amino acid, codon...)
 * Updates the count of accepted values for codon specific parameters and updates
 * the current value to the accepted proposed value for all codon specific parameters.
*/
void PANSEParameter::completeUpdateCodonSpecificParameter()
{
	for (std::string codon : CSPToUpdate)
    {
    	unsigned i = SequenceSummary::codonToIndex(codon);
        numAcceptForCodonSpecificParameters[i]++;

        for (unsigned k = 0u; k < numMutationCategories; k++)
        {
            currentCodonSpecificParameter[alp][k][i] = proposedCodonSpecificParameter[alp][k][i];
        }
        for (unsigned k = 0u; k < numSelectionCategories; k++)
        {
            currentCodonSpecificParameter[lmPri][k][i] = proposedCodonSpecificParameter[lmPri][k][i];
        }
        for (unsigned k = 0u; k < numMutationCategories; k++)
        {
            currentCodonSpecificParameter[nse][k][i] = proposedCodonSpecificParameter[nse][k][i];
        }
    }
    CSPToUpdate.clear();
}


// ----------------------------------------------//
// -------- Partition Function Functions --------//
// ----------------------------------------------//

double PANSEParameter::getPartitionFunction(unsigned mixtureCategory, bool proposed)
{
    if (proposed)
    {
        return partitionFunction_proposed[mixtureCategory];
    }
    return partitionFunction[mixtureCategory];
}


void PANSEParameter::proposePartitionFunction()
{
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        partitionFunction_proposed[i] = std::exp( randNorm( std::log(partitionFunction[i]) , std_partitionFunction) );
    }

}


void PANSEParameter::setPartitionFunction(double newPartitionFunction, unsigned mixtureCategory)
{
    partitionFunction[mixtureCategory] = newPartitionFunction;
    partitionFunction_proposed[mixtureCategory] = newPartitionFunction;
}


double PANSEParameter::getCurrentPartitionFunctionProposalWidth()
{
    return std_partitionFunction;
}


unsigned PANSEParameter::getNumAcceptForPartitionFunction()
{
    return numAcceptForPartitionFunction;
}


void PANSEParameter::updatePartitionFunction()
{
    for (unsigned i = 0u; i < numMixtures; i++)
    {
        partitionFunction[i] = partitionFunction_proposed[i];
    }
    numAcceptForPartitionFunction++;
}

// ----------------------------------------------//
// ---------- Adaptive Width Functions ----------//
// ----------------------------------------------//


void PANSEParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt)
{
  //Gelman BDA 3rd Edition suggests a target acceptance rate of 0.23
  // for high dimensional problems
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
  double adjustFactor; //variable assigned value of either adjustFactorLow or adjustFactorHigh

  factorCriteriaLow = acceptanceTargetLow - diffFactorAdjust;  //below this value weighted sum and factor adjustments are applied
  factorCriteriaHigh = acceptanceTargetHigh + diffFactorAdjust;  //above this value weighted sum and factor adjustments are applied

  adaptiveStepPrev = adaptiveStepCurr;
  adaptiveStepCurr = lastIteration;
  //unsigned samples = adaptiveStepCurr - adaptiveStepPrev;

  my_print("Acceptance rates for Codon Specific Parameters\n");
  my_print("Target range: %-% \n", factorCriteriaLow, factorCriteriaHigh );
  my_print("Adjustment range: < % or > % \n", acceptanceTargetLow, acceptanceTargetHigh );
 //my_print("\tCodon\tElongation Acc.Rat\tNSERate Acc. Rate\tAverage Acc. Rate\n"); //Prop.Width\n";

  my_print("\tCodon\tAcc.Rat\n"); //Prop.Width\n";

  for (unsigned i = 0; i < groupList.size(); i++) //cycle through all of the aa
  {
  	std::string codon = groupList[i];
    unsigned codonIndex = SequenceSummary::codonToIndex(codon);
    double acceptanceLevel_elong = (double)numAcceptForCodonSpecificParameters[codonIndex] / (double)adaptationWidth;
    double acceptanceLevel_nse = (double)numAcceptForNSERates[codonIndex] / (double)adaptationWidth;
    double acceptanceLevel = 0.5 * (acceptanceLevel_elong+acceptanceLevel_nse);
    my_print("\t%:\t%\t%\t%\n", codon.c_str(), acceptanceLevel_elong,acceptanceLevel_nse,acceptanceLevel);
   	//my_print("\t%:\t%\n", codon.c_str(), acceptanceLevel_elong);
   
    traces.updateCodonSpecificAcceptanceRateTrace(codonIndex, acceptanceLevel_elong);
    traces.updateNseRateSpecificAcceptanceRateTrace(codonIndex, acceptanceLevel_nse);
      //Evaluate current acceptance ratio  performance
    if (acceptanceLevel_elong < factorCriteriaLow) acceptanceUnder++;
    else if (acceptanceLevel_elong > factorCriteriaHigh) acceptanceOver++;

    if (adapt)
	{
		if( (acceptanceLevel_elong < acceptanceTargetLow) || (acceptanceLevel_elong > acceptanceTargetHigh) )// adjust proposal width
	  	{

	      // define adjustFactor
	    
	      if (acceptanceLevel_elong < factorCriteriaLow)
	      {
		  	adjustFactor = adjustFactorLow;
		  	std_csp[codonIndex] *= adjustFactor;
		  }
	      else if(acceptanceLevel_elong > factorCriteriaHigh)
	      {
		  	adjustFactor = adjustFactorHigh;
		  	std_csp[codonIndex] *= adjustFactor;
		  }
		 
		}// end if statement


		if ( (acceptanceLevel_nse < acceptanceTargetLow) || (acceptanceLevel_nse > acceptanceTargetHigh))
		{
		 if (acceptanceLevel_nse < factorCriteriaLow)
	      {
		  	adjustFactor = adjustFactorLow;
		  	std_nse[codonIndex] *= adjustFactor;
		  }
	      else if(acceptanceLevel_nse > factorCriteriaHigh)
	      {
		  	adjustFactor = adjustFactorHigh;
		  	std_nse[codonIndex] *= adjustFactor;
		  }
	
		}


	} // end if(adapt)

	numAcceptForCodonSpecificParameters[codonIndex] = 0u;
	numAcceptForNSERates[codonIndex] = 0u;

   }
}




/* adaptPartitionFunctionProposalWidth (NOT EXPOSED)
 * Arguments: adaptionWidth, last iteration (NOT USED), adapt (bool)

 */
void PANSEParameter::adaptPartitionFunctionProposalWidth(unsigned adaptationWidth, bool adapt)
{
    double acceptanceLevel = (double)numAcceptForPartitionFunction / (double)adaptationWidth;
    traces.updatePartitionFunctionAcceptanceRateTrace(acceptanceLevel);
    my_print("Acceptance level for partition function: %\n", acceptanceLevel);
    if (adapt)
    {
        if (acceptanceLevel < 0.2)
            std_partitionFunction *= 0.8;
        if (acceptanceLevel > 0.3)
            std_partitionFunction *= 1.2;
    }
    numAcceptForPartitionFunction = 0u;
}


// -------------------------------------//
// ---------- Other Functions ----------//
// -------------------------------------//


/* getParameterForCategory (RCPP EXPOSED VIA WRAPPER)
 * Arguments: category, parameter type, codon (as a string), where or not proposed or current
 * Gets the value for a given codon specific parameter type and codon based off of if the value needed is the
 * proposed or current one.
*/
double PANSEParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal)
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


PANSEParameter::PANSEParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer) : Parameter(64)
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
  initPANSEParameterSet();

}


PANSEParameter::PANSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer, std::string _mutationSelectionState) :
Parameter(64)
{
  std::vector<std::vector<unsigned>> thetaKMatrix;
  initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
  initPANSEParameterSet();
}

void PANSEParameter::initCovarianceMatrix(SEXP _matrix, std::string codon)
{
	std::vector<double> tmp;
	NumericMatrix matrix(_matrix);

	for (unsigned i = 0u; i < codon.length(); i++)	codon[i] = (char)std::toupper(codon[i]);

	unsigned codonIndex = SequenceSummary::codonToIndexWithReference.find(codon) -> second;
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
	covarianceMatrix[codonIndex] = m;

}



//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//


void PANSEParameter::initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon)
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


void PANSEParameter::initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon)
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

void PANSEParameter::fixAlpha()
{
	fix_alpha = true;
}

void PANSEParameter::fixLambdaPrime()
{
	fix_lp = true;
}

void PANSEParameter::fixNSERate()
{
	fix_nse = true;
}

void PANSEParameter::shareNSERate()
{
	share_nse = true;
}



void PANSEParameter::initNSERateR(double NSERateValue, unsigned mixtureElement, std::string codon)
{
    bool check = checkIndex(mixtureElement, 1, numMixtures);
    if (check)
    {
        mixtureElement--;
        codon[0] = (char)std::toupper(codon[0]);
        codon[1] = (char)std::toupper(codon[1]);
        codon[2] = (char)std::toupper(codon[2]);

        initLambdaPrime(NSERateValue, mixtureElement, codon);
    }
}

void PANSEParameter::initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories,
													std::string paramType)
{
	unsigned value = 0;
	bool check = true;
	if (paramType == "Alpha")
	{
		value = PANSEParameter::alp;
	}
	else if (paramType == "LambdaPrime")
	{
		value = PANSEParameter::lmPri;
	}
    else if (paramType == "NSERate")
    {
        value = PANSEParameter::nse;
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


std::vector<std::vector<double>> PANSEParameter::getProposedAlphaParameter()
{
	return proposedCodonSpecificParameter[alp];
}


std::vector<std::vector<double>> PANSEParameter::getProposedLambdaPrimeParameter()
{
	return proposedCodonSpecificParameter[lmPri];
}


std::vector<std::vector<double>> PANSEParameter::getProposedNSERateParameter()
{
    return proposedCodonSpecificParameter[nse];
}


std::vector<std::vector<double>> PANSEParameter::getCurrentAlphaParameter()
{
	return currentCodonSpecificParameter[alp];
}


std::vector<std::vector<double>> PANSEParameter::getCurrentLambdaPrimeParameter()
{
	return currentCodonSpecificParameter[lmPri];
}


std::vector<std::vector<double>> PANSEParameter::getCurrentNSERateParameter()
{
    return currentCodonSpecificParameter[nse];
}


void PANSEParameter::setProposedAlphaParameter(std::vector<std::vector<double>> alpha)
{
	proposedCodonSpecificParameter[alp] = alpha;
}


void PANSEParameter::setProposedLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime)
{
	proposedCodonSpecificParameter[lmPri] = lambdaPrime;
}


void PANSEParameter::setProposedNSERateParameter(std::vector<std::vector<double>> nseRate)
{
    proposedCodonSpecificParameter[nse] = nseRate;
}


void PANSEParameter::setCurrentAlphaParameter(std::vector<std::vector<double>> alpha)
{
	currentCodonSpecificParameter[alp] = alpha;
}


void PANSEParameter::setCurrentLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime)
{
	currentCodonSpecificParameter[lmPri] = lambdaPrime;
}


void PANSEParameter::setCurrentNSERateParameter(std::vector<std::vector<double>> nseRate)
{
    currentCodonSpecificParameter[nse] = nseRate;
}


// -------------------------------------//
// ---------- Other Functions ----------//
// -------------------------------------//


double PANSEParameter::getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal)
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
		if (paramType == PANSEParameter::alp)
		{
			category = getMutationCategory(mixtureElement); //really alpha here
		}
		else if (paramType == PANSEParameter::lmPri)
		{
			category = getSelectionCategory(mixtureElement);
		}
        else if (paramType == PANSEParameter::nse)
        {
            category = getMutationCategory(mixtureElement);
        }
		rv = getParameterForCategory(category, paramType, codon, proposal);
	}
	return rv;
}
#endif

void PANSEParameter::readAlphaValues(std::string filename)
{
    std::size_t pos;
    std::ifstream currentFile;
    std::string tmpString;
    std::vector <double> tmp;

    tmp.resize(64,1);

    currentFile.open(filename);
    if (currentFile.fail())
        my_printError("Error opening file %\n", filename.c_str());
    else
    {
        currentFile >> tmpString;
        while (currentFile >> tmpString){
            pos = tmpString.find(',');
            if (pos != std::string::npos)
            {
                std::string codon = tmpString.substr(0,3);
                std::string val = tmpString.substr(pos + 1, std::string::npos);
                tmp[SequenceSummary::codonToIndex(codon)] = std::atof(val.c_str());
            }
        }
    }
    currentFile.close();

	unsigned alphaCategories = getNumMutationCategories();

	for (unsigned i = 0; i < alphaCategories; i++)
	{
		currentCodonSpecificParameter[alp][i] = tmp;
		proposedCodonSpecificParameter[alp][i] = tmp;
	}

}

void PANSEParameter::readLambdaValues(std::string filename)
{
    std::size_t pos;
    std::ifstream currentFile;
    std::string tmpString;
    std::vector <double> tmp;

    tmp.resize(64, 0.1);

    currentFile.open(filename);
    if (currentFile.fail())
        my_printError("Error opening file %\n", filename.c_str());
    else
    {
        currentFile >> tmpString;
        while (currentFile >> tmpString){
            pos = tmpString.find(',');
            if (pos != std::string::npos)
            {
                std::string codon = tmpString.substr(0,3);
                std::string val = tmpString.substr(pos + 1, std::string::npos);
                tmp[SequenceSummary::codonToIndex(codon)] = std::atof(val.c_str());
            }
        }
    }

    currentFile.close();
	unsigned lambdaPrimeCategories = getNumSelectionCategories();

	for (unsigned i = 0; i < lambdaPrimeCategories; i++)
	{
		currentCodonSpecificParameter[lmPri][i] = tmp;
		proposedCodonSpecificParameter[lmPri][i] = tmp;
	}

}

void PANSEParameter::readNSEValues(std::string filename)
{
    std::size_t pos;
    std::ifstream currentFile;
    std::string tmpString;
    std::vector <double> tmp;

    tmp.resize(64, 0.1);

    currentFile.open(filename);
    if (currentFile.fail())
        my_printError("Error opening file %\n", filename.c_str());
    else
    {
        currentFile >> tmpString;
        while (currentFile >> tmpString){
            pos = tmpString.find(',');
            if (pos != std::string::npos)
            {
                std::string codon = tmpString.substr(0,3);
                std::string val = tmpString.substr(pos + 1, std::string::npos);
                tmp[SequenceSummary::codonToIndex(codon)] = std::atof(val.c_str());
            }
        }
    }

    currentFile.close();
    unsigned NSERateCategories = getNumMutationCategories();

    for (unsigned i = 0; i < NSERateCategories; i++)
    {
        currentCodonSpecificParameter[nse][i] = tmp;
        proposedCodonSpecificParameter[nse][i] = tmp;
    }

}


std::vector<double> PANSEParameter::oneMixLambda(){
    return currentCodonSpecificParameter[lmPri][0];
}
std::vector<double> PANSEParameter::oneMixAlpha(){
    return currentCodonSpecificParameter[alp][0];
}
std::vector<double> PANSEParameter::oneMixNSE(){
    return currentCodonSpecificParameter[nse][0];
}


void PANSEParameter::setTotalRFPCount(Genome& genome)
{
	Y = genome.getSumRFP();
}

unsigned PANSEParameter::getTotalRFPCount()
{
	return Y;
}
