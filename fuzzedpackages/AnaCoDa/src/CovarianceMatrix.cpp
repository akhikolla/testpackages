#include "include/CovarianceMatrix.h"
#include "include/SequenceSummary.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


CovarianceMatrix::CovarianceMatrix()
{
    /* Initialize with numVariates = 2.
    // Equivalent to calling initCovarianceMatrix(2) */
    numVariates = 2;

    initCovarianceMatrix(numVariates);
}


CovarianceMatrix::CovarianceMatrix(unsigned _numVariates)
{
	numVariates = _numVariates;
	initCovarianceMatrix(_numVariates);
}


CovarianceMatrix::CovarianceMatrix(std::vector <double> &matrix)
{
    numVariates = (int)std::sqrt(matrix.size());
    covMatrix = matrix;
    choleskyMatrix.resize(matrix.size(), 0.0);
}


CovarianceMatrix::CovarianceMatrix(const CovarianceMatrix& other)
{
    numVariates = other.numVariates;
    covMatrix = other.covMatrix;
    choleskyMatrix = other.choleskyMatrix;
}


CovarianceMatrix& CovarianceMatrix::operator=(const CovarianceMatrix& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    numVariates = rhs.numVariates;
    covMatrix = rhs.covMatrix;
	choleskyMatrix = rhs.choleskyMatrix;
    return *this;
}


CovarianceMatrix& CovarianceMatrix::operator+(const CovarianceMatrix& rhs)
{
	std::vector<double> cov = rhs.covMatrix;
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] += cov[i];
	}
	return *this;
}


CovarianceMatrix& CovarianceMatrix::operator*(const double &value)
{
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] *= value;
	}
	return *this;
}


void CovarianceMatrix::operator*=(const double &value)
{
    for (unsigned i = 0; i < covMatrix.size(); i++)
    {
        covMatrix[i] *= value;
    }
}


bool CovarianceMatrix::operator==(const CovarianceMatrix& other) const 
{
    bool match = true;

    if (this->covMatrix != other.covMatrix) { match = false; }
    if (this->choleskyMatrix != other.choleskyMatrix) { match = false; }
    if (this->numVariates != other.numVariates) { match = false; }

    return match;
}


CovarianceMatrix::~CovarianceMatrix()
{
    //dtor
}





//--------------------------------------//
//---------- Matrix Functions ----------//
//--------------------------------------//


void CovarianceMatrix::initCovarianceMatrix(unsigned _numVariates)
{
    numVariates = _numVariates;
    unsigned vectorLength = numVariates * numVariates;
    covMatrix.resize(vectorLength);
    choleskyMatrix.resize(vectorLength);

	double diag_const = 0.01 / (double)numVariates;
    for (unsigned i = 0u; i < vectorLength; i++)
    {
        covMatrix[i] = (i % (numVariates + 1) ? 0.0 : diag_const);
        choleskyMatrix[i] = covMatrix[i];
    }
}


void CovarianceMatrix::setDiag(double val)
{
	for (unsigned i = 0u; i < covMatrix.size(); i++)
	{
		covMatrix[i] = (i % (numVariates + 1) ? covMatrix[i] : val);
	}
}


// adaptation of http://en.wikipedia.org/wiki/Cholesky_decomposition
// http://rosettacode.org/wiki/Cholesky_decomposition#C
void CovarianceMatrix::choleskyDecomposition()
{
    for (unsigned i = 0; i < numVariates; i++)
    {
        for (unsigned j = 0; j < (i + 1); j++)
        {
            double LsubstractSum = 0.0;
            for (unsigned k = 0; k < j; k++)
            {
                LsubstractSum += choleskyMatrix[i * numVariates + k] * choleskyMatrix[j * numVariates + k];
            }
            choleskyMatrix[i * numVariates + j] = (i == j) ? std::sqrt(covMatrix[i * numVariates + i] - LsubstractSum) :
                (1.0 / choleskyMatrix[j * numVariates + j]) * (covMatrix[i * numVariates + j] - LsubstractSum);
        }
    }
}


void CovarianceMatrix::printCovarianceMatrix()
{
    for (unsigned i = 0u; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0)
            my_print("\n");
        my_print("%\t", covMatrix[i]);
    }

    my_print("\n");
}


void CovarianceMatrix::printCholeskyMatrix()
{
    for (unsigned i = 0u; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0)
            my_print("\n");
        my_print("%\t", choleskyMatrix[i]);
    }

    my_print("\n");
}


std::vector<double>* CovarianceMatrix::getCovMatrix()
{
    std::vector<double> *ptr = &covMatrix;
    return ptr;
}


std::vector<double>* CovarianceMatrix::getCholeskyMatrix()
{
    std::vector<double> *ptr = &choleskyMatrix;
    return ptr;
}


int CovarianceMatrix::getNumVariates()
{
    return numVariates;
}


std::vector<double> CovarianceMatrix::transformIidNumbersIntoCovaryingNumbers(std::vector <double> iidNumbers)
{
    std::vector<double> covaryingNumbers;
    for (unsigned i = 0u; i < numVariates; i++)
    {
        double sum = 0.0;
        for (unsigned k = 0u; k < numVariates; k++)
        {
			// testing if [i * numVariates + k] or [k * numVariates + i], first option was default
            sum += choleskyMatrix[k * numVariates + i] * iidNumbers[k];
        }

        covaryingNumbers.push_back(sum);
    }
    return covaryingNumbers;
}


void CovarianceMatrix::calculateSampleCovariance(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned lastIteration)
{
	//order of codonSpecificParameterTrace: paramType, category, numParam, samples
	unsigned numParamTypesInModel = (unsigned)codonSpecificParameterTrace.size();
	std::vector<unsigned> numCategoriesInModelPerParamType(numParamTypesInModel);
	// number of categories can vary between parameter types, see selection shared, mutation shared
	for (unsigned paramType = 0; paramType < numParamTypesInModel; paramType++)
	{
		numCategoriesInModelPerParamType[paramType] = (unsigned)codonSpecificParameterTrace[paramType].size();
	}


	unsigned start = lastIteration - samples;
	
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);

	unsigned IDX = 0;
	for (unsigned paramType1 = 0; paramType1 < numParamTypesInModel; paramType1++)
	{
		unsigned numCategoriesInModel1 = numCategoriesInModelPerParamType[paramType1];
		for (unsigned category1 = 0; category1 < numCategoriesInModel1; category1++)
		{
			for (unsigned param1 = aaStart; param1 < aaEnd; param1++)
			{
				double mean1 = sampleMean(codonSpecificParameterTrace[paramType1][category1][param1], samples, lastIteration);
				for (unsigned paramType2 = 0; paramType2 < numParamTypesInModel; paramType2++)
				{
					unsigned numCategoriesInModel2 = numCategoriesInModelPerParamType[paramType2];
					for (unsigned category2 = 0; category2 < numCategoriesInModel2; category2++)
					{
						for (unsigned param2 = aaStart; param2 < aaEnd; param2++)
						{
							double mean2 = sampleMean(codonSpecificParameterTrace[paramType2][category2][param2], samples, lastIteration);
							double unscaledSampleCov = 0.0;
							for (unsigned i = start; i < lastIteration; i++)
							{
								unscaledSampleCov += (codonSpecificParameterTrace[paramType1][category1][param1][i] - mean1) * (codonSpecificParameterTrace[paramType2][category2][param2][i] - mean2);
							}
							covMatrix[IDX] = unscaledSampleCov / ((double)samples - 1.0);

							IDX++;
						}
					}
				}
			}
		}
	}
}

void CovarianceMatrix::calculateSampleCovarianceForPANSE(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string codon, unsigned samples, unsigned lastIteration)
{
    //order of codonSpecificParameterTrace: paramType, category, numParam, samples
    unsigned numParamTypesInModel = (unsigned)codonSpecificParameterTrace.size();
    std::vector<unsigned> numCategoriesInModelPerParamType(numParamTypesInModel);
    // number of categories can vary between parameter types, see selection shared, mutation shared
    for (unsigned paramType = 0; paramType < numParamTypesInModel; paramType++)
    {
        numCategoriesInModelPerParamType[paramType] = (unsigned)codonSpecificParameterTrace[paramType].size();
    }


    unsigned start = lastIteration - samples;
    
    //unsigned aaStart, aaEnd;
    //SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
    unsigned codonIndex = SequenceSummary::codonToIndex(codon);
    unsigned IDX = 0;
    for (unsigned paramType1 = 0; paramType1 < numParamTypesInModel; paramType1++)
    {
        unsigned numCategoriesInModel1 = numCategoriesInModelPerParamType[paramType1];
        for (unsigned category1 = 0; category1 < numCategoriesInModel1; category1++)
        {
            double mean1 = sampleMean(codonSpecificParameterTrace[paramType1][category1][codonIndex], samples, lastIteration,true);

            for (unsigned paramType2 = 0; paramType2 < numParamTypesInModel; paramType2++)
            {
                unsigned numCategoriesInModel2 = numCategoriesInModelPerParamType[paramType2];
                for (unsigned category2 = 0; category2 < numCategoriesInModel2; category2++)
                {               
                    double mean2 = sampleMean(codonSpecificParameterTrace[paramType2][category2][codonIndex], samples, lastIteration,true);
                    double unscaledSampleCov = 0.0;
                    for (unsigned i = start; i < lastIteration; i++)
                    {
                        unscaledSampleCov += (std::log(codonSpecificParameterTrace[paramType1][category1][codonIndex][i]) - mean1) * (std::log(codonSpecificParameterTrace[paramType2][category2][codonIndex][i]) - mean2);
                    }
                    covMatrix[IDX] = unscaledSampleCov / ((double)samples - 1.0);
                    IDX++;   
                }
            }
        }
    }
}


double CovarianceMatrix::sampleMean(std::vector<float> sampleVector, unsigned samples, unsigned lastIteration,bool log_scale)
{
	double posteriorMean = 0.0;
	unsigned start = lastIteration - samples;
	for (unsigned i = start; i < lastIteration; i++)
	{
        if (log_scale)
        {
            posteriorMean += std::log(sampleVector[i]);
        }
        else
        {
		  posteriorMean += sampleVector[i];
	    }
    }
	return posteriorMean / (double)samples;
}



// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE

void CovarianceMatrix::setCovarianceMatrix(SEXP _matrix)
{
  std::vector<double> tmp;
  NumericMatrix matrix(_matrix);
  unsigned numRows = matrix.nrow();
  covMatrix.resize(numRows * numRows, 0.0);
  numVariates = numRows;
 
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
}





//----------------------------------//
//---------- RCPP Module -----------//
//----------------------------------//


RCPP_MODULE(CovarianceMatrix_mod)
{
  class_<CovarianceMatrix>( "CovarianceMatrix" )

        //Constructors & Destructors:
		.constructor("Empty Constructor")



		//Matrix Functions:
		.method("choleskyDecomposition", &CovarianceMatrix::choleskyDecomposition)
		.method("printCovarianceMatrix", &CovarianceMatrix::printCovarianceMatrix)
		.method("printCholeskyMatrix", &CovarianceMatrix::printCholeskyMatrix)
		.method("setCovarianceMatrix", &CovarianceMatrix::setCovarianceMatrix)
		;
}
#endif
