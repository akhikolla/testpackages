/***************************************************************************
                             SRC/mixmod/Utilities/Util.h  description
    copyright            : (C) MIXMOD Team - 2001-2016
    email                : contact@mixmod.org
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD
    
    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

    All informations available on : http://www.mixmod.org                                                                                               
***************************************************************************/
/** @file Util.h
 *  @brief Constants definitions, various utilities to describe models, and others...
 **/

#ifndef XEM_UTIL_H
#define XEM_UTIL_H

#ifndef WANT_STREAM
#define WANT_STREAM
#endif

#ifndef WANT_MATH
#define WANT_MATH
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
/** Exception Classes*/
#include "mixmod/Utilities/exceptions/Exception.h"
#include "mixmod/Utilities/exceptions/DCVException.h"
#include "mixmod/Utilities/exceptions/DCVonlyInGaussianCaseException.h"
#include "mixmod/Utilities/exceptions/InputException.h"
#include "mixmod/Utilities/exceptions/NumericException.h"
#include "mixmod/Utilities/exceptions/OtherException.h"

// Need matrix output routines
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <ctype.h>
#include <stdint.h>
#include <memory>
#include <limits>
#include <sstream>

#include <typeinfo>
#include <assert.h> //For debugging purpose

// Need matrix applications
//#include "mixmod/Utilities/maths/SelectLibrary.h"
#include <cmath>
#ifndef XEMmathLib
#define XEMmathLib 1 // default is Eigen
#endif

using namespace std;

namespace XEM {

//Macro for throw
#define THROW(Exceptionclass,errortype) throw Exceptionclass(__FILE__,__LINE__,errortype)

/** Declare global variable for NoError */
static InputException defaultException = InputException(noError);
static Exception & NOERROR = defaultException;

class ModelType;
class Matrix;
class Algo;

// Constants definitions

/// Define PI
const double XEMPI = 3.14159265358979323846; // Number pi

// misc
const int DEBUG = 0; // >0 for debug information
                     // 0 : no debug
                     // 1 : param debug
                     // 2 : param and fik, tik, zik debug
const bool DATA_REDUCE = false; // do not reduce binary data (bug with weights, I don't want to try to fix it.....)

extern int VERBOSE;    // 0 = No VERBOSE, 1 = print various execution traces (labels, errors...)
extern int MASSICCC;   // 0 = no output, 1 = output for massiccc purpose (progress file, entropy computation...), (11 = 1 for predict)
enum class IoMode // How to read/Write data in files ? (binary or ascii format)
{
	NUMERIC,
	BINARY
};
extern IoMode IOMODE; // BINARY for exact values (in hexa) in ouputs files for non regression tests; ASCII for human-readable.

//write floating point number into stream (usually a file)
void putDoubleInStream(std::ostream& output, double value, std::string appendChars = "");

// retrieve floating point number from stream (usually a file)
double getDoubleFromStream(std::istream& input);

/// Define number of maximum samples
const int64_t maxNbSample = 1000000;     // Maximum sample size
const int64_t maxPbDimension = 10000;    // Maximum sample dimension
const int64_t minNbIteration = 1;        // Minimum number of iterations
const int64_t minNbIterationForSEM = 50; // Minimum number of iterations for SEM
const int64_t maxNbIteration = 100000;   // Maximum number of iterations
const int64_t defaultNbIteration = 200;  // Default number of iteration
const double minEpsilon = 0;             // Minimum value of eps
const double maxEpsilon = 1;             // Maximum value of eps

const double defaultEpsilon = 1.0E-3; // Default value of eps
const int64_t maxNbNbCluster = 10;    // Maximum size of cluster list
const int64_t maxNbAlgo = 5;          // Maximum number of algorithms
const int64_t defaultNbAlgo = 1;      // Default number of algorithms
const int64_t maxNbStrategy = 10;     // Maximum number of strategies
const int64_t defaultNbStrategy = 1;  // Default number of strategies
const int64_t maxNbModel = 100;       // Maximum number of models
const int64_t defaultNbModel = 1;     // Default number of models
const int64_t maxNbCriterion = 4;     // Maximum number of criteria
const int64_t defaultNbCriterion = 1; // Default number of criteria

const int64_t minNbTryInStrategy = 1;             // min of strategies repeats
const int64_t maxNbTryInStrategy = 100;           // max of strategies repeats
const int64_t defaultNbTryInStrategy = 1;           // number of strategies repeats
const int64_t nbTryInDefaultClusteringStrategy = 1; // number of strategies repeats

const int64_t minNbTryInInit = 1;      // min of repeats in init
const int64_t maxNbTryInInit = 1000;   // max of repeats in init
const int64_t defaultNbTryInInit = 10; // number of repeats in init

const int64_t minNbIterationInInit = 1;                // min number of iterations in init
const int64_t maxNbIterationInInit = 1000;             // max number of iterations in init
const int64_t defaultNbIterationInInit = 5;            // default number of iterations in init
const int64_t defaultNbIterationInInitForSemMax = 100; // default number of iterations in init

const double minEpsilonInInit = 0;         // min number of iterations in init
const double maxEpsilonInInit = 1;         // max number of iterations in init
const double defaultEpsilonInInit = 0.001; // default number of iterations in init

const int64_t maxNbIterationInCEM_INIT = 100; // Maximum number of iterations of CEM in CEM_INIT

const double minOverflow = std::numeric_limits<double>::min();  // Minimum value for overflow
const double minUnderflow = std::numeric_limits<double>::min(); // Minimum value for underflow
const int64_t nbMaxSelection  = 5;              // Maximum number of selection
const int64_t maxNbOutputFiles = 52;            // Maximum number of output Files
const int64_t nbTestOutputFiles = 7;            // Number of output files to compare in test
const double defaultFluryEpsilon = 0.001;       // default value for espilon in flury algorthm
const int64_t maxFluryIter = 7;                 // maximum of number of Flury iterations
const double minDeterminantValue =
		std::numeric_limits<double>::min();     // minimum value of determinant of sigma
const double maxRelativeDiffValueTest = 1.0E-5; // Maximum difference between 2 value in test
const double maxAbsoluteDiffValueTest = 1.0E-8; // Maximum difference between 2 value in test

const int64_t  defaultDCVnumberOfBlocks = 10;   // DCV
const int64_t defaultCVnumberOfBlocks  = 10;    // CV

const double minValueForLLandLLOne = 1.e-10;    // minimum value for LL - LLone

const int64_t int64_t_max = std::numeric_limits<int64_t>::max();

const int64_t nbQualitativeGraphics = 2;
const int64_t nbQuantitativeGraphics = 3;

/*
Notes :
Enumeration types will be called ...Name
Ex : StrategyInitName
 */

enum CVinitBlocks {

	CV_RANDOM = 0, // initialize the CV blocks by random

	CV_DIAG = 1    // initialize the CV blocks by assiging
	/*
	sample 1 : for w=1 to its weight : sample 1 in block w
	sample 2 : for w=1 to its weight : sample w1+1 in block w1+w

	Ex : 1
	//-----
	ind  weight
	1     1
	2     1
	3     1
	...
	ind 1 -> block1
	ind 2 -> block2
	...
	ind V -> blockV
	ind V+1->block1
	...

	Ex 2 :
	//-----
	ind weight
	1     2
	2     2
	3     2
	4     2
	5     2

	if V=4:
	ind 1 -> block1
	ind 1 -> block2
	ind 2 -> block3
	ind 2 -> block4
	ind 3 -> block1
	ind 3 -> block2
	ind 4 -> block3
	ind 4 -> block4
	ind 5 -> block1
	ind 5 -> block2

	->>  block 1 : ind 1 - 3 - 5
		 block 2 : ind 1 - 3 - 5
		 block 3 : ind 2 - 4
		 block 4 : ind 4 - 5
	 */
};
const CVinitBlocks defaultCVinitBlocks = CV_RANDOM;

enum DCVinitBlocks {

	DCV_RANDOM = 0, // initialize the DCV blocks by random

	DCV_DIAG = 1    // initialize the DCV blocks by the same way that CV_DIAG
	/*
	Ex : 1
	//-----
	ind  weight
	1     1
	2     1
	3     1
	...
	ind 1 -> blockTest1
	ind 2 -> blockTest2
	...
	ind V -> blockTestV
	ind V+1->blockTest1
	...

	Ex 2 :
	//-----
	ind weight
	1     10
	2     10
	3     10
	4     10
	5     10

	if V=4:
	ind 1 -> blockTest1
	ind 1 -> blockTest2
	ind 1 -> blockTest3
	ind 1 -> blockTest4
	ind 1 -> blockTest1
	ind 1 -> blockTest2
	ind 1 -> blockTest3
	ind 1 -> blockTest4
	ind 1 -> blockTest1
	ind 1 -> blockTest2
	ind 2 -> blockTest3
	ind 2 -> blockTest4
	ind 2 -> blockTest1
	ind 2 -> blockTest2
	ind 2 -> blockTest3
	ind 2 -> blockTest4
	ind 2 -> blockTest1
	ind 2 -> blockTest2
	ind 2 -> blockTest3
	ind 2 -> blockTest4

	->>  blockTest 1 : ind 1(x3) - 2(x2) - 3(x3) - 4(x2) - 5(x3) - 6(x2) - 7(x3) - 8(x2) - 9(x3) - 10(x2)
	blockLearning 1 : ind 1(x7) - 2(x8) - 3(x7) - 4(x8) - 5(x7) - 6(x8) - 7(x7) - 8(x8) - 9(x7) - 10(x8)
	...
	 */
};
const DCVinitBlocks defaultDCVinitBlocks = DCV_RANDOM;

/** @enum StrategyInitName
	@brief Enumeration of differents strategy initialization
 */
enum StrategyInitName {

	RANDOM = 0,         // Random centers
	USER = 1,           // Initial parameters specified by user
	USER_PARTITION = 2, // Partition specified by user
	SMALL_EM = 3,       // EM strategy for initial parameters
	CEM_INIT = 4,       // initialization with CEM
	SEM_MAX = 5         // initialization with SEM max
};
const StrategyInitName defaultStrategyInitName = SMALL_EM;
const StrategyInitName defaultClusteringStrategyInitName = SMALL_EM;

// Type of convergence for each algorithm

/** @enum AlgoStopName
	@brief Enumeration of differents type of converge of algorithm (stop rule)
 */
enum AlgoStopName {

	NO_STOP_NAME = -1,      // for MAP or M algo
	NBITERATION = 0,        // Number of iterations specified by user
	EPSILON = 1,            // Stationarity of the xml criterion
	NBITERATION_EPSILON = 2 // Number of iterations & xml criterion
};
const AlgoStopName defaultAlgoStopName = NBITERATION_EPSILON;

/** @enum CriterionName
@brief Enumeration of Criterion type
 */
enum CriterionName {

	UNKNOWN_CRITERION_NAME = -1, // Unknown criterion

	BIC = 0,                     // Bayesian information criterion
	CV  = 1,                     // Cross validation criterion
	ICL = 2,                     // Integrated completed likelihood
	NEC = 3,                     // Entropy criterion
	DCV = 4                      // Double Cross validation criterion
};
const CriterionName defaultCriterionName = BIC;
const CriterionName defaultLearnCriterionName = CV;

/** @enum AlgoName
@brief Enumeration of Algo type
 */
enum AlgoName {

	UNKNOWN_ALGO_NAME = -1, // Unknown algorithm
	MAP = 0,                // Maximum a posteriori
	EM = 1,                 // Expectation maximization
	CEM = 2,                // Classification EM
	SEM = 3,                // Stochastic EM
	M = 4                   // Maximization
};
const AlgoName defaultAlgoName = EM;
const AlgoName defaultClusteringAlgoName = EM;

/** @enum FormatFile
	@brief Format of differents data file
 */
const int64_t nbFormatNumeric = 3;
namespace FormatNumeric {

	enum FormatNumericFile {

		txt = 0,  // Format txt (ascii)
		hdf5 = 1, // Format hdf5
		XML = 2   // Format XML
	};
	const FormatNumericFile defaultFormatNumericFile = txt;
}

/** @enum TypePartition
	@brief type of partition
 */
namespace TypePartition {

	enum TypePartition {

		UNKNOWN_PARTITION = 0,
		label = 1,
		partition = 2
	};
	const TypePartition defaultTypePartition = label;
}

struct NumericPartitionFile {

	std::string _fileName;
	FormatNumeric::FormatNumericFile _format;
	TypePartition::TypePartition _type;
};

enum DataType {

	QualitativeData = 0,
	QuantitativeData,
	HeterogeneousData
};

enum ModelGenre {

	QualitativeModel = 0,
	QuantitativeModel,
	HeterogeneousModel
};

bool isKeyword(std::string& name);

/** @struct TWeightedIndividual
	@brief Structure for chain list of differents sample
 */
struct TWeightedIndividual {

	int64_t val; // index of individual
	double weight;
};

/// XEMCVBlock
struct CVBlock {

	int64_t   _nbSample;                          // number of samples in this CV Block
	double _weightTotal;                          // weight Total of this CV Block
	TWeightedIndividual * _tabWeightedIndividual; // array (size=nbSample) of weighted individual
};

/** @enum ModelName
@brief Enumeration of model name
 */
enum ModelName {

	///////////////////////
	//                   //
	//  Gaussian Models  //
	//                   //
	///////////////////////

	// Unknown model type
	UNKNOWN_MODEL_NAME = -1,

	// 28 Gaussian 'Classical' models

	// Spherical Gaussian model: proportion fixed
	Gaussian_p_L_I = 0,
	Gaussian_p_Lk_I ,

	// Spherical Gaussian model: proportion free
	Gaussian_pk_L_I ,
	Gaussian_pk_Lk_I,

	// Diagonal Gaussian model: proportion fixed
	Gaussian_p_L_B,
	Gaussian_p_Lk_B,
	Gaussian_p_L_Bk,
	Gaussian_p_Lk_Bk,

	// Diagonal Gaussian model: proportion free
	Gaussian_pk_L_B,
	Gaussian_pk_Lk_B,
	Gaussian_pk_L_Bk,
	Gaussian_pk_Lk_Bk,

	// Ellipsoidal Gaussian model: proportion fixed
	Gaussian_p_L_C,
	Gaussian_p_Lk_C,
	Gaussian_p_L_D_Ak_D,
	Gaussian_p_Lk_D_Ak_D,
	Gaussian_p_L_Dk_A_Dk,
	Gaussian_p_Lk_Dk_A_Dk,
	Gaussian_p_L_Ck,
	Gaussian_p_Lk_Ck,

	// Ellipsoidal Gaussian model: proportion free
	Gaussian_pk_L_C,
	Gaussian_pk_Lk_C,
	Gaussian_pk_L_D_Ak_D,
	Gaussian_pk_Lk_D_Ak_D,
	Gaussian_pk_L_Dk_A_Dk,
	Gaussian_pk_Lk_Dk_A_Dk,
	Gaussian_pk_L_Ck,
	Gaussian_pk_Lk_Ck,

	//----------------//
	// 16 HD models   //
	//----------------//
	Gaussian_HD_p_AkjBkQkDk,
	Gaussian_HD_p_AkBkQkDk,
	Gaussian_HD_p_AkjBkQkD ,
	Gaussian_HD_p_AjBkQkD ,
	Gaussian_HD_p_AkjBQkD  ,
	Gaussian_HD_p_AjBQkD  ,
	Gaussian_HD_p_AkBkQkD ,
	Gaussian_HD_p_AkBQkD ,

	Gaussian_HD_pk_AkjBkQkDk,
	Gaussian_HD_pk_AkBkQkDk ,
	Gaussian_HD_pk_AkjBkQkD,
	Gaussian_HD_pk_AjBkQkD,
	Gaussian_HD_pk_AkjBQkD,
	Gaussian_HD_pk_AjBQkD ,
	Gaussian_HD_pk_AkBkQkD,
	Gaussian_HD_pk_AkBQkD ,

	////////////////////////
	//                    //
	//  10 Binary Models  //
	//                    //
	////////////////////////

	// proportion fixed
	Binary_p_E ,
	Binary_p_Ek,
	Binary_p_Ej ,
	Binary_p_Ekj,
	Binary_p_Ekjh ,
	// proportion free
	Binary_pk_E ,
	Binary_pk_Ek,
	Binary_pk_Ej ,
	Binary_pk_Ekj,
	Binary_pk_Ekjh ,

	// Heterogeneous model name:proportions free
	Heterogeneous_pk_E_L_B,
	Heterogeneous_pk_E_Lk_B,
	Heterogeneous_pk_E_L_Bk,
	Heterogeneous_pk_E_Lk_Bk,
	Heterogeneous_pk_Ek_L_B,
	Heterogeneous_pk_Ek_Lk_B,
	Heterogeneous_pk_Ek_L_Bk,
	Heterogeneous_pk_Ek_Lk_Bk,
	Heterogeneous_pk_Ej_L_B,
	Heterogeneous_pk_Ej_Lk_B,
	Heterogeneous_pk_Ej_L_Bk,
	Heterogeneous_pk_Ej_Lk_Bk,
	Heterogeneous_pk_Ekj_L_B,
	Heterogeneous_pk_Ekj_Lk_B,
	Heterogeneous_pk_Ekj_L_Bk,
	Heterogeneous_pk_Ekj_Lk_Bk,
	Heterogeneous_pk_Ekjh_L_B,
	Heterogeneous_pk_Ekjh_Lk_B,
	Heterogeneous_pk_Ekjh_L_Bk,
	Heterogeneous_pk_Ekjh_Lk_Bk,
	// Heterogeneous model name:proportions fix
	Heterogeneous_p_E_L_B,
	Heterogeneous_p_E_Lk_B,
	Heterogeneous_p_E_L_Bk,
	Heterogeneous_p_E_Lk_Bk,
	Heterogeneous_p_Ek_L_B,
	Heterogeneous_p_Ek_Lk_B,
	Heterogeneous_p_Ek_L_Bk,
	Heterogeneous_p_Ek_Lk_Bk,
	Heterogeneous_p_Ej_L_B,
	Heterogeneous_p_Ej_Lk_B,
	Heterogeneous_p_Ej_L_Bk,
	Heterogeneous_p_Ej_Lk_Bk,
	Heterogeneous_p_Ekj_L_B,
	Heterogeneous_p_Ekj_Lk_B,
	Heterogeneous_p_Ekj_L_Bk,
	Heterogeneous_p_Ekj_Lk_Bk,
	Heterogeneous_p_Ekjh_L_B,
	Heterogeneous_p_Ekjh_Lk_B,
	Heterogeneous_p_Ekjh_L_Bk,
	Heterogeneous_p_Ekjh_Lk_Bk,

	nbModelName = 54
};

const ModelName defaultGaussianModelName = Gaussian_pk_Lk_C;
const ModelName defaultBinaryModelName = Binary_pk_Ekjh;
const ModelName defaultGaussianHDDAModelName = Gaussian_HD_pk_AkjBkQkD;
const ModelName defaultHeterogeneousModelName = Heterogeneous_pk_Ekjh_Lk_Bk;

// Output mode
enum OutputType {

	BICstandardOutput = 0,       // Standard output mode
	BICnumericStandardOutput,    // Numerical standard output
	BIClabelOutput,              // Label output
	BICparameterOutput,          // Parameter output (numerical)
	BICtikOutput,                // Posterior probabilities output
	BICzikOutput,                // Partition output (notation 1/0)
	BIClikelihoodOutput,         // Log-likelihood, entropy & completed log-likelihood output
	BICnumericLikelihoodOutput,
	BICErrorOutput = 8,          // error code for BIC Criterion for each estimation

	CVstandardOutput,            // CV Standard output mode
	CVnumericStandardOutput,     // CV Numerical standard output
	CVlabelOutput,               // CV Label output
	CVparameterOutput,           // CV Parameter output (numerical)
	CVtikOutput,                 // CV Posterior probabilities output
	CVzikOutput,                 // CV Partition output (notation 1/0)
	CVlikelihoodOutput,          // CV Log-likelihood, entropy & completed log-likelihood output
	CVnumericLikelihoodOutput,
	CVErrorOutput = 17,          // error code for CV Criterion for each estimation

	ICLstandardOutput,           // ICL Standard output mode
	ICLnumericStandardOutput,    // ICL Numerical standard output
	ICLlabelOutput,              // ICL Label output
	ICLparameterOutput,          // ICL Parameter output (numerical)
	ICLtikOutput,                // ICL Posterior probabilities output
	ICLzikOutput,                // ICL Partition output (notation 1/0)
	ICLlikelihoodOutput,         // ICL Log-likelihood, entropy & completed log-likelihood output
	ICLnumericLikelihoodOutput,
	ICLErrorOutput = 26,         // error code for ICL Criterion for each estimation

	NECstandardOutput,           // NEC Standard output mode
	NECnumericStandardOutput,    // NEC Numerical standard output
	NEClabelOutput,              // NEC Label output
	NECparameterOutput,          // NEC Parameter output (numerical)
	NECtikOutput,                // NEC Posterior probabilities output
	NECzikOutput,                // NEC Partition output (notation 1/0)
	NEClikelihoodOutput,         // NEC Log-likelihood, entropy & completed log-likelihood output
	NECnumericLikelihoodOutput,
	NECErrorOutput = 35,         // error code for NEC Criterion for each estimation

	DCVstandardOutput,           // DCV Standard output mode
	DCVnumericStandardOutput,    // DCV Numerical standard output
	DCVlabelOutput,              // DCV Label output
	DCVparameterOutput,          // DCV Parameter output (numerical)
	DCVtikOutput,                // DCV Posterior probabilities output
	DCVzikOutput,                // DCV Partition output (notation 1/0)
	DCVlikelihoodOutput,         // DCV Log-likelihood, entropy & completed log-likelihood output
	DCVnumericLikelihoodOutput,
	DCVErrorOutput = 44,         // error code for DCV Criterion for each estimation

	completeOutput ,             // Complete output mode
	numericCompleteOutput,       // Numerical complete output

	CVlabelClassificationOutput, // label of classification CV method

	errorMixmodOutput,           // error code for mixmod execution
	errorModelOutput,            // error code for NEC Criterion for each estimation

	DCVinfo,                     // double cross validation information
	DCVnumericInfo = 51          // numeric double cross validation information (for validation test)
};

/// compute a^b and throw an error if it's equal to zero
double powAndCheckIfNotNull(double a, double b,
		const Exception & errorType = NumericException("Defaulter", 0, nullDeterminant) );

/// return the nearest int64_t
int64_t Round(double d);

/// convert big char of a string in low char
void ConvertBigtoLowString(std::string & str);

//ModelNameToString
std::string ModelNameToString(const ModelName & modelName);

//StringToModelName
ModelName StringToModelName(const std::string & strModelName);

//get Heterogeneous model name
ModelName getHeterogeneousModelName(const ModelName binaryName, const ModelName gaussianName);
//get binary model name
ModelName getBinaryModelNamefromHeterogeneous(const ModelName HeterogeneousName);
//get gaussian model name
ModelName getGaussianModelNamefromHeterogeneous(const ModelName HeterogeneousName);

// edit modelName
void edit(const ModelName & modelName);

//criterionNameToString
std::string CriterionNameToString(const CriterionName & criterionName);

//StringtoXEMCriterionName
CriterionName StringtoCriterionName(const std::string & str);

// edit CriterionName
void edit(const CriterionName & criterionName);

// AlgoNameToString
std::string AlgoNameToString(const AlgoName & typeAlgo);

//StringToAlgoName
AlgoName StringToAlgoName(const std::string & str);

// edit AlgoName
void edit(const AlgoName & typeAlgo);

//FormatFileToString
std::string FormatNumericFileToString(const FormatNumeric::FormatNumericFile & formatNumericFile);

//StringToFormatFile
FormatNumeric::FormatNumericFile StringToFormatNumericFile(const std::string & strFormatNumericFile);

//TypePartitionToString
std::string TypePartitionToString(const TypePartition::TypePartition & typePartition);

//StringToTypePartition
TypePartition::TypePartition StringToTypePartition(const std::string & strTypePartition);

//StrategyInitNameToString
std::string StrategyInitNameToString(const StrategyInitName & strategyInitName);
//StrategyInitNameToString for 3rd party env (R, Py, etc.)
std::string StrategyInitNameToStringApp(const StrategyInitName & strategyInitName);

//StringToStrategyInitName
StrategyInitName StringToStrategyInitName(const std::string & str);

// edit StrategyInitName
void edit(const StrategyInitName & strategyInitName);

// AlgoStopNameToString
std::string AlgoStopNameToString(const AlgoStopName & algoStopName);

// void AlgoStopName
void edit(const AlgoStopName & algoStopName);

// is modelName has free proportion
bool hasFreeProportion(ModelName modelName);

// is modelName a diagonal Gaussian Model
bool isDiagonal(ModelName modelName);

// is modelName a spherical Gaussian Model
bool isSpherical(ModelName modelName);

// is modelName a general Gaussian Model
bool isGeneral(ModelName modelName);

// is modelName a EDDA (Classical Gaussian)
bool isEDDA(ModelName modelName);

// is modelName a HD (or HDk)
bool isHD(ModelName modelName);

bool isFreeSubDimension(ModelName modelName);

bool isBinary(ModelName modelName);
bool isHeterogeneous(ModelName modelname);

ModelGenre getModelGenre(ModelName modelname);


// traitement sur les tableaux
//----------------------------

// T* copyTab(T * tab, int64_t dim)
template<typename T> T * copyTab(T * tab, int64_t dim) {
	T * res = new T[dim];
	int64_t i;
	for (i = 0; i < dim; i++) {
		res[i] = tab[i];
	}
	return res;
}

// T ** copyTab(T ** tab, int64_t dim1, int64_t dim2)
template<typename T> T ** copyTab(T ** tab, int64_t dim1, int64_t dim2) {
	T ** res = new T*[dim1];
	int64_t i, j;
	for (i = 0; i < dim1; i++) {
		res[i] = new T[dim2];
		for (j = 0 ; j < dim2; j++)
			res[i][j] = tab[i][j];
	}
	return res;
}

// void recopyTab(T * source, T * destination,int64_t dim)
template<typename T> void recopyTab(T * source, T * destination, int64_t dim) {
	int64_t i;
	for (i = 0; i < dim; i++) {
		destination[i] = source[i];
	}
}

inline void recopyTab(int64_t * source, double * destination, int64_t dim) {
	int64_t i;
	for (i = 0; i < dim; i++) {
		destination[i] = source[i];
	}
}

inline void recopyTabToVector(double ** source,
		std::vector<std::vector<double> > & destination, int64_t dim1, int64_t dim2)
{
	destination.resize(dim1);
	int64_t i, j;
	for (i = 0; i < dim1; i++) {
		destination[i].resize(dim2);
		for (j = 0; j < dim2; j++) {
			destination[i][j] = source[i][j];
		}
	}
}

inline void recopyTabToVector(int64_t * source, std::vector<int64_t> & destination, int64_t dim1) {
	destination.resize(dim1);
	int64_t i;
	for (i = 0; i < dim1; i++) {
		destination[i] = source[i];
	}
}

inline void recopyVectorToTab(std::vector<std::vector<double> > source, double **&  destination) {
	int64_t dim1 = source.size();
	int64_t dim2 = source[0].size();
	destination = new double*[dim1];
	for (int64_t i = 0; i < dim1; i++) {
		destination[i] = new double[dim2];
		for (int64_t k = 0; k < dim2; k++) {
			destination[i][k] = source[i][k];
		}
	}
}

inline void recopyVectorToTab(std::vector<int64_t> source, int64_t *&  destination) {
	int64_t dim1 = source.size();
	destination = new int64_t[dim1];
	for (int64_t i = 0; i < dim1; i++) {
		destination[i] = source[i];
	}
}

// void recopyTab(T ** source, T ** destination, int64_t dim1, int64_t dim2)
template<typename T> void recopyTab(T ** source, T ** destination, int64_t dim1, int64_t dim2) {
	int64_t i, j;
	for (i = 0; i < dim1; i++) {
		for (j = 0; j < dim2; j++) {
			destination[i][j] = source[i][j];
		}
	}
}

void editSimpleTab(double * tab, int64_t n, std::ostream & flux, std::string sep = " ", std::string before = " ");
void editSimpleTab(int64_t    * tab, int64_t n, std::ostream & flux);

template<typename T> void editTab(T ** tab, int64_t dim1, int64_t dim2,
		std::ostream & flux, std::string sep = " ", std::string before = "")
{
	T ** p_tab = tab;
	T *  p_tab_i;
	int64_t i, j ;
	for (i = 0; i < dim1; i++) {
		p_tab_i = *p_tab;
		flux << before;
		for (j = 0; j < dim2; j++)
      putDoubleInStream(flux, p_tab_i[j], sep);
		flux << std::endl;
		p_tab++;
	}
}

//Deleters for unique_ptr
template<typename T> 
 struct TabDeleter {
  TabDeleter(int64_t size) : _size(size){};
  int64_t _size;
  void operator()(T** p) {
    for(int64_t i=0;i<_size;i++){
      delete[] p[i];
    };
    delete[] p;

  }
};

template<typename T> 
 struct VectTabDeleter {
 VectTabDeleter(int64_t sizeV, int64_t sizeT) : _sizeV(sizeV), _sizeT(sizeT) {};
  int64_t _sizeV; 
  int64_t _sizeT;  
  void operator()(T*** p) {
    for(int64_t v=0;v<_sizeV;v++){
      for(int64_t i=0;i<_sizeT;i++){
        delete[] p[v][i];
      };
      delete[] p[v];
    }
    delete[] p;

  }
};
 
template<typename T> 
 struct IfChangedDeleter {
  IfChangedDeleter(T* initial) : _initial(initial){};
  T* _initial;
  void operator()(T* p) {
    if(p != _initial) delete p;
  }
};


 
// move on a file until *what* is reached
void moveUntilReach(std::ifstream & fi, std::string what = "datafile");

void readTabFileName(std::ifstream & fi, int64_t nbNbCluster, std::string * tabFileName, std::string & keyWord);

void initToZero(double * tab, int64_t n);
void initToZero(double * tab, int64_t n);

const int64_t SMALL_ENOUGH_TO_USE_SELECTION_SORT = 15;
void echange(double * tab, int64_t i1, int64_t i2);
void echange(int64_t * tab  , int64_t i1, int64_t i2);

void selectionSortWithOrder(double * tabRandom, int64_t * tabOrder, int64_t left, int64_t right);

int64_t partition(double * tabRandom, int64_t * tabOrder, int64_t left, int64_t right);

void quickSortWithOrder(double * tabRandom, int64_t * tabOrder, int64_t left, int64_t right);

int64_t generateRandomIndex(bool * tabIndividualCanBeUsedForInitRandom, double * weight, double totalWeight);

void inputCriterion(std::ifstream & fi, CriterionName & criterionName);

void inputCVinitBlocks(std::ifstream & fi, CVinitBlocks cVinitBlocks);

void inputDCVinitBlocks(std::ifstream & fi, DCVinitBlocks dCVinitBlocks);

 
}

#endif
