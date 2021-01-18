/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/InputException.h  description
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
/** @file XEMInputException.h
 *  @brief Exception class for Input error handling.
 *  @author Parmeet Bhatia
 **/

#ifndef XEM_INPUTEXCEPTION_H
#define XEM_INPUTEXCEPTION_H

#include <string>
#include "mixmod/Utilities/exceptions/Exception.h"

namespace XEM {

class InputException : public Exception {

public:

	Exception * clone() throw ();
	InputException(std::string, int, InputError) throw ();
	InputException(InputError) throw ();
	InputException(const InputException & inputException);
	virtual const char* what() const throw ();
	virtual bool operator==(const Exception&) const throw ();
	virtual void run(std::ostream & flux = std::cout) const throw ();

	virtual ~InputException() throw () {
	}

	static std::map<InputError, const char*> create_map()
	{
		std::map<InputError, const char*> m;

		m.insert(std::make_pair(noError, "No error"));
		m.insert(std::make_pair(nbLinesTooLarge, "Number of lines too large"));
		m.insert(std::make_pair(nbLinesTooSmall, "Number of lines too small"));
		m.insert(std::make_pair(pbDimensionTooLarge, "Dimension size too large"));
		m.insert(std::make_pair(pbDimensionTooSmall, "Dimension size too small"));
		m.insert(std::make_pair(nbCriterionTooLarge, "Number of criterion too large"));
		m.insert(std::make_pair(nbCriterionTooSmall, "Number of criterion too small"));
		m.insert(std::make_pair(wrongCriterionName, "Wrong criterion name/type"));
		m.insert(std::make_pair(nbNbClusterTooLarge, "List of number of clusters too large"));
		m.insert(std::make_pair(nbNbClusterTooSmall, "List of number of clusters too small"));
		m.insert(std::make_pair(nbModelTypeTooLarge, "Number of models too large"));
		m.insert(std::make_pair(nbModelTypeTooSmall, "Number of models too small"));
		m.insert(std::make_pair(wrongModelType, "Wrong model name/type"));
		m.insert(std::make_pair(wrongCVinitType, "Wrong CVinitType"));
		m.insert(std::make_pair(wrongDCVinitType, "Wrong DCVinitType"));
		m.insert(std::make_pair(nbStrategyTypeTooLarge, "Number of strategies too large"));
		m.insert(std::make_pair(nbStrategyTypeTooSmall, "Number of strategies too small"));
		m.insert(std::make_pair(wrongStrategyInitName, "Wrong strategy initialization name"));
		m.insert(std::make_pair(badNbParameterInInit, "Bad parameter number in Init"));
		m.insert(std::make_pair(badNbPartitionInInit, "Bad partition number in Init"));
		m.insert(std::make_pair(errorInitParameter, "Error with USER initialization"));
		m.insert(std::make_pair(nbAlgoTooLarge, "Number of algorithms too large"));
		m.insert(std::make_pair(nbAlgoTooSmall, "Number of algorithms too small"));
		m.insert(std::make_pair(wrongAlgoType, "Wrong algorithm type"));
		m.insert(std::make_pair(nbIterationTooLarge, "Number of iterations too large"));
		m.insert(std::make_pair(nbIterationTooSmall, "Number of iterations too small"));
		m.insert(std::make_pair(epsilonTooSmall, "Value of epsilon too small"));
		m.insert(std::make_pair(epsilonTooLarge, "Value of epsilon too large"));
		m.insert(std::make_pair(wrongDataFileName, "Wrong data file name"));
		m.insert(std::make_pair(wrongLabelFileName, "Wrong label file name"));
		m.insert(std::make_pair(wrongWeightFileName, "Wrong weight file name"));
		m.insert(std::make_pair(wrongParamFileName, "Wrong parameter file name"));
		m.insert(std::make_pair(wrongPartitionFileName, "Wrong label file name"));
		m.insert(std::make_pair(wrongAlgoStopName, "Wrong stopping rules for algorithm"));
		m.insert(std::make_pair(wrongOutputType, "Wrong output mode type"));
		m.insert(std::make_pair(wrongInputFileName, "Wrong input file name"));
		m.insert(std::make_pair(wrongXEMNbParam, "Wrong number of paramaters for Mixmod call"));
		m.insert(std::make_pair(errorNbLines, "Bad writing \"NbLines\" key word"));
		m.insert(std::make_pair(errorPbDimension, "Bad writing \"PbDimension\" key word"));
		m.insert(std::make_pair(errorNbCriterion, "Bad writing \"NbCriterion\" key word"));
		m.insert(std::make_pair(errorListCriterion, "Bad writing \"ListCriterion\" key word"));
		m.insert(std::make_pair(errorNbNbCluster, "Bad writing \"NbNbCluster\" key word"));
		m.insert(std::make_pair(errorListNbCluster, "Bad writing \"ListNbCluster\" key word"));
		m.insert(std::make_pair(errorNbModel, "Bad writing \"NbModel\" key word"));
		m.insert(std::make_pair(errorListModel, "Bad writing \"ListModel\" key word"));
		m.insert(std::make_pair(errorNbStrategy, "Bad writing \"NbStrategy\" key word"));
		m.insert(std::make_pair(errorInitType, "Bad writing \"InitType\" key word"));
		m.insert(std::make_pair(errorInitFile, "Bad writing \"InitFile\" key word"));
		m.insert(std::make_pair(errorNbAlgo, "Bad writing \"NbAlgorithm\" key word"));
		m.insert(std::make_pair(errorAlgo, "Bad writing \"Algorithm\" key word"));
		m.insert(std::make_pair(errorStopRule, "Bad writing \"StopRule\" key word"));
		m.insert(std::make_pair(errorStopRuleValue, "Bad writing \"StopRuleValue\" key word"));
		m.insert(std::make_pair(errorDataFile, "Bad writing \"DataFile\" key word"));
		m.insert(std::make_pair(nbAlgoTypeTooSmall, "number of algoType too small"));
		m.insert(std::make_pair(badStrategyInitName, "strategyInitName incompatible with algoType"));
		m.insert(std::make_pair(errorOpenFile, "Error when opening a output file"));
		m.insert(std::make_pair(errorNbModality, "Minimum number of modality is 2"));
		m.insert(std::make_pair(knownPartitionNeeded, "known partition is needed for M algorithm"));
		m.insert(std::make_pair(badKnownPartition, "bad known partition"));
		m.insert(std::make_pair(endDataFileReach, "the end of data file has been reached before reading all samples] verify nbSample or data file"));
		m.insert(std::make_pair(wrongNbKnownPartition, "Error : wrong number of known Partition in input object"));
		m.insert(std::make_pair(SubDimensionFreeTooLarge, "SubDimensionFree is too large"));
		m.insert(std::make_pair(SubDimensionFreeTooSmall, "SubDimensionFree is too small"));
		m.insert(std::make_pair(SubDimensionEqualTooLarge, "SubDimensionEqual is too large"));
		m.insert(std::make_pair(SubDimensionEqualTooSmall, "SubDimensionEqual is too small"));
		m.insert(std::make_pair(weightTotalIsNotAnInteger, "Error : weightTotal must be an integer"));
		m.insert(std::make_pair(ungivenSubDimension, "Error : sub dimensions are not given for one or several models"));
		m.insert(std::make_pair(wrongNbAlgoWhenMorMAP, "Error : wrong number of algortihms if M or MAP are used"));
		m.insert(std::make_pair(BadInitialsationWhenM, "Error : USER_PARTITION must be the initialisation if M is used"));
		m.insert(std::make_pair(BadInitialsationWhenMAP, "Error : USER must be the initialisation if MAP is used"));
		m.insert(std::make_pair(partitionMustBeComplete, "Error : partition must be complete"));
		m.insert(std::make_pair(inputNotFinalized, "Error : input is not finalized"));
		m.insert(std::make_pair(algorithmMustBeM, "Error : algorithm must be M"));
		m.insert(std::make_pair(knownPartitionAndInitPartitionMustBeEqual, "Error : knownLabel And InitLabel must be equal"));
		m.insert(std::make_pair(nbStrategyMustBe1, "Error : nbStrategy must be equal to 1"));
		m.insert(std::make_pair(wrongNbKnownPartitionOrInitPartition, "Error : wrong number of knownLabel or InitLabel"));
		m.insert(std::make_pair(tooManySampleInInitPartitionAndTooManyClusterNotRepresented, "Error : error in USER_PARTITION initialization : Too many sample in InitPartition and too many cluster not repented"));
		m.insert(std::make_pair(notAvailableForPrediction, "Not available for prediction"));
		m.insert(std::make_pair(wrongValueInMultinomialCase, "wrong value in data set : use 1,2...nbModality"));
		m.insert(std::make_pair(errorInPartitionInput, "Error in partition file : there is not enough lines in the file (nbSample is required)"));
		m.insert(std::make_pair(badNumberOfValuesInLabelInput, "Error in label file : the number of values bust be nbSample"));
		m.insert(std::make_pair(notEnoughValuesInProbaInput, "Error in proba file : there is not enough values in the file (nbSample*nbCluster is required)"));
		m.insert(std::make_pair(badValueInLabelInput, "Error in label file : label must be between 1 and nbCluster"));
		m.insert(std::make_pair(ColumnTypeNotValid, "Bad Format"));
		m.insert(std::make_pair(badStopNameWithSEMAlgo, "Error : bad stop type with SEM : this algortihm must be stopped after a predefined number of iterations"));
		m.insert(std::make_pair(badAlgorithmInHDContext, "Error : bad algorithm in HD context : only M or MAP is available"));
		m.insert(std::make_pair(differentSubDimensionsWithMAP, "Error : given subDimensions in init file and input file are different"));
		m.insert(std::make_pair(wrongSubDimension, "Error : Wrong sub dimension type for given model"));
		m.insert(std::make_pair(missingRequiredInputs, "Error : Missing required inputs (data, nbSample, pbDimension, tabNbCluster, nbNbCluster)"));
		m.insert(std::make_pair(wrongCriterionPositionInSet, "Wrong criterion position in set"));
		m.insert(std::make_pair(wrongCriterionPositionInGet, "Wrong criterion position in get"));
		m.insert(std::make_pair(wrongCriterionPositionInInsert, "Wrong criterion position in insert"));
		m.insert(std::make_pair(wrongCriterionPositionInRemove, "Wrong criterion position in remove"));
		m.insert(std::make_pair(wrongModelPositionInSet, "Wrong model position in set"));
		m.insert(std::make_pair(wrongModelPositionInGet, "Wrong model position in get"));
		m.insert(std::make_pair(wrongModelPositionInInsert, "Wrong model position in insert"));
		m.insert(std::make_pair(wrongModelPositionInRemove, "Wrong model position in remove"));
		m.insert(std::make_pair(wrongModelPositionInSetSubDimensionEqual, "Wrong model position in set sub dimension equal"));
		m.insert(std::make_pair(wrongModelPositionInSetSubDimensionFree, "Wrong model position in set sub dimension free"));
		m.insert(std::make_pair(wrongModelInSetSubDimensionEqual, "SetSubDimensionEqual could not be called with this model"));
		m.insert(std::make_pair(wrongModelInSetSubDimensionFree, "SetSubDimensionFree could not be called with this model"));
		m.insert(std::make_pair(badSetKnownPartition, "Error in setKnownPartition (impossible if nbNbCluster>1)"));
		m.insert(std::make_pair(wrongStrategyPositionInSetOrGetMethod, "Wrong strategy position in set or get method"));
		m.insert(std::make_pair(badInitPart, "Bad Initialization Partition"));
		m.insert(std::make_pair(nbTryInStrategyTooSmall, "Number of tries in strategy too small"));
		m.insert(std::make_pair(nbTryInStrategyTooLarge, "Number of tries in strategy too large"));
		m.insert(std::make_pair(nbTryInInitTooSmall, "Number of tries in init too small"));
		m.insert(std::make_pair(nbTryInInitTooLarge, "Number of tries in init too large"));
		m.insert(std::make_pair(nbIterationInInitTooSmall, "Number of iterations in init too small"));
		m.insert(std::make_pair(nbIterationInInitTooLarge, "Number of iterations in init too large"));
		m.insert(std::make_pair(epsilonInInitTooSmall, "Epsilon in init too small"));
		m.insert(std::make_pair(epsilonInInitTooLarge, "Epsilon in init too large"));
		m.insert(std::make_pair(wrongNbStrategyTryValue, "Wrong number of tries in strategy"));
		m.insert(std::make_pair(badSetNbTry, "Number of tries in strategy could not change"));
		m.insert(std::make_pair(badSetNbTryInInit, "Number of tries in init could not change"));
		m.insert(std::make_pair(badSetNbIterationInInit, "Number of iterations in init could not change"));
		m.insert(std::make_pair(badSetEpsilonInInit, "Epsilon in init could not change"));
		m.insert(std::make_pair(badSetStopNameInInit, "Stop name could not change in this context"));
		m.insert(std::make_pair(badCriterion, "Bad Criterion"));
		m.insert(std::make_pair(badAlgo, "Bad Algorithm"));
		m.insert(std::make_pair(badAlgoStop, "Bad Algorithm Stop Name"));
		m.insert(std::make_pair(DAInput, "XEMDAInput not implemented"));
		m.insert(std::make_pair(wrongModelName, "Wrong model Name"));
		m.insert(std::make_pair(knownPartitionNotAvailable, "known Partition is not available"));
		m.insert(std::make_pair(tooManyWeightColumnDescription, "Too many WeightColumnDescription"));
		m.insert(std::make_pair(badDataDescription, "Bad Data Description"));
		m.insert(std::make_pair(badLabelDescription, "Bad Label Description"));
		m.insert(std::make_pair(errorInColumnDescription, "Bad size of Column Description"));
		m.insert(std::make_pair(errorInXEMInputSelector, "Bad size of Column Description"));
		m.insert(std::make_pair(wrongIndexInGetMethod, "wrong index in get method"));
		m.insert(std::make_pair(badInputType, "Bad Input type"));
		m.insert(std::make_pair(wrongKnownPartitionPositionInSet, "wrong known partition position in set"));
		m.insert(std::make_pair(wrongKnownPartitionPositionInRemove, "wrong known partition position in remove"));
		m.insert(std::make_pair(HDModelsAreNotAvailableInClusteringContext, "HD Models Are Not Available In Clustering Context"));

		return m;
	}

	static std::map<InputError, const char*> mapErrorMsg;

protected:

	InputError _errorType;
};

}

#endif /* XEMINPUTEXCEPTION_H_ */
