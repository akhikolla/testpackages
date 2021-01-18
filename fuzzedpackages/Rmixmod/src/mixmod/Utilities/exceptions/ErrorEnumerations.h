/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/ErrorEnumerations.h  description
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
/** @file ErrorEnumerations.h
 *  @brief Enumerations for various types of errors.
 *  @author Parmeet Bhatia
 **/

#ifndef ERRORENUMERATIONS_H
#define ERRORENUMERATIONS_H

namespace XEM {

enum InputError {

	noError = 0,
	nbLinesTooLarge,
	nbLinesTooSmall,
	pbDimensionTooLarge,
	pbDimensionTooSmall,
	nbCriterionTooLarge,
	nbCriterionTooSmall,
	wrongCriterionName,
	nbNbClusterTooLarge,
	nbNbClusterTooSmall,
	nbModelTypeTooLarge,
	nbModelTypeTooSmall,
	wrongModelType,
	wrongCVinitType,
	wrongDCVinitType,
	nbStrategyTypeTooLarge,
	badNbParameterInInit,
	badNbPartitionInInit,
	nbStrategyTypeTooSmall,
	wrongStrategyInitName,
	errorInitParameter,
	nbAlgoTooLarge,
	nbAlgoTooSmall,
	wrongAlgoType,
	nbIterationTooLarge,
	nbIterationTooSmall,
	epsilonTooSmall,
	epsilonTooLarge,
	wrongDataFileName,
	wrongWeightFileName,
	wrongParamFileName,
	wrongLabelFileName,
	wrongPartitionFileName,
	wrongAlgoStopName,
	wrongOutputType,
	wrongInputFileName,
	wrongXEMNbParam,
	errorNbLines,
	errorPbDimension,
	errorNbCriterion,
	errorListCriterion,
	errorNbNbCluster,
	errorListNbCluster,
	errorNbModel,
	errorListModel,
	errorNbStrategy,
	errorInitType,
	errorInitFile,
	errorNbAlgo,
	errorAlgo,
	errorStopRule,
	errorStopRuleValue,
	errorDataFile,
	nbAlgoTypeTooSmall,
	badStrategyInitName,
	errorOpenFile,
	errorNbModality,
	knownPartitionNeeded,
	badKnownPartition,
	endDataFileReach,
	wrongNbKnownPartition,
	SubDimensionFreeTooLarge,
	SubDimensionFreeTooSmall,
	SubDimensionEqualTooLarge,
	SubDimensionEqualTooSmall,
	weightTotalIsNotAnInteger,
	ungivenSubDimension,
	inputNotFinalized,
	wrongNbAlgoWhenMorMAP,
	BadInitialsationWhenM,
	BadInitialsationWhenMAP,
	partitionMustBeComplete,
	algorithmMustBeM,
	knownPartitionAndInitPartitionMustBeEqual,
	nbStrategyMustBe1,
	wrongNbKnownPartitionOrInitPartition,
	tooManySampleInInitPartitionAndTooManyClusterNotRepresented,
	notAvailableForPrediction,
	ColumnTypeNotValid,
	errorInPartitionInput,
	wrongValueInMultinomialCase,
	badNumberOfValuesInLabelInput,
	notEnoughValuesInProbaInput,
	badValueInLabelInput,

	badStopNameWithSEMAlgo,
	badAlgorithmInHDContext,
	differentSubDimensionsWithMAP,
	wrongSubDimension,
	missingRequiredInputs,
	// the following errors should not be occur because they must have been see by MVCControler
	wrongCriterionPositionInSet,
	wrongCriterionPositionInGet,
	wrongCriterionPositionInInsert,
	wrongCriterionPositionInRemove,
	wrongModelPositionInSet,
	wrongModelPositionInGet,
	wrongModelPositionInInsert,
	wrongModelPositionInRemove,
	wrongModelPositionInSetSubDimensionEqual,
	wrongModelPositionInSetSubDimensionFree,
	wrongModelInSetSubDimensionEqual,
	wrongModelInSetSubDimensionFree,
	wrongKnownPartitionPositionInSet,
	wrongKnownPartitionPositionInRemove,
	badSetKnownPartition,
	wrongStrategyPositionInSetOrGetMethod,
	nbTryInStrategyTooSmall,
	nbTryInStrategyTooLarge,
	nbTryInInitTooSmall,
	nbTryInInitTooLarge,
	epsilonInInitTooSmall,
	epsilonInInitTooLarge,
	nbIterationInInitTooSmall,
	nbIterationInInitTooLarge,
	wrongNbStrategyTryValue,
	badInitPart,
	badSetNbTry,
	badSetNbTryInInit,
	badSetNbIterationInInit,
	badSetEpsilonInInit,

	badSetStopNameInInit,
	badCriterion,
	badAlgo,
	badAlgoStop,
	badInputType,
	DAInput,
	wrongModelName,
	knownPartitionNotAvailable,
	tooManyWeightColumnDescription,
	badDataDescription,
	badLabelDescription,
	errorInColumnDescription,
	errorInXEMInputSelector,
	wrongIndexInGetMethod,
	HDModelsAreNotAvailableInClusteringContext
};

enum DCVError {

	wrongDCVinitBlocks,
	wrongDCVnumberOfBlocks,
	DCVmustBeDIAG,
	forbiddenCallToGetBestCVModel
};

enum DCVonlyInGaussianCaseError {

	allCVCriterionErrorForAnEstimationInDCVContext,
	NbDCVBlocksTooSmall
};

enum NumericError {

	int64_t_max_error,
	CEM_INIT_error,
	SEM_MAX_error,
	SMALL_EM_error,
	tabNkNotInteger,
	sumFiNullAndfkTPrimNull,
	sumFiNullInMultinomialCase,
	nonPositiveDefiniteMatrix,
	nullDeterminant,
	randomProblem,
	nullLikelihood,
	noProbability,
	pbNEC,
	nullNk,

	numericError,
	errorSigmaConditionNumber,
	minDeterminantSigmaValueError,
	minDeterminantWValueError,
	minDeterminantDiagWkValueError,

	minDeterminantDiagWValueError,
	minDeterminantBValueError,
	minDeterminantRValueError,
	minDeterminantWkValueError,
	minDeterminantShapeValueError,
	minDeterminantDiagQtmpValueError
};

enum OtherError {

	badFormat,
	nullPointerError,
	wrongMatrixType,
	wrongConstructorType,
	nonImplementedMethod,
	badBinaryParameterClass,
	internalMixmodError,
	FunctionNotYetImplemented,
	AllTriesGotErros,
	AllModelsGotErros,
    xmlFeaturesNotAvailable,
	UnknownReason
};

}

#endif /* ERRORENUMERATIONS_H_ */
