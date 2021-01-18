#ifndef STANDALONE
#include "include/base/Parameter.h"
#include "include/ROC/ROCParameter.h"
#include "include/PA/PAParameter.h"
#include "include/PANSE/PANSEParameter.h"
#include "include/FONSE/FONSEParameter.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(Trace)
RCPP_EXPOSED_CLASS(Genome)
RCPP_EXPOSED_CLASS(CovarianceMatrix)

//' @name readPhiValue
//' @title readPhiValue
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Read synthesis rate values from file. File should be two column file <gene_id,phi> and is expected to have a header row
//' @param filename name of file to be read


//' @name setGroupList
//' @title setGroupList
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Set amino acids (ROC, FONSE) or codons (PA, PANSE) for which parameters will be estimated. Note that non-default groupLists are still in beta testing and should be used with caution.
//' @param List of strings epresenting groups for parameters to be estimated. Should be one letter amino acid (ROC, FONSE) or list of sense codons (PA, PANSE). 

//' @name getGroupList
//' @title getGroupList
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get amino acids (ROC, FONSE) or codons (PA, PANSE) for which parameters will be estimated
//' @return returns list of amino acids or codons


//' @name getTraceObject
//' @title getTraceObject
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get Trace object stored by a Parameter object. Useful for plotting certain parameter traces.
//' @return Trace object

//' @name initializeSynthesisRateByGenome
//' @title initializeSynthesisRateByGenome
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Initialize synthesis rates using SCUO values calcuated from the genome
//' @param genome a Genome object

//' @name initializeSynthesisRateByRandom
//' @title initializeSynthesisRateByRandom
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Initialize synthesis rates by drawing a from a lognormal distribution with mean = -(sd_phi)^2/2 and sd = sd_phi
//' @param sd_phi a positive value which will be the standard deviation of the lognormal distribution

//' @name initializeSynthesisRateByList
//' @title initializeSynthesisRateByList
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Initialize synthesis rates with values passed in as a list
//' @param expression a list of values to use as initial synthesis rate values. Should be same size as number of genes in genome.

//' @name getSynthesisRate
//' @title getSynthesisRate
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get current synthesis rates for all genes and all mixtures 
//' @return 2 by 2 vector of numeric values

//' @name getSynthesisRatePosteriorMeanForGene
//' @title getSynthesisRatePosteriorMeanForGene
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get posterior mean synthesis rate value for a gene
//' @param samples number of samples over which to calculate mean
//' @param geneIndex corresponding index of gene in genome for which posterior mean synthesis rate will be calculated. Should be a number between 1 and length(genome) 
//' @param log_scale Calculate posterior mean on log scale
//' @return posterior mean synthesis rate for gene

//' @name getSynthesisRatePosteriorVarianceForGene
//' @title getSynthesisRatePosteriorVarianceForGene
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get synthesis rate variance for a gene
//' @param samples number of samples over which to calculate variance
//' @param geneIndex corresponding index of gene in genome for which synthesis rate variance will be calculated. Should be a number between 1 and length(genome) 
//' @param unbiased Should calculate variance using unbiased (N-1) or biased (N) correction
//' @param log_scale Calculate variance on log scale
//' @return posterior mean synthesis rate for gene

//' @name getEstimatedMixtureAssignmentForGene
//' @title getEstimatedMixtureAssignmentForGene
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get estimated mixture assignment for gene
//' @param samples number of samples over which to calculate mixture assignment
//' @param geneIndex corresponding index of gene in genome. Should be a number between 1 and length(genome). 
//' @return returns value between 1 and n, where n is number of mixtures

//' @name getEstimatedMixtureAssignmentProbabilitiesForGene
//' @title getEstimatedMixtureAssignmentProbabilitiesForGene
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get estimated mixture assignment probabilities for gene
//' @param samples number of samples over which to calculate mixture assignment probabilities
//' @param geneIndex corresponding index of gene in genome. Should be a number between 1 and length(genome). 
//' @return returns vector of probabilities representing mixture probabilities for gene

//' @name getStdDevSynthesisRatePosteriorMean
//' @title getStdDevSynthesisRatePosteriorMean
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate posterior mean of standard deviation parameter of lognormal describing distribution of synthesis rates
//' @param samples number of samples over which to calculate posterior mean
//' @param mixture mixture index to use. Should be number between 0 and n-1, where n is number of mixtures
//' @return returns posterior mean for standard deviation of lognormal distribution of synthesis rates

//' @name getStdDevSynthesisRateVariance
//' @title getStdDevSynthesisRateVariance
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate variance of standard deviation parameter of lognormal describing distribution of synthesis rates
//' @param samples number of samples over which to calculate variance
//' @param mixture mixture index to use. Should be number between 0 and n-1, where n is number of mixtures
//' @param unbiased If TRUE, should calculate variance using unbiased (N-1). Otherwise, used biased (N) correction
//' @return returns variance for standard deviation of lognormal distribution of synthesis rates

//' @name getCodonSpecificPosteriorMeanForCodon
//' @title getCodonSpecificPosteriorMeanForCodon
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate codon-specific parameter (CSP) posterior mean
//' @param mixtureElement mixture to calculate CSP posterior mean. Should be between 1 and n, where n is number of mixtures.
//' @param samples number of samples to use for calculating posterior mean
//' @param codon codon to calculate CSP
//' @param paramType CSP to calculate posterior mean for. 0: Mutation (ROC,FONSE) or Alpha (PA, PANSE). 1: Selection (ROC,FONSE), Lambda (PANSE), Lambda^prime (PA). 2: NSERate (PANSE) 
//' @param withoutReference If model uses reference codon, then ignore this codon (fixed at 0). Should be TRUE for ROC and FONSE. Should be FALSE for PA and PANSE.
//' @param log_scale If true, calculate posterior mean on log scale. Should only be used for PA and PANSE.
//' @return posterior mean value for CSP

//' @name getCodonSpecificPosteriorVarianceForCodon
//' @title getCodonSpecificPosteriorVarianceForCodon
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate codon-specific parameter (CSP) variance
//' @param mixtureElement mixture to calculate CSP variance. Should be between 1 and n, where n is number of mixtures.
//' @param samples number of samples to use for calculating variance
//' @param codon codon to calculate CSP
//' @param paramType CSP to calculate variance for. 0: Mutation (ROC,FONSE) or Alpha (PA, PANSE). 1: Selection (ROC,FONSE), Lambda (PANSE), Lambda^prime (PA). 2: NSERate (PANSE) 
//' @param unbiased If TRUE, should calculate variance using unbiased (N-1). Otherwise, used biased (N) correction
//' @param withoutReference If model uses reference codon, then ignore this codon (fixed at 0). Should be TRUE for ROC and FONSE. Should be FALSE for PA and PANSE.
//' @param log_scale If true, calculate posterior mean on log scale. Should only be used for PA and PANSE.
//' @return variance over trace for CSP

//' @name getCodonSpecificQuantilesForCodon
//' @title getCodonSpecificQuantilesForCodon
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate quantiles of CSP traces
//' @param mixtureElement mixture to calculate CSP variance. Should be between 1 and n, where n is number of mixtures.
//' @param samples number of samples to use for calculating variance
//' @param codon codon to calculate CSP
//' @param paramType CSP to calculate variance for. 0: Mutation (ROC,FONSE) or Alpha (PA, PANSE). 1: Selection (ROC,FONSE), Lambda (PANSE), Lambda^prime (PA). 2: NSERate (PANSE) 
//' @param probs vector of two doubles between 0 and 1 indicating range over which to calculate quantiles. <0.0275, 0.975> would give 95\% quantiles.
//' @param withoutReference If model uses reference codon, then ignore this codon (fixed at 0). Should be TRUE for ROC and FONSE. Should be FALSE for PA and PANSE.
//' @param log_scale If true, calculate posterior mean on log scale. Should only be used for PA and PANSE.
//' @return vector representing lower and upper bound of quantile

//' @name getNoiseOffsetPosteriorMean
//' @title getNoiseOffsetPosteriorMean
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate posterior mean of standard deviation parameter of lognormal describing distribution of synthesis rates
//' @param index mixture index to use. Should be number between 0 and n-1, where n is number of mixtures
//' @param samples number of samples over which to calculate posterior mean
//' @return returns posterior mean for standard deviation of lognormal distribution of synthesis rates

//' @name getNoiseOffsetVariance
//' @title getNoiseOffsetVariance
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Calculate variance of noise offset parameter used when fitting model with empirical estimates of synthesis rates (ie. withPhi fits)
//' @param index mixture index to use. Should be number between 0 and n-1, where n is number of mixtures
//' @param samples number of samples over which to calculate variance
//' @param unbiased If TRUE, should calculate variance using unbiased (N-1). Otherwise, used biased (N) correction
//' @return returns variance for noise offset


//' @name fixSphi
//' @title fixSphi
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Fix the value of s_phi (standard deviation of lognormal for synthesis rates) at its current value


//' @name initMutationCategories
//' @title initMutationCategories
//' @description Initialize values for mutation CSP. File should be of comma-separated with header. Three columns should be of order Amino_acid,Codon,Value
//' @param files list of files containing starting values. Number of files should equal the number of categories.
//' @param numCategories number of mutation categories (should be less than or equal to number of mixtures)
//' @param fix Can use this parameter to fix mutation at current values (won't change over course of MCMC run)

//' @name initSelectionCategories
//' @title initSelectionCategories
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Initialize values for selection CSP. File should be of comma-separated with header. Three columns should be of order Amino_acid,Codon,Value
//' @param files list of files containing starting values. Number of files should equal the number of categories.
//' @param numCategories number of mutation categories (should be less than or equal to number of mixtures)
//' @param fix Can use this parameter to fix selection at current values (won't change over course of MCMC run)

//' @name fixDM
//' @title fixDM
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Fix the value of mutation its current value

//' @name fixDEta
//' @title fixDEta
//' @description Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Fix the value of selection its current value


RCPP_MODULE(Parameter_mod)
{
	class_<Parameter>("Parameter")

		.constructor()
		//Initialization and Restart Functions:
		.method("initializeSynthesisRateByGenome", &Parameter::initializeSynthesisRateByGenome)
		.method("initializeSynthesisRateByList", &Parameter::initializeSynthesisRateByList)
		.method("initializeSynthesisRateByRandom", &Parameter::initializeSynthesisRateByRandom)
		//checkIndex is not listed/exposed since it is only called from the other R functions
		.method("readPhiValues", &Parameter::readPhiValues) //Not a R wrapper



		//Mixture Definition Matrix and Category Functions:
		.method("getMutationCategoryForMixture", &Parameter::getMutationCategoryForMixture)
		.method("getSelectionCategoryForMixture", &Parameter::getSelectionCategoryForMixture)
		.method("getSynthesisRateCategoryForMixture", &Parameter::getSynthesisRateCategoryForMixture)
		.method("getCategories", &Parameter::getCategories)
        .method("setCategories", &Parameter::setCategories)
        .method("setCategoriesForTrace", &Parameter::setCategoriesForTrace)
		//setNumMutationCategories and setNumSelectionCategories are taken care of in the
		//properties section below.



		//Group List Functions:
		.method("getGroupList", &Parameter::getGroupList)
		.method("setGroupList", &Parameter::setGroupList)

		//Trace Functions:
		.method("getTraceObject", &Parameter::getTraceObject) //TODO: only used in R?
		.method("setTraceObject", &Parameter::setTraceObject)

		//Synthesis Rate Functions:
		.method("getSynthesisRate", &Parameter::getSynthesisRateR)
		.method("getCurrentSynthesisRateForMixture", &Parameter::getCurrentSynthesisRateForMixture)

		//Iteration Functions:
		.method("getLastIteration", &Parameter::getLastIteration) //Not a R wrapper
		.method("setLastIteration", &Parameter::setLastIteration) //Not a R wrapper

		//Posterior, Variance, and Estimates Functions:
		.method("getSynthesisRatePosteriorMeanForGene",
		        &Parameter::getSynthesisRatePosteriorMeanForGene)
		.method("getSynthesisRateVarianceForGene",
		        &Parameter::getSynthesisRateVarianceForGene)
		.method("getEstimatedMixtureAssignmentForGene", &Parameter::getEstimatedMixtureAssignmentForGene,
		        "returns the mixture assignment for a given gene") //Not a R wrapper
		.method("getEstimatedMixtureAssignmentProbabilitiesForGene", &Parameter::getEstimatedMixtureAssignmentProbabilitiesForGene,
		        "returns the probabilities assignment for a given gene") //Not a R wrapper
		.method("getStdDevSynthesisRatePosteriorMean", &Parameter::getStdDevSynthesisRatePosteriorMean) //Not a R wrapper
		.method("getCodonSpecificPosteriorMean", &Parameter::getCodonSpecificPosteriorMeanForCodon)
		.method("getStdDevSynthesisRateVariance", &Parameter::getStdDevSynthesisRateVariance)
		.method("getCodonSpecificVariance", &Parameter::getCodonSpecificVarianceForCodon)
        .method("getCodonSpecificQuantile", &Parameter::getCodonSpecificQuantileForCodon)
        .method("getExpressionQuantile", &Parameter::getExpressionQuantileForGene)

        //TODO: these functions are not wrapped! May not run correctly
		.method("getNoiseOffsetPosteriorMean", &Parameter::getNoiseOffsetPosteriorMean)
		.method("getNoiseOffsetVariance", &Parameter::getNoiseOffsetVariance)

		// Noise
		.method("setInitialValuesForSepsilon", &Parameter::setInitialValuesForSepsilon)//Not a R wrapper


		.method("setNumObservedSynthesisRateSets", &Parameter::setNumObservedPhiSets)//Not a R wrapper

		//Other Functions:
		.method("getMixtureAssignment", &Parameter::getMixtureAssignmentR)
		.method("setMixtureAssignment", &Parameter::setMixtureAssignmentR)
		.method("getMixtureAssignmentForGene", &Parameter::getMixtureAssignmentForGeneR)
		.method("setMixtureAssignmentForGene", &Parameter::setMixtureAssignmentForGene)
		//setNumMixtureElements it taken care in the properties section below

		//Other Functions:
		.method("calculateSelectionCoefficients", &Parameter::calculateSelectionCoefficientsR)


		.method("fixSphi",&Parameter::fixStdDevSynthesis)



		//Used for getters and setters
		.property("numMutationCategories", &Parameter::getNumMutationCategories, &Parameter::setNumMutationCategories)
		.property("numSelectionCategories", &Parameter::getNumSelectionCategories, &Parameter::setNumSelectionCategories)
		.property("numMixtures", &Parameter::getNumMixtureElements, &Parameter::setNumMixtureElements)
;



		//only temporary testing
		function("randNorm", &Parameter::randNorm);
		function("randLogNorm", &Parameter::randLogNorm);
		function("randExp", &Parameter::randExp);
		function("randGamma", &Parameter::randGamma);
		function("randUnif", &Parameter::randUnif);
		function("densityNorm", &Parameter::densityNorm);
		function("densityLogNorm", &Parameter::densityLogNorm);


	class_<ROCParameter>( "ROCParameter" )
		.derives<Parameter>("Parameter")


		//Constructors & Destructors:
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()


		//Initialization, Restart, Index Checking:
		.method("initCovarianceMatrix", &ROCParameter::initCovarianceMatrix)
		.method("initMutationCategories", &ROCParameter::initMutationCategories) //Not an R wrapper
		.method("initSelectionCategories", &ROCParameter::initSelectionCategories) //Not an R wrapper
		.method("getCovarianceMatrixForAA", &ROCParameter::getCovarianceMatrixForAA) //Not an R wrapper
		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)

		//Prior Functions:
		.method("getMutationPriorStandardDeviation", &ROCParameter::getMutationPriorStandardDeviation)
		.method("getMutationPriorMean", &ROCParameter::getMutationPriorMean)
		.method("setMutationPriorStandardDeviation", &ROCParameter::setMutationPriorStandardDeviationR)
		.method("setMutationPriorMean", &ROCParameter::setMutationPriorMeanR)
		.method("setProposeByPrior", &ROCParameter::setProposeByPrior)

		//Posterior, Variance, and Estimates Functions:
		
		//CSP Functions:
		//Listed in the properties section below. NOTE: these getter/setters are ONLY
		//used in R

		.method("fixDM",&ROCParameter::fixDM)//Not a R wrapper
		.method("fixDEta",&ROCParameter::fixDEta)//Not a R wrapper


		//Other Functions:
		.property("proposedMutationParameter", &ROCParameter::getProposedMutationParameter,
		        &ROCParameter::setProposedMutationParameter) //R Specific
		.property("proposedSelectionParameter", &ROCParameter::getProposedSelectionParameter,
		        &ROCParameter::setProposedSelectionParameter) //R Specific
		.property("currentMutationParameter", &ROCParameter::getCurrentMutationParameter,
		        &ROCParameter::setCurrentMutationParameter) //R Specific
		.property("currentSelectionParameter", &ROCParameter::getCurrentSelectionParameter,
		        &ROCParameter::setCurrentSelectionParameter) //R Specific
		;


	class_<PAParameter>("PAParameter")
		.derives<Parameter>("Parameter")



		//Constructors & Destructors:
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()



		//Initialization, Restart, Index Checking:
		.method("initAlpha", &PAParameter::initAlphaR)
		.method("initLambdaPrime", &PAParameter::initLambdaPrimeR)
		.method("initMutationSelectionCategories", &PAParameter::initMutationSelectionCategoriesR)


		//CSP Functions:
		//Listed in the properties section below. NOTE: these getter/setters are ONLY
		//used in R

		//Other Functions:
		.method("getParameterForCategory", &PAParameter::getParameterForCategoryR)



		.property("proposedAlphaParameter", &PAParameter::getProposedAlphaParameter,
		        &PAParameter::setProposedAlphaParameter) //R Specific
		.property("proposedLambdaPrimeParameter", &PAParameter::getProposedLambdaPrimeParameter,
		        &PAParameter::setProposedLambdaPrimeParameter) //R Specific
		.property("currentAlphaParameter", &PAParameter::getCurrentAlphaParameter,
		        &PAParameter::setCurrentAlphaParameter) //R Specific
		.property("currentLambdaPrimeParameter", &PAParameter::getCurrentLambdaPrimeParameter,
		        &PAParameter::setCurrentLambdaPrimeParameter) //R Specific
		;

	class_<PANSEParameter>("PANSEParameter")
		.derives<Parameter>("Parameter")



		//Constructors & Destructors:
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()



		//Initialization, Restart, Index Checking:
		.method("initAlpha", &PANSEParameter::initAlphaR)
		.method("initLambdaPrime", &PANSEParameter::initLambdaPrimeR)
		.method("initNSERate", &PANSEParameter::initNSERateR)
		.method("initMutationSelectionCategories", &PANSEParameter::initMutationSelectionCategoriesR)
		.method("fixAlpha",&PANSEParameter::fixAlpha)
		.method("fixLambdaPrime",&PANSEParameter::fixLambdaPrime)
		.method("fixNSERate",&PANSEParameter::fixNSERate)
		.method("shareNSERate",&PANSEParameter::shareNSERate)

		//CSP Functions:
		//Listed in the properties section below. NOTE: these getter/setters are ONLY
		//used in R
		.method("initCovarianceMatrix", &PANSEParameter::initCovarianceMatrix)
		//Other Functions:
		.method("getParameterForCategory", &PANSEParameter::getParameterForCategoryR)
		.method("setPartitionFunction", &PANSEParameter::setPartitionFunction)
		.method("setTotalRFPCount", &PANSEParameter::setTotalRFPCount)
		.property("proposedAlphaParameter", &PANSEParameter::getProposedAlphaParameter,
		        &PANSEParameter::setProposedAlphaParameter) //R Specific
		.property("proposedLambdaPrimeParameter", &PANSEParameter::getProposedLambdaPrimeParameter,
		        &PANSEParameter::setProposedLambdaPrimeParameter) //R Specific
        .property("proposedNSERateParameter", &PANSEParameter::getProposedNSERateParameter,
                &PANSEParameter::setProposedNSERateParameter) //R Specific
		.property("currentAlphaParameter", &PANSEParameter::getCurrentAlphaParameter,
		        &PANSEParameter::setCurrentAlphaParameter) //R Specific
		.property("currentLambdaPrimeParameter", &PANSEParameter::getCurrentLambdaPrimeParameter,
		        &PANSEParameter::setCurrentLambdaPrimeParameter) //R Specific
        .property("currentNSERateParameter", &PANSEParameter::getCurrentNSERateParameter,
                &PANSEParameter::setCurrentNSERateParameter) //R Specific
		;

	class_<FONSEParameter>("FONSEParameter")
		.derives<Parameter>("Parameter")



		//Constructors & Destructors:
		.constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool,double>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string,double>()



		//Initialization, Restart, Index Checking:
		.method("initCovarianceMatrix", &FONSEParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &FONSEParameter::getCovarianceMatrixForAA) //Not an R wrapper
		.method("initMutation", &FONSEParameter::initMutation)
		.method("initMutationCategories", &FONSEParameter::initMutationCategories)
		.method("initSelection", &FONSEParameter::initSelection)
		.method("initSelectionCategories", &FONSEParameter::initSelectionCategories)

		.method("fixedInitiationCost",&FONSEParameter::fixedInitiationCost)

		//Prior Functions:
		.method("getMutationPriorStandardDeviation", &FONSEParameter::getMutationPriorStandardDeviation)
		.method("setMutationPriorStandardDeviation", &FONSEParameter::setMutationPriorStandardDeviation)

		.method("fixDM",&FONSEParameter::fixDM)//Not a R wrapper
		.method("fixDOmega",&FONSEParameter::fixDOmega)//Not a R wrapper


		// .property("proposedMutationParameter", &FONSEParameter::getProposedMutationParameter,
		//         &FONSEParameter::setProposedMutationParameter) //R Specific
		// .property("proposedSelectionParameter", &FONSEParameter::getProposedSelectionParameter,
		//         &FONSEParameter::setProposedSelectionParameter) //R Specific
		.property("currentMutationParameter", &FONSEParameter::getCurrentMutationParameter,
		        &FONSEParameter::setCurrentMutationParameter) //R Specific
		.property("currentSelectionParameter", &FONSEParameter::getCurrentSelectionParameter,
		        &FONSEParameter::setCurrentSelectionParameter) //R Specific
		.property("mutation_prior_sd", &FONSEParameter::getMutationPriorStandardDeviation,
		        &FONSEParameter::setMutationPriorStandardDeviation)
		;
}
#endif
