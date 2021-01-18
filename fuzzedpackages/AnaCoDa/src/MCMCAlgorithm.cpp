#include "include/MCMCAlgorithm.h"
#include <vector>
#include <random>

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


// //C++ runs only.
// #ifdef STANDALONE

// #endif


//Open MP
#ifdef _OPENMP
//#ifndef __APPLE__
#include <omp.h>
#include <thread>
#endif


//--------------------------------------------------//
//----------- Constructors & Destructors -----------//
//--------------------------------------------------//


/*MCMCAlgorithm constructor (RCPP EXPOSED)
* Arguments: None
* Sets up the object with the specified default values. Every step is a sample in
* this case, so every iteration trace values will be stored. All parameters are estimated.
*/
MCMCAlgorithm::MCMCAlgorithm() : samples(1000), thinning(1), adaptiveWidth(100 * thinning), estimateSynthesisRate(true),
	estimateCodonSpecificParameter(true), estimateHyperParameter(true)
{
	posteriorTrace.resize(samples + 1); // +1 for storing initial evaluation
	likelihoodTrace.resize(samples + 1);

    writeRestartFile = false;
	multipleFiles = false;
	fileWriteInterval = 1u;
	lastConvergenceTest = 0u;

	estimateMixtureAssignment = true;
	stepsToAdapt = -1;
}


/* MCMCAlgorithm constructor (RCPP EXPOSED)
* Arguments: number of samples wanted, thinning of iterations to get samples, how many iterations adaptive width
* can be changed, bool for where the synthesis rate, codon specific, and hyper parameters are estimated (respectively)
* Sets up the object with the given parameters. Adaptive with is actually set to adaptive width * thinning. NOTE:
* hyper parameters normally affect synthesis rate parameters, so either both should be turned on or neither.
*/
MCMCAlgorithm::MCMCAlgorithm(unsigned _samples, unsigned _thinning, unsigned _adaptiveWidth, bool _estimateSynthesisRate,
							 bool _estimateCodonSpecificParameter, bool _estimateHyperParameter) : samples(_samples),
							 thinning(_thinning), adaptiveWidth(_adaptiveWidth * thinning),
							 estimateSynthesisRate(_estimateSynthesisRate),
							 estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
							 estimateHyperParameter(_estimateHyperParameter)
{
	posteriorTrace.resize(samples + 1);// +1 for storing initial evaluation
	likelihoodTrace.resize(samples + 1);
	writeRestartFile = false;
	multipleFiles = false;
	fileWriteInterval = 1u;
	lastConvergenceTest = 0u;
	estimateMixtureAssignment = true;
	stepsToAdapt = -1;
}


/* MCMCAlgorithm destructor (NOT EXPOSED)
* Arguments: None
* Standard destructor for the object.
*/
MCMCAlgorithm::~MCMCAlgorithm()
{
	//dtor
}





//----------------------------------------------------//
//---------- Acceptance Rejection Functions ----------//
//----------------------------------------------------//

/* acceptRejectSynthesisRateLevelForAllGenes (NOT EXPOSED)
* Arguments: reference to a genome and a model. which iteration (step) is currently being estimated
* Based on the loglikelihood probabilities proposed for a gene (with and without reverse jump), it is decided
* whether or not a synthesis rate is accepted for that gene or not. Mixture assignment is also updated when
* applicable.
*/
double MCMCAlgorithm::acceptRejectSynthesisRateLevelForAllGenes(Genome& genome, Model& model, int iteration)
{
    
	double loglikelihood = 0.0;
	double logPosterior = 0.0;
    //double logPosterior2 = 0.0; //currently unused
	int numGenes = genome.getGenomeSize();

	unsigned numSynthesisRateCategories = model.getNumSynthesisRateCategories();
	unsigned numMixtures = model.getNumMixtureElements();
	std::vector <double> dirichletParameters(numMixtures, 0);

	//initialize parameter's size
	for (unsigned i = 0u; i < numGenes; i++)
	{
		Gene *gene = &genome.getGene(i);

		/*
			 Since some values returned by calculatelogPosteriorRatioPerGene are very small (~ -1100),
			 exponentiation leads to 0.
			 To solve this problem, we adjust the value by a constant c.
			 I choose to use the average value across all mixtures.
			 We justify this by
			 P = Sum(p_i*f(...))
			 => f' = c*f
			 => ln(f') = ln(c) + ln(f)
			 => ln(P) = ln( Sum(p_i*f'(...)) )
			 => ln(P) = ln(P') - ln(c)
			 Note that we use the inverse sign because our values of ln(f) and ln(f') are negative.
		 */

		double maxValue = -1.0e+20;
		double maxValue2 = -1.0e+20;
		unsigned mixtureIndex = 0u;

		std::vector <double> unscaledLogProb_curr(numSynthesisRateCategories, 0.0);
		std::vector <double> unscaledLogProb_prop(numSynthesisRateCategories, 0.0);

		//Added by Alex
		std::vector <double> unscaledLogLike_curr(numSynthesisRateCategories, 0.0);
		std::vector <double> unscaledLogLike_prop(numSynthesisRateCategories, 0.0);

		std::vector <double> unscaledLogPost_curr(numSynthesisRateCategories, 0.0);
		std::vector <double> unscaledLogPost_prop(numSynthesisRateCategories, 0.0);

		std::vector <double> unscaledLogProb_curr_singleMixture(numMixtures, 0.0);
		std::vector <double> unscaledLogProb_prop_singleMixture(numMixtures, 0.0);
		std::vector <double> probabilities(numMixtures, 0.0);

		for (unsigned k = 0u; k < numSynthesisRateCategories; k++)
		{
			// logProbabilityRatio contains the logProbabilityRatio in element 0,
			// the current unscaled probability in element 1 and the proposed unscaled probability in element 2
			std::vector<unsigned> mixtureElements = model.getMixtureElementsOfSelectionCategory(k);
			for (unsigned n = 0u; n < mixtureElements.size(); n++)
			{
				unsigned mixtureElement = mixtureElements[n];
				double logProbabilityRatio[7];
				model.calculateLogLikelihoodRatioPerGene(*gene, i, mixtureElement, logProbabilityRatio); //has to return likelihood, not just ratio, no priors

				// log posterior with and without rev. jump probability
				unscaledLogProb_curr[k] += logProbabilityRatio[1]; // with rev. jump prob.
				unscaledLogProb_prop[k] += logProbabilityRatio[2]; // with rev. jump prob.

				// TODO: unscaledLogPost_curr and unscaledLogPost_prop have been made redundant by modification to code adding unscaledLogProb_curr_singleMixture
				// unscaledLogProb_prop_singleMixture. Remove.
				unscaledLogPost_curr[k] += logProbabilityRatio[3]; // without rev. jump prob.
				unscaledLogPost_prop[k] += logProbabilityRatio[4]; // without rev. jump prob.
				unscaledLogLike_curr[k] += logProbabilityRatio[5]; //current logLikelihood
				unscaledLogLike_prop[k] += logProbabilityRatio[6]; //proposed logLikelihood

				unscaledLogProb_curr_singleMixture[mixtureIndex] = logProbabilityRatio[3];
				unscaledLogProb_prop_singleMixture[mixtureIndex] = logProbabilityRatio[4];

				maxValue = unscaledLogProb_curr_singleMixture[mixtureIndex] > maxValue ?
						   unscaledLogProb_curr_singleMixture[mixtureIndex] : maxValue;

				maxValue2 = unscaledLogProb_prop_singleMixture[mixtureIndex] > maxValue2 ?
							unscaledLogProb_prop_singleMixture[mixtureIndex] : maxValue2;
				//maxValue2 = unscaledLogPost_prop[k] > maxValue2 ?
				//		unscaledLogPost_prop[k] : maxValue2;

				mixtureIndex++;
			}
		}
		maxValue2 = maxValue > maxValue2 ? maxValue : maxValue2;

		// adjust the the unscaled probabilities by the constant c
		// ln(f') = ln(c) + ln(f)
		// calculate ln(P) = ln( Sum(p_i*f'(...)) ) and obtain normalizing constant for new p_i
		double normalizingProbabilityConstant = 0.0;

		for (unsigned k = 0u; k < numMixtures; k++)
		{
			unscaledLogProb_curr_singleMixture[k] -= maxValue;
			probabilities[k] = model.getCategoryProbability(k) * std::exp(unscaledLogProb_curr_singleMixture[k]);
			normalizingProbabilityConstant += probabilities[k];
            if (std::isnan(probabilities[k]))
            {
                my_print("\n\n\n");
                my_print("Gene: %, %\n", i, gene->getId());
                my_print("Mixture: %\n", k);
                my_print("unscaled Mix. Prop. of gene: %\n", probabilities[k]);
                my_print("Mix. cat. prob.: %\n", model.getCategoryProbability(k));
                my_print("curr logLik - maxVal: %\n", unscaledLogProb_curr_singleMixture[k]);
                my_print("exp(curr logLik - maxVal): %\n", std::exp(unscaledLogProb_curr_singleMixture[k]));
                my_print("Max Val. %\n", maxValue);
                my_print("\n\n\n");
				break;
            }
            unscaledLogProb_curr_singleMixture[k] += maxValue;
		}
		// normalize probabilities
		for (unsigned k = 0u; k < numMixtures; k++)
		{
			probabilities[k] = probabilities[k] / normalizingProbabilityConstant;
		}

		double currGeneLogPost = 0.0;
		for (unsigned k = 0u; k < numSynthesisRateCategories; k++)
		{
			// that we loop over numSynthesisRateCategories instead of numMixtures is not a problem
			// in the case of selectionShared, since we sum above over each mixture sharing a category.

			// We do not need to add std::log(model.getCategoryProbability(k)) since it will cancel in the ratio!
			double currLogPost = unscaledLogProb_curr[k];
			double propLogPost = unscaledLogProb_prop[k];
			std::vector<unsigned> mixtureElements = model.getMixtureElementsOfSelectionCategory(k);
            double alpha = -Parameter::randExp(1);
         

			if ( (alpha < (propLogPost - currLogPost)) && estimateSynthesisRate )
			{
				model.updateSynthesisRate(i,k);

				for (unsigned n = 0u; n < mixtureElements.size(); n++)
				{
					unsigned element = mixtureElements[n];
					currGeneLogPost += probabilities[element] * std::exp(unscaledLogProb_prop_singleMixture[element] - maxValue2);
				}
				loglikelihood += unscaledLogLike_prop[k]; //unscaled
			}
			else
			{
				for (unsigned n = 0u; n < mixtureElements.size(); n++)
				{
					unsigned element = mixtureElements[n];
					currGeneLogPost += probabilities[element] * std::exp(unscaledLogProb_curr_singleMixture[element] - maxValue2);
				}
				loglikelihood += unscaledLogLike_curr[k]; //unscaled
			}

            if (std::isnan(logPosterior))
            {
                my_print("\n\n\n");
                my_print("logPosterior: %\n", logPosterior);
                my_print("currGeneLogPost: %\n", currGeneLogPost);
                my_print("Gene: %, %\n", i, gene->getId());
                my_print("Mixture: %\n", k);
                my_print("Mix. Prop.: %\n", probabilities[k]);
                my_print("cur logLik: %\n", unscaledLogPost_curr[k]);
                my_print("prop logLik: %\n", unscaledLogPost_prop[k]);
                my_print("Accepted?: % < % , %\n", alpha, (propLogPost - currLogPost), alpha < (propLogPost - currLogPost));
                my_print("\n\n\n");
            }
		}
		logPosterior += std::log(currGeneLogPost) + maxValue2;

		// Get category in which the gene is placed in.
		// If we use multiple sequence observation (like different mutants),
		// randMultinom needs a parameter N to place N observations in numMixture buckets
		unsigned categoryOfGene = Parameter::randMultinom(probabilities, numMixtures);
		if (estimateMixtureAssignment)
			model.setMixtureAssignment(i, categoryOfGene);

		dirichletParameters[categoryOfGene] += 1;
		if ((iteration % thinning) == 0)
		{
			model.updateSynthesisRateTrace(iteration/thinning, i);
			model.updateMixtureAssignmentTrace(iteration/thinning, i);
		}

	}

	// take all priors into account
	//loglikelihood = logPosterior;
	logPosterior += model.calculateAllPriors();
    if (std::isnan(logPosterior))
    {
        my_print("\n\n\n");
        my_print("logPosterior NaN after addition of prior values!\n");
        my_print("\n\n\n");
    }
	std::vector <double> newMixtureProbabilities(numMixtures, 0);
	Parameter::randDirichlet(dirichletParameters, numMixtures, newMixtureProbabilities);
	for (unsigned k = 0u; k < numMixtures; k++)
	{
		model.setCategoryProbability(k, newMixtureProbabilities[k]);
	}
	if ((iteration % thinning) == 0)
	{
		model.updateMixtureProbabilitiesTrace(iteration/thinning);
		likelihoodTrace[iteration/thinning] = loglikelihood;
	}
	return logPosterior;
}



/* acceptRejectCodonSpecificParameter (NOT EXPOSED)
 * Arguments: reference to a genome and a model. which iteration (step) is currently being estimated
 * Calculates the logLikelihood for each grouping based on codon specific parameters. If this is greater
 * than a random number from the exponential distribution we update the parameters from proposed to
 * current. Update the trace when applicable.
*/
void MCMCAlgorithm::acceptRejectCodonSpecificParameter(Genome& genome, PANSEModel& model, int iteration)
{
	std::vector<double> acceptanceRatioForAllMixtures(5,0.0);
	unsigned size = model.getGroupListSize();

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine e(seed);

	std::vector<unsigned> groups(size);
	std::iota(groups.begin(),groups.end(),0);
	std::shuffle ( groups.begin(), groups.end(),e);

	std::string param_1,param_2;
	std::uniform_int_distribution<int> uni(0,1);

	int x = uni(e);

	if (x == 0)
	{
		param_1 = "Elongation";
		param_2 = "NSE";
	}
	else
	{
		param_2 = "Elongation";
		param_1 = "NSE";
	}
	
	for (unsigned i = 0; i < size; i++)
	{

		std::string grouping = model.getGrouping(groups[i]);
		model.calculateLogLikelihoodRatioPerGroupingPerCategory(grouping, genome, acceptanceRatioForAllMixtures,param_1);
    	double threshold = -Parameter::randExp(1);
 		if (threshold < acceptanceRatioForAllMixtures[0] && std::isfinite(acceptanceRatioForAllMixtures[0]) && !std::isnan(acceptanceRatioForAllMixtures[2]))
		{	
			if (std::isnan(acceptanceRatioForAllMixtures[0]))
			{
				my_print("ERROR: Accepted proposed value that results in NaN\n");
			}
			
			// moves proposed codon specific parameters to current codon specific parameters
			model.updateCodonSpecificParameter(grouping,param_1);
		}	
	}

	std::shuffle ( groups.begin(), groups.end(),e);
	for (unsigned i = 0; i < size; i++)
	{

		std::string grouping = model.getGrouping(groups[i]);
		// calculate likelihood ratio for every Category for current AA (ROC, FONSE) or Codon (PA, PANSE)
		model.calculateLogLikelihoodRatioPerGroupingPerCategory(grouping, genome, acceptanceRatioForAllMixtures,param_2);
		double threshold = -Parameter::randExp(1);
		if (threshold < acceptanceRatioForAllMixtures[0] && std::isfinite(acceptanceRatioForAllMixtures[0]) && !std::isnan(acceptanceRatioForAllMixtures[2]))
		{	
			// moves proposed codon specific parameters to current codon specific parameters
			if (std::isnan(acceptanceRatioForAllMixtures[0]))
			{
				my_print("ERROR: Accepted proposed value that results in NaN\n");
			}
			model.updateCodonSpecificParameter(grouping,param_2);
			if ((iteration % thinning) == 0)
			{
				likelihoodTrace[(iteration / thinning)] = acceptanceRatioForAllMixtures[2];//will be 0
				posteriorTrace[(iteration / thinning)] = acceptanceRatioForAllMixtures[4];//will be 0
			}
		}
		else
		{
			if ((iteration % thinning) == 0)
			{
				likelihoodTrace[(iteration / thinning)] = acceptanceRatioForAllMixtures[1];
				posteriorTrace[(iteration / thinning)] = acceptanceRatioForAllMixtures[3];
			}	
		}
	}
	

	if ((iteration % thinning) == 0)
	{
		for (unsigned i = 0;i < size; i++)
		{
			std::string grouping = model.getGrouping(i);
			model.updateCodonSpecificParameterTrace(iteration/thinning, grouping);
		}
	}
}



/* acceptRejectCodonSpecificParameter (NOT EXPOSED)
 * Arguments: reference to a genome and a model. which iteration (step) is currently being estimated
 * Calculates the logLikelihood for each grouping based on codon specific parameters. If this is greater
 * than a random number from the exponential distribution we update the parameters from proposed to
 * current. Update the trace when applicable.
*/
void MCMCAlgorithm::acceptRejectCodonSpecificParameter(Genome& genome, Model& model, int iteration)
{
	std::vector<double> acceptanceRatioForAllMixtures(5,0.0);
	unsigned size = model.getGroupListSize();

	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// std::default_random_engine e(seed);

	// std::random_device rd;
 //    std::mt19937 g(rd());

	// std::vector<unsigned> groups(size);
	// std::iota(groups.begin(),groups.end(),0);
	// std::shuffle ( groups.begin(), groups.end(),g);
	
	for (unsigned i = 0; i < size; i++)
	{

		std::string grouping = model.getGrouping(i);
		// calculate likelihood ratio for every Category for current AA (ROC, FONSE) or Codon (PA, PANSE)
		model.calculateLogLikelihoodRatioPerGroupingPerCategory(grouping, genome, acceptanceRatioForAllMixtures);
    	double threshold = -Parameter::randExp(1);
 		if (threshold < acceptanceRatioForAllMixtures[0] && std::isfinite(acceptanceRatioForAllMixtures[0]) && !std::isnan(acceptanceRatioForAllMixtures[2]))
		{	
			
			
			// moves proposed codon specific parameters to current codon specific parameters
			if (std::isnan(acceptanceRatioForAllMixtures[0]))
			{
				my_print("ERROR: Accepted proposed value that results in NaN\n");
			}
			model.updateCodonSpecificParameter(grouping);
			if ((iteration % thinning) == 0)
			{
				likelihoodTrace[(iteration / thinning)] += acceptanceRatioForAllMixtures[2];//will be 0
				posteriorTrace[(iteration / thinning)] += acceptanceRatioForAllMixtures[4];//will be 0
			}
			
		}
		else
		{
			if ((iteration % thinning) == 0)
			{
				likelihoodTrace[(iteration / thinning)] += acceptanceRatioForAllMixtures[1];
				posteriorTrace[(iteration / thinning)] += acceptanceRatioForAllMixtures[3];

			}
			
		}
	
		// if ((iteration % thinning) == 0)
		// {
		// 	model.updateCodonSpecificParameterTrace(iteration/thinning, grouping);
		// }
	}
	if ((iteration % thinning) == 0)
	{
		for (unsigned i = 0;i < size; i++)
		{
			std::string grouping = model.getGrouping(i);
			model.updateCodonSpecificParameterTrace(iteration/thinning, grouping);
		}
	}
}





/* acceptRejectHyperParameter (NOT EXPOSED)
 * Arguments: reference to a genome and a model. which iteration (step) is currently being estimated
 * Calculates the logLikelihood for hyper parameters. If the calculated value is greater than a random number
 * from the exponential distribution we update the parameters from proposed to current. Update the trace when applicable.
*/
void MCMCAlgorithm::acceptRejectHyperParameter(Genome &genome, Model& model, unsigned iteration)
{
	std::vector <double> logProbabilityRatios;
	model.calculateLogLikelihoodRatioForHyperParameters(genome, iteration, logProbabilityRatios);
	for (unsigned i = 0; i < logProbabilityRatios.size(); i++)
	{
		//my_print("%\n",logProbabilityRatios[i]);
		if (!std::isfinite(logProbabilityRatios[i]))
            my_print("logProbabilityRatio % not finite!\n", i);

		if (-Parameter::randExp(1) < logProbabilityRatios[i])
			model.updateHyperParameter(i);
	}

	if ((iteration % thinning) == 0)
		model.updateHyperParameterTraces(iteration/thinning);
}





//------------------------------------//
//---------- MCMC Functions ----------//
//------------------------------------//

/* run (RCPP EXPOSED)
 * Arguments: reference to a genome and a model. number of cores to run on (unless running on a MAC). Number of
 * iterations to allow initial conditions to vary.
 * Runs the MCMC algorithm for the set number of iterations. Run can terminate early if Geweke score meets
 * the set criteria.
*/
void MCMCAlgorithm::run(Genome& genome, Model& model, unsigned numCores, unsigned divergenceIterations)
{
#ifdef _OPENMP
//#ifndef __APPLE__
	omp_set_num_threads(numCores);
#endif

	// Replace with reportSample?
	unsigned reportStep = (100u < thinning) ? thinning : 100u;

	// Allows to diverge from initial conditions (divergenceIterations controls the divergence).
	// This allows for varying initial conditions for better exploration of the parameter space.
	varyInitialConditions(genome, model, divergenceIterations);

	unsigned maximumIterations = samples * thinning;
	// initialize everything

	//model.setNumPhiGroupings(genome.getGene(0).getObservedSynthesisRateValues().size());
	model.initTraces(samples + 1, genome.getGenomeSize(),(estimateSynthesisRate||estimateMixtureAssignment)); //Samples + 2 so we can store the starting and ending values.
	// starting the MCMC
	
	model.updateTracesWithInitialValues(genome);
	if (stepsToAdapt == -1)
		stepsToAdapt = maximumIterations;
	my_print("Type: %\n",model.type);
	my_print("Starting MCMC\n");
	my_print("\tEstimate Codon Specific Parameters? % \n", (estimateCodonSpecificParameter ? "TRUE" : "FALSE") );
	my_print("\tEstimate Hyper Parameters? % \n", (estimateHyperParameter ? "TRUE" : "FALSE") );
	my_print("\tEstimate Synthesis rates? % \n", (estimateSynthesisRate ? "TRUE" : "FALSE") );
	my_print("\tStarting MCMC with % iterations\n", maximumIterations);
	my_print("\tAdapting will stop after % steps\n", stepsToAdapt);

	// set the last iteration to the max iterations,
	// this way if the MCMC doesn't exit based on Geweke score, it will use the max iteration for posterior means
	model.setLastIteration(samples);
	for (int iteration = 1u; iteration <= maximumIterations; iteration++)
	{
		if (writeRestartFile)
		{
			if ((iteration) % fileWriteInterval == 0u)
			{
				my_print("Begin saving restart file(s) at sample (iteration): % (%)\n",  (iteration / thinning), iteration);

				if (multipleFiles)
				{
					std::ostringstream oss;
					oss << file << "_" << (iteration) / thinning;
					std::string tmp = oss.str();
					model.writeRestartFile(tmp);
				}
				else
				{
					model.writeRestartFile(file);
				}
			}
		}
		if ((iteration) % reportStep == 0u)
		{
            #ifndef STANDALONE
            Rcpp::checkUserInterrupt();
            #endif
	   		my_print("Status at thinned sample (iteration): % (%)\n",  (iteration / thinning), iteration);
			my_print("\t current logPosterior: % \n", posteriorTrace[(iteration/thinning) - 1]);
			my_print("\t current logLikelihood: %\n", likelihoodTrace[(iteration/thinning) - 1]);
			if (iteration > stepsToAdapt)
				my_print("No longer adapting\n");

			model.printHyperParameters();
			for (unsigned i = 0u; i < model.getNumMixtureElements(); i++)
			{
				my_print("\t current Mixture element probability for element %: %\n", i, model.getCategoryProbability(i));
			}
		}
		if (estimateCodonSpecificParameter)
		{
			model.proposeCodonSpecificParameter();
			acceptRejectCodonSpecificParameter(genome, model, iteration);
			
            //TODO:Probably do a nan check
			if ((iteration % adaptiveWidth) == 0u)
				model.adaptCodonSpecificParameterProposalWidth(adaptiveWidth, iteration / thinning, iteration <= stepsToAdapt);
		}
		// update hyper parameter
		if (estimateHyperParameter)
		{
			model.updateGibbsSampledHyperParameters(genome);
			model.proposeHyperParameters();
			acceptRejectHyperParameter(genome, model, iteration);
            //TODO:Probably do a nan check
			if ((iteration % adaptiveWidth) == 0u)
				model.adaptHyperParameterProposalWidths(adaptiveWidth, iteration <= stepsToAdapt);

		}
		// update expression level values
		if (estimateSynthesisRate || estimateMixtureAssignment)
		{
			model.proposeSynthesisRateLevels();
			double logPost = acceptRejectSynthesisRateLevelForAllGenes(genome, model, iteration);
			if ((iteration % thinning) == 0u)
			{
				posteriorTrace[(iteration / thinning)] = logPost;
				if (std::isnan(logPost))
				{
					my_printError("ERROR: LogPosterior is NaN, exiting at iteration %\n", iteration);
					model.setLastIteration(iteration / thinning);
					return;
				}
			}
			if ((iteration % adaptiveWidth) == 0u)
				model.adaptSynthesisRateProposalWidth(adaptiveWidth, iteration <= stepsToAdapt);
  		}



		if ((iteration % (50 * adaptiveWidth)) == 0u)
		{
			double gewekeScore = calculateGewekeScore(iteration/thinning);

			my_print("##################################################\n");
			my_print("Geweke Score after % iterations: %\n", iteration, gewekeScore);
			my_print("##################################################\n");

			if (std::abs(gewekeScore) < 1.96)
			{
				my_print("Stopping run based on convergence after % iterations\n\n", iteration);

				// Comment out this break to keep the run from stopping on convergence
				//model.setLastIteration(iteration/thinning);
				//break;
			}
		}
	} // end MCMC loop

	if (writeRestartFile)
	{
		std::ostringstream oss;
		oss << file << "_final";
		std::string tmp = oss.str();
		model.writeRestartFile(tmp);
	}
	my_print("leaving MCMC loop\n");
}

/* run (RCPP EXPOSED)
 * Arguments: reference to a genome and a model. number of cores to run on (unless running on a MAC). Number of
 * iterations to allow initial conditions to vary.
 * Runs the MCMC algorithm for the set number of iterations. Run can terminate early if Geweke score meets
 * the set criteria.
*/
void MCMCAlgorithm::run_PANSE(Genome& genome, PANSEModel& model, unsigned numCores, unsigned divergenceIterations)
{
#ifdef _OPENMP
//#ifndef __APPLE__
	omp_set_num_threads(numCores);
#endif

	// Replace with reportSample?
	unsigned reportStep = (100u < thinning) ? thinning : 100u;

	// Allows to diverge from initial conditions (divergenceIterations controls the divergence).
	// This allows for varying initial conditions for better exploration of the parameter space.
	varyInitialConditions(genome, model, divergenceIterations);

	unsigned maximumIterations = samples * thinning;
	// initialize everything

	//model.setNumPhiGroupings(genome.getGene(0).getObservedSynthesisRateValues().size());
	model.initTraces(samples + 1, genome.getGenomeSize(),(estimateSynthesisRate||estimateMixtureAssignment)); //Samples + 2 so we can store the starting and ending values.
	// starting the MCMC
	
	model.updateTracesWithInitialValues(genome);
	if (stepsToAdapt == -1)
		stepsToAdapt = maximumIterations;
	my_print("Type: %\n",model.type);
	my_print("Starting MCMC\n");
	my_print("\tEstimate Codon Specific Parameters? % \n", (estimateCodonSpecificParameter ? "TRUE" : "FALSE") );
	my_print("\tEstimate Hyper Parameters? % \n", (estimateHyperParameter ? "TRUE" : "FALSE") );
	my_print("\tEstimate Synthesis rates? % \n", (estimateSynthesisRate ? "TRUE" : "FALSE") );
	my_print("\tStarting MCMC with % iterations\n", maximumIterations);
	my_print("\tAdapting will stop after % steps\n", stepsToAdapt);

	// set the last iteration to the max iterations,
	// this way if the MCMC doesn't exit based on Geweke score, it will use the max iteration for posterior means
	model.setLastIteration(samples);
	for (int iteration = 1u; iteration <= maximumIterations; iteration++)
	{
		if (writeRestartFile)
		{
			if ((iteration) % fileWriteInterval == 0u)
			{
				my_print("Begin saving restart file(s) at sample (iteration): % (%)\n",  (iteration / thinning), iteration);

				if (multipleFiles)
				{
					std::ostringstream oss;
					oss << file << "_" << (iteration) / thinning;
					std::string tmp = oss.str();
					model.writeRestartFile(tmp);
				}
				else
				{
					model.writeRestartFile(file);
				}
			}
		}
		if ((iteration) % reportStep == 0u)
		{
            #ifndef STANDALONE
            Rcpp::checkUserInterrupt();
            #endif
	   		my_print("Status at thinned sample (iteration): % (%)\n",  (iteration / thinning), iteration);
			my_print("\t current logPosterior: % \n", posteriorTrace[(iteration/thinning) - 1]);
			my_print("\t current logLikelihood: %\n", likelihoodTrace[(iteration/thinning) - 1]);
			if (iteration > stepsToAdapt)
				my_print("No longer adapting\n");

			model.printHyperParameters();
			for (unsigned i = 0u; i < model.getNumMixtureElements(); i++)
			{
				my_print("\t current Mixture element probability for element %: %\n", i, model.getCategoryProbability(i));
			}
		}
		if (estimateCodonSpecificParameter)
		{
			model.proposeCodonSpecificParameter();
			acceptRejectCodonSpecificParameter(genome, model, iteration);
			
            //TODO:Probably do a nan check
			if ((iteration % adaptiveWidth) == 0u)
				model.adaptCodonSpecificParameterProposalWidth(adaptiveWidth, iteration / thinning, iteration <= stepsToAdapt);
		}
		// update hyper parameter
		if (estimateHyperParameter)
		{
			model.updateGibbsSampledHyperParameters(genome);
			model.proposeHyperParameters();
			acceptRejectHyperParameter(genome, model, iteration);
            //TODO:Probably do a nan check
			if ((iteration % adaptiveWidth) == 0u)
				model.adaptHyperParameterProposalWidths(adaptiveWidth, iteration <= stepsToAdapt);

		}
		// update expression level values
		if (estimateSynthesisRate || estimateMixtureAssignment)
		{
			model.proposeSynthesisRateLevels();
			double logPost = acceptRejectSynthesisRateLevelForAllGenes(genome, model, iteration);
			if ((iteration % thinning) == 0u)
			{
				posteriorTrace[(iteration / thinning)] = logPost;
				if (std::isnan(logPost))
				{
					my_printError("ERROR: LogPosterior is NaN, exiting at iteration %\n", iteration);
					model.setLastIteration(iteration / thinning);
					return;
				}
			}
			if ((iteration % adaptiveWidth) == 0u)
				model.adaptSynthesisRateProposalWidth(adaptiveWidth, iteration <= stepsToAdapt);
  		}



		if ((iteration % (50 * adaptiveWidth)) == 0u)
		{
			double gewekeScore = calculateGewekeScore(iteration/thinning);

			my_print("##################################################\n");
			my_print("Geweke Score after % iterations: %\n", iteration, gewekeScore);
			my_print("##################################################\n");

			if (std::abs(gewekeScore) < 1.96)
			{
				my_print("Stopping run based on convergence after % iterations\n\n", iteration);

				// Comment out this break to keep the run from stopping on convergence
				//model.setLastIteration(iteration/thinning);
				//break;
			}
		}
	} // end MCMC loop

	if (writeRestartFile)
	{
		std::ostringstream oss;
		oss << file << "_final";
		std::string tmp = oss.str();
		model.writeRestartFile(tmp);
	}
	my_print("leaving MCMC loop\n");
}





/* varyInitialConditions (NOT EXPOSED)
 * Arguments: reference to a genome and a model. Number of iterations that the model can diverge from initial
 * conditions.
 * Allows to diverge from initial conditions (divergenceIterations controls the divergence).
 * This allows for varying initial conditions for better exploration of the parameter space.
 * This functions does not take the likelihood into account, the random walk is only constrained by the prior
 * for each parameter. If no prior (flat prior) is present, an unconstrained random walk is performed.
*/
void MCMCAlgorithm::varyInitialConditions(Genome& genome, Model& model, unsigned divergenceIterations)
{
	// Currently unused variables
	//double previous_post = 0.0;
	//double current_post = 0.0;

	// NOTE: IF PRIORS ARE ADDED, TAKE INTO ACCOUNT HERE!
	my_print("Allowing divergence from initial conditions for % iterations.\n\n", divergenceIterations);
	// divergence from initial conditions is not stored in trace

	// how many steps do you want to walk "away" from the initial conditions
	for (unsigned iteration = 0u; iteration < divergenceIterations; iteration++)
	{
		if (estimateCodonSpecificParameter)
		{
			model.proposeCodonSpecificParameter();
		}
		if (estimateHyperParameter)
		{
			model.proposeHyperParameters();
		}
		if (estimateSynthesisRate)
		{
			model.proposeSynthesisRateLevels();
		}
		// no prior on codon specific parameters -> just accept everything
		if (estimateCodonSpecificParameter)
		{
			unsigned size = model.getGroupListSize();
			for (unsigned i = 0; i < size; i++)
			{
				std::string grouping = model.getGrouping(i);
				model.updateCodonSpecificParameter(grouping);
			}
			model.completeUpdateCodonSpecificParameter();
		}
		// no prior on hyper parameters -> just accept everything
		if (estimateHyperParameter)
		{
			model.updateAllHyperParameter();
		}
		// prior on phi values -> take prior into account, but only the prior no likelihood
		if (estimateSynthesisRate)
		{
			unsigned numGenes = genome.getGenomeSize();
			unsigned numSynthesisRateCategories = model.getNumSynthesisRateCategories();
			for (unsigned i = 0u; i < numGenes; i++)
			{
				for (unsigned k = 0u; k < numSynthesisRateCategories; k++)
				{
					// map from mixture to category and obtain corresponding phi value
					unsigned expressionCategory = model.getSynthesisRateCategory(k);
					double phiValue = model.getSynthesisRate(i, expressionCategory, false);
					double phiValue_proposed = model.getSynthesisRate(i, expressionCategory, true);

					unsigned mixture = model.getMixtureAssignment(k);
					mixture = model.getSynthesisRateCategory(mixture);
					double stdDevSynthesisRate = model.getStdDevSynthesisRate(mixture, false);
					double mPhi = (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2);

					// accept/ reject based on prior ratio
					double logPhiProbability = Parameter::densityLogNorm(phiValue, mPhi, stdDevSynthesisRate, true);
					double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, mPhi, stdDevSynthesisRate, true);
					if ( -Parameter::randExp(1) < (logPhiProbability_proposed - logPhiProbability) )
						model.updateSynthesisRate(i, k);
				}
			}
		}

		// update Gibbs sampled parameter.
		if (estimateHyperParameter)
		{
			model.updateGibbsSampledHyperParameters(genome);
		}
	}
}


/* calculateGewekeScore (NOT EXPOSED)
 * Arguments: current iteration
 * Calculate the Geweke score based of the last test and the posterior means. If this is the first
 * convergence test, lastConvergenceTest should be equal to 0.
*/
double MCMCAlgorithm::calculateGewekeScore(unsigned current_iteration)
{
	double posteriorMean1 = 0.0;
	double posteriorMean2 = 0.0;
	double posteriorVariance1 = 0.0;
	double posteriorVariance2 = 0.0;

	unsigned end1 = (unsigned)std::round( (current_iteration - lastConvergenceTest) * 0.1) + lastConvergenceTest;
	unsigned start2 = (unsigned)std::round(current_iteration - (current_iteration * 0.5));

	double numSamples1 = (double) (end1 - lastConvergenceTest);
	double numSamples2 = std::round(current_iteration * 0.5);

	// calculate mean and and variance of first part of likelihood trace
	for (unsigned i = lastConvergenceTest; i < end1; i++)
	{
		posteriorMean1 += posteriorTrace[i];
	}
	posteriorMean1 = posteriorMean1 / numSamples1;
	for (unsigned i = lastConvergenceTest; i < end1; i++)
	{
		posteriorVariance1 += (posteriorTrace[i] - posteriorMean1) * (posteriorTrace[i] - posteriorMean1);
	}
	posteriorVariance1 = posteriorVariance1 / numSamples1;


	// calculate mean and and variance of last part of likelihood trace
	for (unsigned i = start2; i < current_iteration; i++)
	{
		posteriorMean2 += posteriorTrace[i];
	}
	posteriorMean2 = posteriorMean2 / numSamples2;
	for (unsigned i = start2; i < current_iteration; i++)
	{
		posteriorVariance2 += (posteriorTrace[i] - posteriorMean2) * (posteriorTrace[i] - posteriorMean2);
	}
	posteriorVariance2 = posteriorVariance2 / numSamples2;

	lastConvergenceTest = current_iteration;
	// Geweke score
	return (posteriorMean1 - posteriorMean2) / std::sqrt( ( posteriorVariance1 / numSamples1 ) + ( posteriorVariance2 / numSamples2 ) );
}


/* isEstimateSynthesisRate (NOT EXPOSED)
 * Arguments: None
 * Return the boolean value for if synthesis rate should be estimated in this run.
*/
bool MCMCAlgorithm::isEstimateSynthesisRate()
{
	return estimateSynthesisRate;
}


/* isEstimateCodonSpecificParameter (NOT EXPOSED)
 * Arguments: None
 * Return the boolean value for if codon specific parameters should be estimated in this run.
*/
bool MCMCAlgorithm::isEstimateCodonSpecificParameter()
{
	return estimateCodonSpecificParameter;
}


/* isEstimateHyperParameter (NOT EXPOSED)
 * Arguments: None
 * Return the boolean value for if hyper parameters should be estimated in this run.
*/
bool MCMCAlgorithm::isEstimateHyperParameter()
{
	return estimateHyperParameter;
}


/* isEstimateMixtureAssignment (NOT EXPOSED)
 * Arguments: None
 * Return the boolean value for if mixture assignment should be estimated in this run.
*/
bool MCMCAlgorithm::isEstimateMixtureAssignment()
{
	return estimateMixtureAssignment;
}


/* setEstimateSynthesisRate (NOT EXPOSED)
 * Arguments: boolean describing if the parameter should be estimated.
 * Sets estimateSynthesisRate to the specified value.
*/
void MCMCAlgorithm::setEstimateSynthesisRate(bool in)
{
	estimateSynthesisRate = in;
}


/* setEstimateCodonSpecificParameter (NOT EXPOSED)
 * Arguments: boolean describing if the parameter should be estimated.
 * Sets estimateCodonSpecificParameter to the specified value.
*/
void MCMCAlgorithm::setEstimateCodonSpecificParameter(bool in)
{
	estimateCodonSpecificParameter = in;
}


/* setEstimateHyperParameter (NOT EXPOSED)
 * Arguments: boolean describing if the parameter should be estimated.
 * Sets estimateHyperParameter to the specified value.
*/
void MCMCAlgorithm::setEstimateHyperParameter(bool in)
{
	estimateHyperParameter = in;
}


/* setEstimateMixtureAssignment (RCPP EXPOSED)
 * Arguments: boolean describing if the parameter should be estimated.
 * Sets estimateMixtureAssignment to the specified value.
*/
void MCMCAlgorithm::setEstimateMixtureAssignment(bool in)
{
	estimateMixtureAssignment = in;
}


/* setRestartFileSettings (RCPP EXPOSED)
 * Arguments: base of filename to write to, interval to write file, bool telling if one or many files should be written.
 * Filename is the base of the file name to be written to. The interval will be placed in front to specify the step
 * it was written at and allow for multiple files to be written if multiple is set to true. The file write interval
 * is actually the interval specified times the thinning.
*/
//' @name setRestartFileSettings
//' @title setRestartFileSettings
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Set restart file output name and frequency prior to running MCMC
//' @param filename name of restart file
//' @param interval number of samples (ie. iterations * thinning) between writing new restart file
//' @param multiple if true, will output a new restart file at each interval (file name will include sample it was written at)

void MCMCAlgorithm::setRestartFileSettings(std::string filename, unsigned interval, bool multiple)
{
	file = filename.substr(0,  filename.find_last_of("."));
	file = file + ".rst";
	fileWriteInterval = interval * thinning;
	multipleFiles = multiple;
	writeRestartFile = true;
}


/* setStepsToAdapt (RCPP EXPOSED)
 * Arguments: steps (unsigned)
 * Will set the specified steps to adapt for the run if the value is less than samples * thinning (aka, the number
 * of steps the run will last).The default parameter passed in as -1 uses the full iterations.
*/
//' @name setStepsToAdapt
//' @title setStepsToAdapt
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Set number of iterations (total iterations = samples * thinning) to allow proposal widths to adapt
//' @param steps a postive value
void MCMCAlgorithm::setStepsToAdapt(unsigned steps)
{
	if (steps <= samples * thinning)
		stepsToAdapt = steps;
	else
		my_printError("ERROR: Cannot set steps - value must be smaller than samples times thinning (maxIterations)\n");
}


/* getStepsToAdapt (RCPP EXPOSED)
 * Arguments: None
 * Return the value of stepsToAdapt
*/
//' @name getStepsToAdapt
//' @title getStepsToAdapt
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Return number of iterations (total iterations = samples * thinning) to allow proposal widths to adapt
//' @return number of sample steps to adapt
int MCMCAlgorithm::getStepsToAdapt()
{
	return stepsToAdapt;
}


/* getLogPosteriorTrace (RCPP EXPOSED)
 * Arguments: None
 * Return the log posterior trace.
*/

//' @name getLogPosteriorTrace 
//' @title getLogPosteriorTrace 
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Returns the logPosterior trace
//' @return vector representing logPosterior trace
std::vector<double> MCMCAlgorithm::getLogPosteriorTrace()
{
	return posteriorTrace;
}


/* getLogLikelihoodTrace (RCPP EXPOSED)
 * Arguments: None
 * Return the log likelihoodTrace trace.
*/

//' @name getLogLikelihoodTrace
//' @title getLogLikelihoodTrace
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Returns the logLikelihood trace
//' @return vector representing logLikelihood trace
std::vector<double> MCMCAlgorithm::getLogLikelihoodTrace()
{
	return likelihoodTrace;
}


/* getLogLikelihoodPosteriorMean (RCPP EXPOSED)
 * Arguments: number of samples to use to calculate the posterior.
 * Calculates the mean for the specified number of samples from the end of the trace
 * for the loglikelihood trace.
*/


//' @name getLogPosteriorMean
//' @title getLogPosteriorMean
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Calculate the mean log posterior probability over the last n samples
//' @param samples postive value less than total length of the MCMC trace
//' @return mean logPosterior
double MCMCAlgorithm::getLogPosteriorMean(unsigned _samples)
{
	double posteriorMean = 0.0;
	unsigned traceLength = (unsigned)posteriorTrace.size();

	if (_samples > traceLength)
	{
		my_printError("Warning in MCMCAlgorithm::getLogLikelihoodPosteriorMean throws: Number of anticipated samples");
		my_printError("(%) is greater than the length of the available trace (%). Whole trace is used for posterior estimate! \n", _samples, traceLength);
	}
	unsigned start = traceLength - _samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += posteriorTrace[i];
	}

	return posteriorMean / (double)_samples;
}


/* acf (NOT EXPOSED)
 * Arguments:
 *
*/
std::vector<double> MCMCAlgorithm::acf(std::vector<double>& x, int nrows, int ncols, int lagmax, bool correlation, bool demean)
{
	if (demean)
	{
		double sum = 0.0;
		for (unsigned i = 0u; i < x.size(); i++) sum += x[i];
		double mean = sum / (double)x.size();
		for (unsigned i = 0u; i < x.size(); i++) x[i] = x[i] - mean;
	}

	std::vector<double> acf(lagmax, 1.0);
	int d1 = lagmax + 1, d2 = ncols*d1;

	for (int u = 0; u < ncols; u++)
	{
		for (int v = 0; v < ncols; v++)
		{
			for (int lag = 0; lag <= lagmax; lag++)
			{
				double sum = 0.0; int nu = 0;
				for (int i = 0; i < nrows-lag; i++)
				{
					nu++;
					sum += x[i + lag + nrows*u] * x[i + nrows*v];
				}
				acf[lag + d1*u + d2*v] = sum/(nu + lag);
			}
		}
	}
	if (correlation)
	{
		if (nrows == 1)
		{
			for (int u = 0; u < ncols; u++)
			{
				acf[0 + d1*u + d2*u] = 1.0;
			}
		}
		else
		{
			double *se = new double[ncols]();
			for (int u = 0; u < ncols; u++)
			{
				se[u] = sqrt(acf[0 + d1*u + d2*u]);
			}
			for (int u = 0; u < ncols; u++)
			{
				for (int v = 0; v < ncols; v++)
				{
					for (int lag = 0; lag <= lagmax; lag++) // ensure correlations remain in  [-1,1] :
					{
						double a = acf[lag + d1*u + d2*v] / (se[u]*se[v]);
						acf[lag + d1*u + d2*v] = (a > 1.) ? 1. : ((a < -1.) ? -1. : a);
					}
				}
			}
			// delete [] se; //TODO: add once testable
		}
	}
	return acf;
}


/* solveToeplitzMatrix (NOT EXPOSED)
 * Arguments:
 *
*/
std::vector<std::vector<double>> MCMCAlgorithm::solveToeplitzMatrix(int lr, std::vector<double> r, std::vector<double> g)
{
// TODO switch f from 2d index to 1d index
	//choleskyMatrix[i * numVariates + k]
	//      solves Toeplitz matrix equation toep(r)f=g(1+.)
	//		by Levinson's algorithm
	//      a is a workspace of size lr, the number of equations
	std::vector<double> f(lr*lr, 0.0);
	std::vector<double> var(lr, 0.0);
	std::vector<std::vector<double>> returnVec(2);

	double* a = new double[lr]();

	unsigned l1,l2,k;
	double v, d, q, hold;
	v = r[0];
	d = r[1];
	a[0] = 1.0;
	f[0] = g[1]/v;
	q = f[0]*r[1];
	var[0] = (1 - f[0]*f[0])*r[0];

	if (lr == 1) return returnVec;
	for (unsigned l = 1; l < lr; l++)
	{
		a[l] = -d/v;
		if (l > 1)
		{
			l1 = (l - 2)/2;
			l2 = l1 + 1;
			for (unsigned j = 1; j < l2; j++)
			{
				hold = a[j];
				k = l - j + 1;
				a[j] = a[j] + a[l]*a[k];
				a[k] = a[k] + a[l]*hold;
			}
			if (2*l1 != l - 2) a[l2+1] = a[l2+1]*(1.0 + a[l]);
		}
		v = v + a[l]*d;
		f[l*lr + l] = (g[l+1] - q)/v;
		for (unsigned j = 0; j < (l-1); j++)
		{
			f[l*lr + j] = f[(l-1)*lr +j] + f[l*lr + l]*a[l-j+1];
		}
		//  estimate the innovations variance
		var[l] = var[l-1] * (1 - f[l*lr + l]*f[l*lr + l]);
		if (l == lr) return returnVec;
		d = 0.0;
		q = 0.0;
		for (unsigned i = 0; i < l; i++)
		{
			k = l-i+2;
			d = d + a[i]*r[k];
			q = q + f[l*lr + i]*r[k];
		}
	}
	//delete [] a; //TODO: add once testable

	returnVec[0] = f;
	returnVec[1] = var;
	return returnVec;
}





// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE


//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


/* getSamples (RCPP EXPOSED)
 * Arguments: None
 * Return samples.
*/
//' @name getSamples
//' @title getSamples 
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Return number of samples set for MCMCAlgorithm object
//' @return number of samples used during MCMC
unsigned MCMCAlgorithm::getSamples()
{
    return samples;
}


/* getThinning (RCPP EXPOSED)
 * Arguments: None
 * Return thinning.
*/
//' @name getThinning
//' @title getThinning 
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Return thinning value, which is the number of iterations (total iterations = samples * thinning) not being kept
//' @return thinning value used during MCMC
unsigned MCMCAlgorithm::getThinning()
{
    return thinning;
}


/* getAdaptiveWidth (RCPP EXPOSED)
 * Arguments: None
 * Return adaptiveWidth.
*/
//' @name getAdaptiveWidth
//' @title getAdaptiveWidth
//' @description Return sample adaptiveWidth value, which is the number of samples (not iterations) between adapting parameter proposal widths
//' @return number of sample steps between adapting proposal widths
unsigned MCMCAlgorithm::getAdaptiveWidth()
{
    return adaptiveWidth;
}


/* setSamples (RCPP EXPOSED)
 * Arguments: None
 * Change samples.
*/
//' @name setSamples
//' @title setSamples
//' @description Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Set number of samples set for MCMCAlgorithm object
//' @param _samples postive value

void MCMCAlgorithm::setSamples(unsigned _samples)
{
    samples = _samples;
}


/* setThinning(RCPP EXPOSED)
 * Arguments: None
 * Change thinning.
*/
//' @name setThinning
//' @title setThinning
//' @description  Set thinning value, which is the number of iterations (total iterations = samples * thinning) not being kept
//' @param _thinning postive value

void MCMCAlgorithm::setThinning(unsigned _thinning)
{
    thinning = _thinning;
}


/* setAdaptiveWidth (RCPP EXPOSED)
 * Arguments: None
 * Change adaptiveWidth.
*/
//' @name setAdaptiveWidth
//' @title setAdaptiveWidth
//' @description  Method of MCMC class (access via mcmc$<function name>, where mcmc is an object initialized by initializeMCMCObject). Set sample adaptiveWidth value, which is the number of samples (not iterations) between adapting parameter proposal widths
//' @param _adaptiveWidth postive value

void MCMCAlgorithm::setAdaptiveWidth(unsigned _adaptiveWidth)
{
    adaptiveWidth = _adaptiveWidth;
}


/* setLogPosteriorTrace (RCPP EXPOSED)
 * Arguments: None
 * Change log poster trace.
*/
void MCMCAlgorithm::setLogPosteriorTrace(std::vector<double> _posteriorTrace)
{
    posteriorTrace = _posteriorTrace;
}


/* setLogLikelihoodTrace (RCPP EXPOSED)
 * Arguments: None
 * Change log likelihood trace.
*/
void MCMCAlgorithm::setLogLikelihoodTrace(std::vector<double> _likelihoodTrace)
{
	likelihoodTrace = _likelihoodTrace;
}


//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//

RCPP_EXPOSED_CLASS(Genome)
RCPP_EXPOSED_CLASS(ROCParameter)
RCPP_EXPOSED_CLASS(ROCModel)
RCPP_EXPOSED_CLASS(PANSEModel)
RCPP_EXPOSED_CLASS(PANSEParameter)
RCPP_EXPOSED_CLASS(Model)


RCPP_MODULE(MCMCAlgorithm_mod)
{
	class_<MCMCAlgorithm>( "MCMCAlgorithm" )

		//Constructors & Destructors:
		.constructor("empty constructor")
		.constructor <unsigned, unsigned, unsigned, bool, bool, bool>()



		//MCMC Functions:
		.method("run", &MCMCAlgorithm::run)
		.method("run_PANSE", &MCMCAlgorithm::run_PANSE)
		.method("setEstimateMixtureAssignment", &MCMCAlgorithm::setEstimateMixtureAssignment)
		.method("setRestartFileSettings", &MCMCAlgorithm::setRestartFileSettings)
		.method("getLogPosteriorTrace", &MCMCAlgorithm::getLogPosteriorTrace)
		.method("getLogLikelihoodTrace", &MCMCAlgorithm::getLogLikelihoodTrace)
		.method("getLogPosteriorMean", &MCMCAlgorithm::getLogPosteriorMean)



		//Other Functions:
        .method("getSamples", &MCMCAlgorithm::getSamples)
        .method("getThinning", &MCMCAlgorithm::getThinning)
        .method("getAdaptiveWidth", &MCMCAlgorithm::getAdaptiveWidth)
        .method("setSamples", &MCMCAlgorithm::setSamples)
        .method("setThinning", &MCMCAlgorithm::setThinning)
        .method("setAdaptiveWidth", &MCMCAlgorithm::setAdaptiveWidth)
        .method("setLogPosteriorTrace", &MCMCAlgorithm::setLogPosteriorTrace)
        .method("setLogLikelihoodTrace", &MCMCAlgorithm::setLogLikelihoodTrace)
        .method("setStepsToAdapt", &MCMCAlgorithm::setStepsToAdapt)
        .method("getStepsToAdapt", &MCMCAlgorithm::getStepsToAdapt)
		;



	//MCMC Functions:
	//function("TestACF", &MCMCAlgorithm::acf); //TEST THAT ONLY!
	//function("solveToeplitzMatrix", &MCMCAlgorithm::solveToeplitzMatrix); //TEST THAT ONLY!
}
#endif
