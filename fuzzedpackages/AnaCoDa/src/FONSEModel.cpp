#include "include/FONSE/FONSEModel.h"


//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//


FONSEModel::FONSEModel(bool _withPhi, bool _fix_sEpsilon) : Model()
{
	parameter = 0;
	withPhi = _withPhi;
	fix_sEpsilon = _fix_sEpsilon;

}


FONSEModel::~FONSEModel()
{
	//dtor
	//delete parameter;
}


double FONSEModel::calculateLogLikelihoodRatioPerAA(Gene& gene, std::string grouping, double *mutation, double *selection, double phiValue,double a1_value)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double logLikelihood = 0.0;

	std::vector <unsigned> *positions;
	std::vector <double> codonProb(6, 0);

	//Find the maximum index
	unsigned minIndexVal = 0u;
	for (unsigned i = 1; i < (numCodons - 1); i++)
	{
		if (selection[minIndexVal] > selection[i])
		{
			minIndexVal = i;
		}
	}

	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(grouping, aaStart, aaEnd, false);
	for (unsigned i = aaStart, k = 0; i < aaEnd; i++, k++)
	{
		positions = gene.geneData.getCodonPositions(i);
		for (unsigned j = 0; j < positions->size(); j++)
		{
			calculateLogCodonProbabilityVector(numCodons, positions->at(j), minIndexVal, mutation, selection, phiValue, a1_value, codonProb);
			logLikelihood += codonProb[k];
		}
		//positions->clear();
	}
 	return logLikelihood;
}


double FONSEModel::calculateMutationPrior(std::string grouping, bool proposed)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping, true);
	double mutation[5];

	double priorValue = 0.0;

	unsigned numMutCat = parameter->getNumMutationCategories();
	double mutation_prior_sd = parameter->getMutationPriorStandardDeviation();
	for (unsigned i = 0u; i < numMutCat; i++)
	{
		parameter->getParameterForCategory(i, FONSEParameter::dM, grouping, proposed, mutation);
		for (unsigned k = 0u; k < numCodons; k++)
		{
			priorValue += Parameter::densityNorm(mutation[k], 0.0, mutation_prior_sd, true);
		}
	}
	return priorValue;
}



//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void FONSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;
	std::string curAA;
	double mutation[5];
	double selection[5];

	//SequenceSummary *sequenceSummary = gene.getSequenceSummary(); //currently unused

	// get correct index for everything
	unsigned mutationCategory = parameter->getMutationCategory(k);
	unsigned selectionCategory = parameter->getSelectionCategory(k);
	unsigned expressionCategory = parameter->getSynthesisRateCategory(k);

	double phiValue = parameter->getSynthesisRate(geneIndex, expressionCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, expressionCategory, true);
	double a1_value = getInitiationCost(false);

	SequenceSummary *sequenceSummary = gene.getSequenceSummary();

#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection,curAA) reduction(+:likelihood,likelihood_proposed)
#endif
	for (unsigned i = 0u; i < getGroupListSize(); i++)
	{
		curAA = getGrouping(i);
		if (sequenceSummary->getAACountForAA(i) == 0) continue;
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, curAA, false, mutation);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, curAA, false, selection);
		likelihood += calculateLogLikelihoodRatioPerAA(gene, curAA, mutation, selection, phiValue, a1_value);
		likelihood_proposed += calculateLogLikelihoodRatioPerAA(gene, curAA, mutation, selection, phiValue_proposed, a1_value);
	}

	unsigned mixture = getMixtureAssignment(geneIndex);
	mixture = getSynthesisRateCategory(mixture);
	double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(mixture, false);
	double mPhi = (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5); // X * 0.5 = X / 2
	//double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(selectionCategory, false);
	double logPhiProbability = Parameter::densityLogNorm(phiValue, mPhi, stdDevSynthesisRate, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, mPhi, stdDevSynthesisRate, true);
	if (withPhi)
	{
		for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
		{
			double obsPhi = gene.getObservedSynthesisRate(i);
			if (obsPhi > -1.0)
			{
				double logObsPhi = std::log(obsPhi);
				logPhiProbability += Parameter::densityNorm(logObsPhi, std::log(phiValue) + getNoiseOffset(i), getObservedSynthesisNoise(i), true);
				logPhiProbability_proposed += Parameter::densityNorm(logObsPhi, std::log(phiValue_proposed) + getNoiseOffset(i), getObservedSynthesisNoise(i), true);
			}
		}
	}
	double currentLogPosterior = (likelihood + logPhiProbability);
	double proposedLogPosterior = (likelihood_proposed + logPhiProbability_proposed);
	logProbabilityRatio[0] = (proposedLogPosterior - currentLogPosterior) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogPosterior - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogPosterior - std::log(phiValue);
	logProbabilityRatio[3] = currentLogPosterior;
	logProbabilityRatio[4] = proposedLogPosterior;
	logProbabilityRatio[5] = likelihood;
	logProbabilityRatio[6] = likelihood_proposed;
}


void FONSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures)
{
	unsigned numGenes = genome.getGenomeSize();
	//int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;
	double posterior_proposed,posterior;

	double mutation[5];
	double selection[5];
	double mutation_proposed[5];
	double selection_proposed[5];

	std::string curAA;

	Gene *gene;
	SequenceSummary *sequenceSummary;
	unsigned aaIndex = SequenceSummary::AAToAAIndex(grouping);
	double a1_value = getInitiationCost(false);
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, curAA, gene, sequenceSummary) reduction(+:likelihood,likelihood_proposed)
#endif
	for (unsigned i = 0u; i < numGenes; i++)
	{
		gene = &genome.getGene(i);
		sequenceSummary = gene->getSequenceSummary();
		if (sequenceSummary->getAACountForAA(aaIndex) == 0) continue;

		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned mutationCategory = parameter->getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get phi value, calculate likelihood conditional on phi
		double phiValue = parameter->getSynthesisRate(i, expressionCategory, false);

		// get current mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, grouping, false, mutation);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, grouping, false, selection);

		// get proposed mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, grouping, true, mutation_proposed);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, grouping, true, selection_proposed);
		likelihood += calculateLogLikelihoodRatioPerAA(*gene, grouping, mutation, selection, phiValue,a1_value);
		likelihood_proposed += calculateLogLikelihoodRatioPerAA(*gene, grouping, mutation_proposed, selection_proposed, phiValue,a1_value);
	}

	posterior_proposed = likelihood_proposed + calculateMutationPrior(grouping, true);
	posterior = likelihood + calculateMutationPrior(grouping, false);
	logAcceptanceRatioForAllMixtures[0] = (posterior_proposed - posterior);
	logAcceptanceRatioForAllMixtures[1] = likelihood;
	logAcceptanceRatioForAllMixtures[2] = likelihood_proposed;
	logAcceptanceRatioForAllMixtures[3] = posterior;
	logAcceptanceRatioForAllMixtures[4] = posterior_proposed;

}


void FONSEModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{
	double mutation[5];
	double selection[5];
	std::string curAA;
	Gene *gene;

	double lpr_sphi = 0.0;
	double lpr_a1 = 0.0;
	double lpr = 0.0;
	if (withPhi)
	{
		// one for each noiseOffset, one for stdDevSynthesisRate, one for initiation_cost a1
		logProbabilityRatio.resize(getNumPhiGroupings() + 2);
	}
	else
	{
		logProbabilityRatio.resize(2);
	}

	double a1_current = getInitiationCost(false);
	double a1_proposed = getInitiationCost(true);

	unsigned selectionCategory = getNumSynthesisRateCategories();
	std::vector<double> currentStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> currentMphi(selectionCategory, 0.0);
	std::vector<double> proposedStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> proposedMphi(selectionCategory, 0.0);
	for (unsigned i = 0u; i < selectionCategory; i++)
	{
		currentStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, false);
		currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) * 0.5);
		proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
		proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) * 0.5);
		// take the Jacobian into account for the non-linear transformation from logN to N distribution
		lpr_sphi -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
	}

	lpr_a1 -= (std::log(a1_current) - std::log(a1_proposed));

#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene,mutation, selection, curAA) reduction(+:lpr_sphi,lpr_a1)
#endif
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		gene = &genome.getGene(i);
		unsigned mixtureElement = getMixtureAssignment(i);
		unsigned expressionCategory = getSynthesisRateCategory(mixtureElement);
		double phi = getSynthesisRate(i, expressionCategory, false);
		lpr_sphi += Parameter::densityLogNorm(phi, proposedMphi[expressionCategory], proposedStdDevSynthesisRate[expressionCategory], true)
			   - Parameter::densityLogNorm(phi, currentMphi[expressionCategory], currentStdDevSynthesisRate[expressionCategory], true);

		unsigned mutationCategory = parameter->getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter->getSelectionCategory(mixtureElement);

		for (unsigned j = 0u; j < getGroupListSize(); j++)
		{
			curAA = getGrouping(j);

			parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, curAA, false, mutation);
			parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, curAA, false, selection);

			lpr_a1 += (calculateLogLikelihoodRatioPerAA(*gene, curAA, mutation, selection, phi, a1_proposed) 
				- calculateLogLikelihoodRatioPerAA(*gene, curAA, mutation, selection, phi, a1_current));
		}
	}
	logProbabilityRatio[0] = lpr_sphi;
	logProbabilityRatio[1] = lpr_a1;
	if (withPhi)
	{
		for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
		{
			
			lpr = 0.0;
			double noiseOffset = getNoiseOffset(i, false);
			double noiseOffset_proposed = getNoiseOffset(i, true);
			double observedSynthesisNoise = getObservedSynthesisNoise(i);
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
			for (unsigned j = 0u; j < genome.getGenomeSize(); j++)
			{
				unsigned mixtureAssignment = getMixtureAssignment(j);
				mixtureAssignment = getSynthesisRateCategory(mixtureAssignment);
				double logPhi = std::log(getSynthesisRate(j, mixtureAssignment, false));
				double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
				if (obsPhi > -1.0)
				{
					double logObsPhi = std::log(obsPhi);
					double proposed = Parameter::densityNorm(logObsPhi, logPhi + noiseOffset_proposed, observedSynthesisNoise, true);
					double current = Parameter::densityNorm(logObsPhi, logPhi + noiseOffset, observedSynthesisNoise, true);
					lpr += proposed - current;
				}
			}
			logProbabilityRatio[i+2] = lpr;
		}
	}

}





//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void FONSEModel::initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
	parameter->initAllTraces(samples, num_genes, estimateSynthesisRate);
}


void FONSEModel::writeRestartFile(std::string filename)
{
	return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double FONSEModel::getCategoryProbability(unsigned i)
{
	return parameter->getCategoryProbability(i);
}


unsigned FONSEModel::getMutationCategory(unsigned mixture)
{
	return parameter->getMutationCategory(mixture);
}


unsigned FONSEModel::getSelectionCategory(unsigned mixture)
{
	return parameter->getSelectionCategory(mixture);
}


unsigned FONSEModel::getSynthesisRateCategory(unsigned mixture)
{
	return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> FONSEModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
	return parameter->getMixtureElementsOfSelectionCategory(k);
}


double FONSEModel::getInitiationCost(bool proposed)
{
	return parameter->getInitiationCost(proposed);
}



//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned FONSEModel::getGroupListSize()
{
	return parameter->getGroupListSize();
}


std::string FONSEModel::getGrouping(unsigned index)
{
	return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double FONSEModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double FONSEModel::getCurrentStdDevSynthesisRateProposalWidth()
{
	return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}


void FONSEModel::updateStdDevSynthesisRate()
{
	parameter->updateStdDevSynthesisRate();
}


double FONSEModel::getCurrentInitiationCostProposalWidth()
{
	return parameter->getCurrentInitiationCostProposalWidth();
}

void FONSEModel::updateInitiationCost()
{
	parameter->updateInitiationCost();
}





//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double FONSEModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
	return parameter->getSynthesisRate(index, mixture, proposed);
}


void FONSEModel::updateSynthesisRate(unsigned i, unsigned k)
{
	parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned FONSEModel::getLastIteration()
{
	return parameter->getLastIteration();
}


void FONSEModel::setLastIteration(unsigned iteration)
{
	parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void FONSEModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
	parameter->updateStdDevSynthesisRateTrace(sample);
}

void FONSEModel::updateInitiationCostParameterTrace(unsigned sample)
{
	parameter->updateInitiationCostParameterTrace(sample);
}

void FONSEModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
	parameter->updateSynthesisRateTrace(sample, i);
}


void FONSEModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
	parameter->updateMixtureAssignmentTrace(sample, i);
}


void FONSEModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
	parameter->updateMixtureProbabilitiesTrace(sample);
}


void FONSEModel::updateCodonSpecificParameterTrace(unsigned sample, std::string grouping)
{
	parameter->updateCodonSpecificParameterTrace(sample, grouping);
}


void FONSEModel::updateHyperParameterTraces(unsigned sample)
{
	updateStdDevSynthesisRateTrace(sample);
	updateInitiationCostParameterTrace(sample);
	if (withPhi)
	{
		updateNoiseOffsetTrace(sample);
    	updateObservedSynthesisNoiseTrace(sample);
    }
}


void FONSEModel::updateTracesWithInitialValues(Genome & genome)
{
	std::vector <std::string> groupList = parameter->getGroupList();

	for (unsigned i = 0; i < genome.getGenomeSize(); i++)
	{
		parameter->updateSynthesisRateTrace(0, i);
		parameter->updateMixtureAssignmentTrace(0, i);
	}

	for (unsigned i = 0; i < groupList.size(); i++)
	{
		parameter->updateCodonSpecificParameterTrace(0, getGrouping(i));
	}
}





//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


void FONSEModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void FONSEModel::adaptInitiationCostProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptInitiationCostProposalWidth(adaptiveWidth, adapt);
}

void FONSEModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void FONSEModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
	parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void FONSEModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
	adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
	adaptInitiationCostProposalWidth(adaptiveWidth,adapt);
	if (withPhi)
		adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void FONSEModel::proposeCodonSpecificParameter()
{
	parameter->proposeCodonSpecificParameter();
}


void FONSEModel::proposeHyperParameters()
{
	parameter->proposeHyperParameters();
	if (withPhi)
	{
		parameter->proposeNoiseOffset();
	}
}


void FONSEModel::proposeSynthesisRateLevels()
{
	parameter->proposeSynthesisRateLevels();
}


unsigned FONSEModel::getNumPhiGroupings()
{
	return parameter->getNumObservedPhiSets();
}


unsigned FONSEModel::getMixtureAssignment(unsigned index)
{
	return parameter->getMixtureAssignment(index);
}


unsigned FONSEModel::getNumMixtureElements()
{
	return parameter->getNumMixtureElements();
}


unsigned FONSEModel::getNumSynthesisRateCategories()
{
	return parameter->getNumSynthesisRateCategories();
}


void FONSEModel::setNumPhiGroupings(unsigned value)
{
	parameter->setNumObservedPhiSets(value);
}


void FONSEModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
	parameter->setMixtureAssignment(i, catOfGene);
}


void FONSEModel::setCategoryProbability(unsigned mixture, double value)
{
	parameter->setCategoryProbability(mixture, value);
}


void FONSEModel::updateCodonSpecificParameter(std::string grouping)
{
	parameter->updateCodonSpecificParameter(grouping);
}

void FONSEModel::completeUpdateCodonSpecificParameter()
{
    parameter->completeUpdateCodonSpecificParameter();
}

//Noise offset functions

double FONSEModel::getNoiseOffset(unsigned index, bool proposed)
{
	return parameter->getNoiseOffset(index, proposed);
}


double FONSEModel::getObservedSynthesisNoise(unsigned index)
{
	return parameter->getObservedSynthesisNoise(index);
}


double FONSEModel::getCurrentNoiseOffsetProposalWidth(unsigned index)
{
	return parameter->getCurrentNoiseOffsetProposalWidth(index);
}


void FONSEModel::updateNoiseOffset(unsigned index)
{
	parameter->updateNoiseOffset(index);
}


void FONSEModel::updateNoiseOffsetTrace(unsigned sample)
{
	parameter->updateNoiseOffsetTraces(sample);
}


void FONSEModel::updateObservedSynthesisNoiseTrace(unsigned sample)
{
	parameter->updateObservedSynthesisNoiseTraces(sample);
}


void FONSEModel::adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}



void FONSEModel::updateGibbsSampledHyperParameters(Genome &genome)
{
  // estimate s_epsilon by sampling from a gamma distribution and transforming it into an inverse gamma sample
	
	if (withPhi)
	{
		if(!fix_sEpsilon)
		{
			double shape = ((double)genome.getGenomeSize() - 1.0) / 2.0;
			for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
			{
				double rate = 0.0; //Prior on s_epsilon goes here?
				unsigned mixtureAssignment;
				double noiseOffset = getNoiseOffset(i);
				for (unsigned j = 0; j < genome.getGenomeSize(); j++)
				{
					mixtureAssignment = parameter->getMixtureAssignment(j);
					double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
					if (obsPhi > -1.0)
					{
						double sum = std::log(obsPhi) - noiseOffset - std::log(parameter->getSynthesisRate(j, mixtureAssignment, false));
						rate += (sum * sum);
					}else{
						// missing observation.
						shape -= 0.5;
						//Reduce shape because initial estimate assumes there are no missing observations
					}
				}
				rate /= 2.0;
				double rand = parameter->randGamma(shape, rate);

				// Below the gamma sample is transformed into an inverse gamma sample
				// According to Gilchrist et al (2015) Supporting Materials p. S6
				// The sample 1/T is supposed to be equal to $s_\epsilon^2$.
				double sepsilon = std::sqrt(1.0/rand);
				parameter->setObservedSynthesisNoise(i, sepsilon);
			}
		}
	}
}



void FONSEModel::updateAllHyperParameter()
{
	updateStdDevSynthesisRate();
	updateInitiationCost();
	if (withPhi)
	{
		for (unsigned i =0; i < parameter->getNumObservedPhiSets(); i++)
		{
			updateNoiseOffset(i);
		}
	}
}

void FONSEModel::updateHyperParameter(unsigned hp)
{
	// NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
	if (hp == 0)
	{
		updateStdDevSynthesisRate();
	}
	else if (hp == 1)
	{
		updateInitiationCost();
	}		
	else if (hp > 1 and withPhi)
	{	
		//subtract off 2 because the first two parameters withh be the updateStdDevSynthesisRate
		updateNoiseOffset(hp - 2);
	}
}


void FONSEModel::simulateGenome(Genome & genome)
{
	unsigned codonIndex;
	std::string curAA;

	std::string tmpDesc = "Simulated Gene";

	for (unsigned geneIndex = 0; geneIndex < genome.getGenomeSize(); geneIndex++) //loop over all genes in the genome
	{
		Gene gene = genome.getGene(geneIndex);
		SequenceSummary sequenceSummary = gene.geneData;
		std::string tmpSeq = "ATG"; //Always will have the start amino acid
		unsigned mixtureElement = getMixtureAssignment(geneIndex);
		unsigned mutationCategory = getMutationCategory(mixtureElement);
		unsigned selectionCategory = getSelectionCategory(mixtureElement);
		unsigned synthesisRateCategory = getSynthesisRateCategory(mixtureElement);
		double phi = getSynthesisRate(geneIndex, synthesisRateCategory, false);
		double a1_value = getInitiationCost(false);
		std::string geneSeq = gene.getSequence();
		for (unsigned position = 1; position < (geneSeq.size() / 3); position++)
		{
			std::string codon = geneSeq.substr((position * 3), 3);
			curAA = SequenceSummary::codonToAA(codon);

			//TODO: Throw an error here instead
			if (curAA == "X") {
				if (position < (geneSeq.size() / 3) - 1) my_print("Warning: Internal stop codon found in gene % at position %. Ignoring and moving on.\n", gene.getId(), position);
				continue;
			}

			unsigned numCodons = SequenceSummary::GetNumCodonsForAA(curAA);

			double* codonProb = new double[numCodons](); //size the arrays to the proper size based on # of codons.
			double* mutation = new double[numCodons - 1]();
			double* selection = new double[numCodons - 1]();


			if (curAA == "M" || curAA == "W")
			{
				codonProb[0] = 1;
			}
			else
			{
				getParameterForCategory(mutationCategory, FONSEParameter::dM, curAA, false, mutation);
				getParameterForCategory(selectionCategory, FONSEParameter::dOmega, curAA, false, selection);
				calculateCodonProbabilityVector(numCodons, position, mutation, selection, phi, a1_value,codonProb);
			}


			codonIndex = Parameter::randMultinom(codonProb, numCodons);
			unsigned aaStart, aaEnd;
			SequenceSummary::AAToCodonRange(curAA, aaStart, aaEnd, false);  //need the first spot in the array where the codons for curAA are
			codon = sequenceSummary.indexToCodon(aaStart + codonIndex);//get the correct codon based off codonIndex
			tmpSeq += codon;
		}
		std::string codon = sequenceSummary.indexToCodon((unsigned)Parameter::randUnif(61.0, 64.0)); //randomly choose a stop codon, from range 61-63
		tmpSeq += codon;
		Gene simulatedGene(tmpSeq, tmpDesc, gene.getId());
		genome.addGene(simulatedGene, true);
	}
}


void FONSEModel::printHyperParameters()
{
	for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
	{
		my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, getStdDevSynthesisRate(i));
	}
	my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());
	my_print("\t current initiation cost a_1 estimate: %\n",getInitiationCost(false));
	my_print("\t current initiation cost a_1 proposal width: %\n", getCurrentInitiationCostProposalWidth());

}


/* getParameter (RCPP EXPOSED)
* Arguments: None
*
* Returns the FONSEParameter of the model.
*/
FONSEParameter* FONSEModel::getParameter()
{
	return parameter;
}

void FONSEModel::setParameter(FONSEParameter &_parameter)
{
	parameter = &_parameter;
}


double FONSEModel::calculateAllPriors()
{
	double priorRatio = 0.0;
	unsigned size = getGroupListSize();

	for (unsigned i = 0; i < size; i++)
	{
		std::string grouping = getGrouping(i);
		priorRatio += calculateMutationPrior(grouping, false);
	}

	// add more priors if necessary.

	return priorRatio;
}

//Calculates the log probability of each codon for an amino acid and puts them in a vector.
void FONSEModel::calculateLogCodonProbabilityVector(unsigned numCodons, unsigned position, unsigned minIndexValue,
												 double *mutation, double *selection, double phi, double a1_value, std::vector <double> &codonProb)
{
	double denominator;

	/* log(c_i) = \Delta M - (\phi * \beta(i) * \Delta \omega),                 *
	 * where \beta(i) = a_1 + (i * a_2)                                         *
	 *                                                                          *
	 * Right now a_1 and a_2 are set to 4.0. However, we are planning on making *
	 * them hyperparameters in the future, since they are constant for the      *
	 * entire genome.                                                           */

	 // If the min(selection) is less than zero than we have to adjust the reference codon.
	 // If the reference codon is the min value (0) then we do not have to adjust the reference codon.
	 // This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
	 // causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexValue] < 0.0)
	{
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = -(mutation[i] - mutation[minIndexValue]) - (phi * (a1_value + (4.0 * position)) * (selection[i] - selection[minIndexValue]));
			denominator += std::exp(codonProb[i]);
		}
		//Alphabetically, the last codon is the reference codon.
		codonProb[numCodons - 1] = (mutation[minIndexValue]) + (phi * (a1_value + (4.0 * position)) * selection[minIndexValue]);
		denominator += std::exp(codonProb[numCodons - 1]);
	}
	else
	{
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = -(mutation[i]) - (phi * (a1_value + (4.0 * position)) * selection[i]);
			denominator += std::exp(codonProb[i]);
		}
		//Again, the last codon is the reference codon
		codonProb[numCodons - 1] = 0.0;
	}

	//Here we take the log of the denominator (the summation term) so that we can finish calculating
	//the log probabilities simple by subtracting the log of the denominator from each element.
	denominator = std::log(denominator);
	for (unsigned i = 0; i < numCodons; i++)
	{
		codonProb[i] = codonProb[i] - denominator;
	}
}


//Since the simulateGenome function utilizes the codon probability vector, but doesn't deal with the log values,
//this function simply returns the vector with each codon's probability.
void FONSEModel::calculateCodonProbabilityVector(unsigned numCodons, unsigned position,
													double *mutation, double *selection, double phi, double a1_value, double codonProb[])
{
	double denominator;
	unsigned minIndexValue = 0u;
	for (unsigned i = 1; i < (numCodons - 1); i++)
	{
		if (selection[minIndexValue] < selection[i])
		{
			minIndexValue = i;
		}
	}

	/* c_i = exp[\Delta M - (\phi * \beta(i) * \Delta \omega)],                 *
	* where \beta(i) = a_1 + (i * a_2)                                         *
	*                                                                          *
	* Right now a_1 and a_2 are set to 4.0. However, we are planning on making *
	* them hyperparameters in the future, since they are constant for the      *
	* entire genome.
	* Note codon position, indexing starts at 0,thus we don't need (i-1)
	* \Delta \omega includes q Ne terms                                        */


	// If the min(selection) is less than zero than we have to adjust the reference codon.
	// If the reference codon is the min value (0) then we do not have to adjust the reference codon.
	// This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
	// causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexValue] < 0.0) {
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp(-(mutation[i] - mutation[minIndexValue]) - (phi * (a1_value + (4.0 * position)) * (selection[i] - selection[minIndexValue])));
			denominator += codonProb[i];
		}
		//Alphabetically, the last codon is the reference codon.
		codonProb[numCodons - 1] = std::exp((mutation[minIndexValue]) + (phi * (a1_value + (4.0 * position)) * selection[minIndexValue]));
		denominator += codonProb[numCodons - 1];
	}
	else
	{
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp(-(mutation[i]) - (phi * (a1_value + (4.0 * position)) * selection[i]));
			denominator += codonProb[i];
		}
		//Again, the last codon is the reference codon
		codonProb[numCodons - 1] = 1.0;
	}

	//As is found in ROCModel.cpp, multiplication is a faster operation than division so we
	//save time here by dividing once and then multiplying numCodons times instead of dividing
	//numCodons times.
	denominator = 1 / denominator;
	for (unsigned i = 0; i < numCodons; i++)
	{
		codonProb[i] *= denominator;
	}
}



void FONSEModel::getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal, double* returnValue)
{
	parameter->getParameterForCategory(category, param, aa, proposal, returnValue);
}








// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//


#ifndef STANDALONE


std::vector<double> FONSEModel::CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi, double a1_value, unsigned position)
{
	unsigned numCodons = mutation.size() + 1;
	double* _mutation = &mutation[0];
	double* _selection = &selection[0];
	double* codonProb = new double[numCodons]();
	calculateCodonProbabilityVector(numCodons, position, _mutation, _selection, phi, a1_value, codonProb);
	std::vector<double> returnVector(codonProb, codonProb + numCodons);
	return returnVector;
}


#endif
