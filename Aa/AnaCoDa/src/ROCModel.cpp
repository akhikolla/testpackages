#include "include/ROC/ROCModel.h"


//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//


ROCModel::ROCModel(bool _withPhi, bool _fix_sEpsilon) : Model()
{
	parameter = 0;
	withPhi = _withPhi;
	fix_sEpsilon = _fix_sEpsilon;
}


ROCModel::~ROCModel()
{
	//dtor
	//delete parameter;
}

double ROCModel::calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue)
{
	double logLikelihood = 0.0;
	// calculate codon probabilities
	double logCodonProbabilities[6];
	calculateLogCodonProbabilityVector(numCodons, mutation, selection, phiValue, logCodonProbabilities);

	// calculate likelihood for current AA for this combination of selection and mutation category
	for (unsigned i = 0; i < numCodons; i++)
	{
		if (codonCount[i] == 0) continue;
		logLikelihood += logCodonProbabilities[i] * codonCount[i];
	}
	return logLikelihood;
}


double ROCModel::calculateMutationPrior(std::string grouping, bool proposed)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping, true);
	double mutation[5],mutation_mean[5],mutation_sd[5];

	double priorValue = 0.0;

	unsigned numMutCat = parameter->getNumMutationCategories();
	for (unsigned i = 0u; i < numMutCat; i++)
	{
		parameter->getParameterForCategory(i, ROCParameter::dM, grouping, proposed, mutation);
		parameter->getMutationPriorStandardDeviationForCategoryForGroup(i,grouping,mutation_sd);
		parameter->getMutationPriorMeanForCategoryForGroup(i,grouping,mutation_mean);
		for (unsigned k = 0u; k < numCodons; k++)
			priorValue += Parameter::densityNorm(mutation[k], mutation_mean[k], mutation_sd[k], true);
	}
	return priorValue;
}


void ROCModel::obtainCodonCount(SequenceSummary *sequenceSummary, std::string curAA, int codonCount[])
{
	unsigned aaStart, aaEnd;
	SequenceSummary::AAToCodonRange(curAA, aaStart, aaEnd, false);
	// get codon counts for AA
	unsigned j = 0u;
	for (unsigned i = aaStart; i < aaEnd; i++, j++)
		codonCount[j] = sequenceSummary->getCodonCountForCodon(i);
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void ROCModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;

	SequenceSummary *sequenceSummary = gene.getSequenceSummary();

	// get correct index for everything
	unsigned mutationCategory = parameter->getMutationCategory(k);
	unsigned selectionCategory = parameter->getSelectionCategory(k);
	unsigned expressionCategory = parameter->getSynthesisRateCategory(k);

	double phiValue = parameter->getSynthesisRate(geneIndex, expressionCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, expressionCategory, true);

	double mutation[5];
	double selection[5];
	int codonCount[6];
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, codonCount) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
	for (unsigned i = 0u; i < getGroupListSize(); i++)
	{
		std::string curAA = getGrouping(i);

		// skip amino acids which do not occur in current gene. Avoid useless calculations and multiplying by 0
		if (sequenceSummary->getAACountForAA(i) == 0) continue;

		// get number of codons for AA (total number not parameter->count)
		unsigned numCodons = sequenceSummary->GetNumCodonsForAA(curAA);
		// get mutation and selection parameter->for gene
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
		// get codon occurrence in sequence
		obtainCodonCount(sequenceSummary, curAA, codonCount);

		logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue_proposed);
	}
	unsigned mixture = getMixtureAssignment(geneIndex);
	mixture = getSynthesisRateCategory(mixture);
	double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(mixture, false);
	double mPhi = (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5); // X * 0.5 = X / 2
	double logPhiProbability = Parameter::densityLogNorm(phiValue, mPhi, stdDevSynthesisRate, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, mPhi, stdDevSynthesisRate, true);

	// TODO: make this work for more than one phi value, or for genes that don't have phi values
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

	double currentLogPosterior = (logLikelihood + logPhiProbability);
	double proposedLogPosterior = (logLikelihood_proposed + logPhiProbability_proposed);

	//TODO: Can't see where logProbabilityRatio[0], which is acceptance ratio with reverse jump, is used. Consider deleting. 
	logProbabilityRatio[0] = (proposedLogPosterior - currentLogPosterior) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogPosterior - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogPosterior - std::log(phiValue);
	logProbabilityRatio[3] = currentLogPosterior;
	logProbabilityRatio[4] = proposedLogPosterior;
	logProbabilityRatio[5] = logLikelihood;
	logProbabilityRatio[6] = logLikelihood_proposed;

}


void ROCModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures)
{
	int numGenes = genome.getGenomeSize();
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	//my_print("Current grouping: %\n",grouping);
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;
	double posterior, posterior_proposed;
	double mutation[5];
	double selection[5];
	double mutation_proposed[5];
	double selection_proposed[5];

	int codonCount[6];
	Gene *gene;
	SequenceSummary *sequenceSummary;
	unsigned aaIndex = SequenceSummary::AAToAAIndex(grouping);
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, codonCount, gene, sequenceSummary) reduction(+:likelihood,likelihood_proposed)
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
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, grouping, false, mutation);
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, grouping, false, selection);
		// get proposed mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, grouping, true, mutation_proposed);
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, grouping, true, selection_proposed);
		obtainCodonCount(sequenceSummary, grouping, codonCount);
		likelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation_proposed, selection_proposed, phiValue);
	}
	bool dm_fixed = parameter -> isDMFixed();
	if (!dm_fixed)
	{
		posterior_proposed = likelihood_proposed + calculateMutationPrior(grouping, true);
		posterior = likelihood + calculateMutationPrior(grouping, false);
	}
	else
    {
		posterior_proposed = likelihood_proposed;
		posterior = likelihood;
	}
	logAcceptanceRatioForAllMixtures[0] = (posterior_proposed - posterior);
	logAcceptanceRatioForAllMixtures[1] = likelihood;
	logAcceptanceRatioForAllMixtures[2] = likelihood_proposed;
	logAcceptanceRatioForAllMixtures[3] = posterior;
	logAcceptanceRatioForAllMixtures[4] = posterior_proposed;
}


void ROCModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> &logProbabilityRatio)
{
	double lpr = 0.0;
	unsigned selectionCategory = getNumSynthesisRateCategories();
	std::vector<double> currentStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> currentMphi(selectionCategory, 0.0);
	std::vector<double> proposedStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> proposedMphi(selectionCategory, 0.0);

	//Calculating reverse jump probabilities due to asymmetry of logNormal
	for (unsigned i = 0u; i < selectionCategory; i++)
	{
		currentStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, false);
		currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) * 0.5);
		proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
		proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) * 0.5);
		// take the Jacobian into account for the non-linear transformation from logN to N distribution
		lpr -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
	}

	if (withPhi)
	{
		// one for each noiseOffset, and one for stdDevSynthesisRate
		logProbabilityRatio.resize(getNumPhiGroupings() + 1);
	}
	else
		logProbabilityRatio.resize(1);

#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		unsigned mixture = getMixtureAssignment(i);
		mixture = getSynthesisRateCategory(mixture);
		double phi = getSynthesisRate(i, mixture, false);
		if (!std::isfinite(phi))
			my_printError("Error: Phi value for gene % is not finite (%)!", i, phi);

		lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedStdDevSynthesisRate[mixture], true)
			   - Parameter::densityLogNorm(phi, currentMphi[mixture], currentStdDevSynthesisRate[mixture], true);
	}

	// TODO: USE CONSTANTS INSTEAD OF 0
	logProbabilityRatio[0] = lpr;

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
			logProbabilityRatio[i+1] = lpr;
		}
	}
}





//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void ROCModel::initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
	parameter -> initAllTraces(samples, num_genes,estimateSynthesisRate);
}



void ROCModel::writeRestartFile(std::string filename)
{
	return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double ROCModel::getCategoryProbability(unsigned i)
{
	return parameter->getCategoryProbability(i);
}


unsigned ROCModel::getMutationCategory(unsigned mixture)
{
	return parameter ->getMutationCategory(mixture);
}


unsigned ROCModel::getSelectionCategory(unsigned mixture)
{
	return parameter ->getSelectionCategory(mixture);
}


unsigned ROCModel::getSynthesisRateCategory(unsigned mixture)
{
	return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> ROCModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
	return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned ROCModel::getGroupListSize()
{
	return parameter->getGroupListSize();
}


std::string ROCModel::getGrouping(unsigned index)
{
	return parameter -> getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double ROCModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double ROCModel::getCurrentStdDevSynthesisRateProposalWidth()
{
	return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}


void ROCModel::updateStdDevSynthesisRate()
{
	parameter->updateStdDevSynthesisRate();
}





//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double ROCModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
	return parameter->getSynthesisRate(index, mixture, proposed);
}


void ROCModel::updateSynthesisRate(unsigned i, unsigned k)
{
	parameter->updateSynthesisRate(i,k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned ROCModel::getLastIteration()
{
	return parameter->getLastIteration();
}


void ROCModel::setLastIteration(unsigned iteration)
{
	parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void ROCModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
	parameter->updateStdDevSynthesisRateTrace(sample);
}


void ROCModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
	parameter->updateSynthesisRateTrace(sample, i);
}


void ROCModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
	parameter->updateMixtureAssignmentTrace(sample, i);
}


void ROCModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
	parameter->updateMixtureProbabilitiesTrace(sample);
}


void ROCModel::updateCodonSpecificParameterTrace(unsigned sample, std::string grouping)
{
	parameter->updateCodonSpecificParameterTrace(sample,grouping);
}


void ROCModel::updateHyperParameterTraces(unsigned sample)
{
	updateStdDevSynthesisRateTrace(sample);
	if (withPhi)
	{
		updateNoiseOffsetTrace(sample);
		updateObservedSynthesisNoiseTrace(sample);
	}
}


void ROCModel::updateTracesWithInitialValues(Genome &genome)
{
	std::vector <std::string> groupList = parameter->getGroupList();

	// for (unsigned i = 0; i < genome.getGenomeSize(); i++)
	// {
	// 	parameter->updateSynthesisRateTrace(0, i);
	// 	parameter->updateMixtureAssignmentTrace(0, i);
	// }

	for (unsigned i = 0; i < groupList.size(); i++)
		parameter->updateCodonSpecificParameterTrace(0, getGrouping(i));
}





//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


void ROCModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void ROCModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void ROCModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
	parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void ROCModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
	adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
	if (withPhi)
		adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}


//Noise offset functions

double ROCModel::getNoiseOffset(unsigned index, bool proposed)
{
	return parameter->getNoiseOffset(index, proposed);
}


double ROCModel::getObservedSynthesisNoise(unsigned index)
{
	return parameter->getObservedSynthesisNoise(index);
}


double ROCModel::getCurrentNoiseOffsetProposalWidth(unsigned index)
{
	return parameter->getCurrentNoiseOffsetProposalWidth(index);
}


void ROCModel::updateNoiseOffset(unsigned index)
{
	parameter->updateNoiseOffset(index);
}


void ROCModel::updateNoiseOffsetTrace(unsigned sample)
{
	parameter->updateNoiseOffsetTraces(sample);
}


void ROCModel::updateObservedSynthesisNoiseTrace(unsigned sample)
{
	parameter->updateObservedSynthesisNoiseTraces(sample);
}


void ROCModel::adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}



void ROCModel::updateGibbsSampledHyperParameters(Genome &genome)
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




//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void ROCModel::proposeCodonSpecificParameter()
{
	parameter->proposeCodonSpecificParameter();
}


void ROCModel::proposeHyperParameters()
{
	parameter->proposeStdDevSynthesisRate();
	if (withPhi)
		parameter->proposeNoiseOffset();
}


void ROCModel::proposeSynthesisRateLevels()
{
	parameter->proposeSynthesisRateLevels();
}


unsigned ROCModel::getNumPhiGroupings()
{
	return parameter->getNumObservedPhiSets();
}


unsigned ROCModel::getMixtureAssignment(unsigned index)
{
	return parameter->getMixtureAssignment(index);
}


unsigned ROCModel::getNumMixtureElements()
{
	return parameter->getNumMixtureElements();
}


unsigned ROCModel::getNumSynthesisRateCategories()
{
	return parameter->getNumSynthesisRateCategories();
}


void ROCModel::setNumPhiGroupings(unsigned value)
{
	parameter->setNumObservedPhiSets(value);
}


void ROCModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
	parameter->setMixtureAssignment(i, catOfGene);
}


void ROCModel::setCategoryProbability(unsigned mixture, double value)
{
	parameter->setCategoryProbability(mixture, value);
}


void ROCModel::updateCodonSpecificParameter(std::string grouping)
{
	parameter->updateCodonSpecificParameter(grouping);
}

void ROCModel::completeUpdateCodonSpecificParameter()
{
    parameter->completeUpdateCodonSpecificParameter();
}


void ROCModel::updateAllHyperParameter()
{
	updateStdDevSynthesisRate();
	for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
	{
		updateNoiseOffset(i);
	}
}


void ROCModel::updateHyperParameter(unsigned hp)
{
	if (hp == 0)
		updateStdDevSynthesisRate();
	else if (hp > 0 && withPhi)
		updateNoiseOffset(hp - 1);
}


void ROCModel::simulateGenome(Genome &genome)
{
	unsigned codonIndex;

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

		std::string geneSeq = gene.getSequence();
		for (unsigned position = 1; position < (geneSeq.size() / 3); position++)
		{
			std::string codon = geneSeq.substr((position * 3), 3);
			std::string aa = SequenceSummary::codonToAA(codon);

			if (aa == "X") {
				if (position < (geneSeq.size() / 3) - 1) my_print("Warning: Internal stop codon found in gene % at position %. Ignoring and moving on.\n", gene.getId(), position);
				continue;
			}

			unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa);

			double* codonProb = new double[numCodons](); //size the arrays to the proper size based on # of codons.
			double* mutation = new double[numCodons - 1]();
			double* selection = new double[numCodons - 1]();

			if (aa == "M" || aa == "W")
				codonProb[0] = 1;
			else
			{
				getParameterForCategory(mutationCategory, ROCParameter::dM, aa, false, mutation);
				getParameterForCategory(selectionCategory, ROCParameter::dEta, aa, false, selection);
				calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb);
			}
			codonIndex = Parameter::randMultinom(codonProb, numCodons);
			unsigned aaStart, aaEnd;
			SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, false); //need the first spot in the array where the codons for curAA are
			codon = sequenceSummary.indexToCodon(aaStart + codonIndex);//get the correct codon based off codonIndex
			tmpSeq += codon;
		}
		std::string codon =	sequenceSummary.indexToCodon((unsigned)Parameter::randUnif(61.0, 64.0)); //randomly choose a stop codon, from range 61-63
		tmpSeq += codon;
		Gene simulatedGene(tmpSeq, tmpDesc, gene.getId());
		genome.addGene(simulatedGene, true);
	}
}


void ROCModel::printHyperParameters()
{
	for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
		my_print("\t current stdDevSynthesisRate estimate for selection category %: %\n", i, getStdDevSynthesisRate(i, false));

	my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());

	if (withPhi)
	{
		my_print("\t current noiseOffset estimates:");

		for (unsigned i = 0; i < getNumPhiGroupings(); i++)
			my_print(" %", getNoiseOffset(i, false));

		my_print("\n\t current noiseOffset proposal widths:");

		for (unsigned i = 0; i < getNumPhiGroupings(); i++)
			my_print(" %", getCurrentNoiseOffsetProposalWidth(i));

		my_print("\n\t current observedSynthesisNoise estimates:");

		for (unsigned i = 0; i < getNumPhiGroupings(); i++)
			my_print(" %", getObservedSynthesisNoise(i));

		my_print("\n");
	}
}

/* getParameter (RCPP EXPOSED)
* Arguments: None
*
* Returns the ROCParameter of the model.
*/
ROCParameter ROCModel::getParameter()
{
	return *parameter;
}

void ROCModel::setParameter(ROCParameter &_parameter)
{
	parameter = &_parameter;
}


double ROCModel::calculateAllPriors()
{
	double prior = 0.0;
	unsigned size = getGroupListSize();

	for (unsigned i = 0; i < size; i++)
	{
		std::string grouping = getGrouping(i);
		prior += calculateMutationPrior(grouping, false);
	}

	// add more priors if necessary.

	return prior;
}


void ROCModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
{
	// calculate numerator and denominator for codon probabilities
	unsigned minIndexVal = 0u;
	double denominator;
	for (unsigned i = 1u; i < (numCodons - 1); i++)
	{
		if (selection[minIndexVal] > selection[i])
			minIndexVal = i;
	}

	// if the min(selection) is less than zero than we have to adjust the reference codon.
	// if the reference codon is the min value (0) then we do not have to adjust the reference codon.
	// This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
	// causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexVal] < 0.0)
	{
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp( -(mutation[i] - mutation[minIndexVal]) - ((selection[i] - selection[minIndexVal]) * phi) );
			//codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
			denominator += codonProb[i];
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = std::exp(mutation[minIndexVal] + selection[minIndexVal] * phi);
		denominator += codonProb[numCodons - 1];
	}
	else
	{
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
			denominator += codonProb[i];
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = 1.0;
	}

	denominator = 1 / denominator; // Multiplication is faster than division: replace multiple divisions below by one up here.
	// normalize codon probabilities
	for (unsigned i = 0u; i < numCodons; i++)
		codonProb[i] = codonProb[i] * denominator; // denominator is 1/denominator. see above
}


//Function that calculates the log of each codon's probability instead of just the probability.
void ROCModel::calculateLogCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
{
	// calculate numerator and denominator for codon probabilities
	unsigned minIndexVal = 0u;
	double denominator;
	for (unsigned i = 1u; i < (numCodons - 1); i++)
	{
		if (selection[minIndexVal] > selection[i])
			minIndexVal = i;
	}

	// if the min(selection) is less than zero than we have to adjust the reference codon.
	// if the reference codon is the min value (0) then we do not have to adjust the reference codon.
	// This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
	// causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexVal] < 0.0)
	{
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = -(mutation[i] - mutation[minIndexVal]) - ((selection[i] - selection[minIndexVal]) * phi);
			//codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
			denominator += std::exp(codonProb[i]);
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = mutation[minIndexVal] + selection[minIndexVal] * phi;
		denominator += std::exp(codonProb[numCodons - 1]);
	}
	else
	{
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++)
		{
			codonProb[i] = -mutation[i] - (selection[i] * phi);
			denominator += std::exp(codonProb[i]);
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = 0.0;
	}

	denominator = std::log(denominator);
	for (unsigned i = 0u; i < numCodons; i++)
		codonProb[i] = codonProb[i] - denominator;
}


void ROCModel::getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal, double* returnValue)
{
	parameter -> getParameterForCategory(category, param, aa, proposal, returnValue);
}












// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//


#ifndef STANDALONE


std::vector<double> ROCModel::CalculateProbabilitiesForCodons(std::vector<double> mutation,
		std::vector<double> selection, double phi)
{
	unsigned numCodons = mutation.size() + 1;
	double* _mutation = &mutation[0];
	double* _selection = &selection[0];
	double* codonProb = new double[numCodons]();
	calculateCodonProbabilityVector(numCodons, _mutation, _selection, phi, codonProb);
	std::vector<double> returnVector(codonProb, codonProb + numCodons);
	return returnVector;
}

#endif
