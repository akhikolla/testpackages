#include "include/PANSE/PANSEModel.h"

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//

PANSEModel::PANSEModel(unsigned _RFPCountColumn, bool _withPhi, bool _fix_sEpsilon) : Model()
{
    parameter = NULL;
    RFPCountColumn = _RFPCountColumn - 1;
    withPhi = _withPhi;
    fix_sEpsilon = _fix_sEpsilon;
    //ctor
}


PANSEModel::~PANSEModel()
{
    //dtor
    //TODO: call Parent's deconstructor
    //delete parameter;
}


void PANSEModel::calculateZ(std::string grouping,Genome& genome,std::vector<double> &Z,std::string param)
{

    double currAlpha,currLambda,currNSERate,propAlpha,propLambda,propNSERate;
    double currZ = 0;
    double propZ = 0;
    Gene* gene;
    prob_successful.resize(getGroupListSize(),1000);
#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambda,currNSERate,propAlpha,propLambda,propNSERate) reduction(+:currZ,propZ)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        double prop_prob_successful = 1000;
        gene = &genome.getGene(i);
        unsigned codonIndex;
        std::string codon;

        unsigned mixtureElement = parameter->getMixtureAssignment(i);
        // how is the mixture element defined. Which categories make it up
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        

        double phi = parameter->getSynthesisRate(i, synthesisRateCategory, false);
        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();

        double propSigma = 0;
        double currSigma = 0;
        double currZ_gene = 0;
        double propZ_gene = 0;
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            codonIndex = positions[positionIndex];
            codon = gene->geneData.indexToCodon(codonIndex);
            
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            currLambda = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
            
            if(codon == grouping)
            {
                if (param == "Elongation")
                {
                    propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true);
                    propLambda = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, true);
                    propZ_gene += std::exp(propSigma) * propAlpha/propLambda;
                    if (prop_prob_successful > 500)
                    {
                        prop_prob_successful = elongationProbabilityLog(propAlpha, propLambda,1/currNSERate);
                        if (prop_prob_successful > 0.0)
                        {
                            prop_prob_successful = 0.0;
                        }
                    }
                    
                }
                else
                {
                    propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
                    propZ_gene += std::exp(propSigma) * currAlpha/currLambda;
                    if (prop_prob_successful > 500)
                    {
                        prop_prob_successful = elongationProbabilityLog(currAlpha, currLambda,1/propNSERate);
                        if (prop_prob_successful > 0.0)
                        {
                            prop_prob_successful = 0.0;
                        }
                    }
                }
                propSigma = propSigma + prop_prob_successful;
           }
           else
           {
                propZ_gene += std::exp(propSigma) * currAlpha/currLambda;
                if (prob_successful[codonIndex] > 500)
                {
                    prob_successful[codonIndex] = elongationProbabilityLog(currAlpha, currLambda,1/currNSERate);
                    if (prob_successful[codonIndex] > 0.0)
                    {
                        prob_successful[codonIndex] = 0.0;
                    }
                }
                propSigma = propSigma + prob_successful[codonIndex];  
            }

            currZ_gene += std::exp(currSigma) * currAlpha/currLambda;
            if (prob_successful[codonIndex] > 500)
            {
                prob_successful[codonIndex] = elongationProbabilityLog(currAlpha, currLambda,1/currNSERate);
                if (prob_successful[codonIndex] > 0.0)
                {
                    prob_successful[codonIndex] = 0.0;
                }
            }

            currSigma = currSigma + prob_successful[codonIndex];
        }
        currZ += phi * currZ_gene;
        propZ += phi * propZ_gene;
    }
    prob_successful.clear();
    Z[0] = currZ;
    Z[1] = propZ;
}




double PANSEModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
        unsigned currRFPObserved, double phiValue, double prevSigma, double lgamma_currAlpha, double log_currLambdaPrime, double log_phi,double lgamma_rfp_alpha)
{
    double term1 = lgamma_rfp_alpha - lgamma_currAlpha;//std::lgamma(currAlpha);
    double term2 = log_phi + std::log(prevSigma) - std::log(currLambdaPrime + (phiValue * prevSigma));
    double term3 = log_currLambdaPrime - std::log(currLambdaPrime + (phiValue * prevSigma));

    term2 *= currRFPObserved;
    term3 *= currAlpha;


    double rv = term1 + term2 + term3;
    return rv;
}


//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void PANSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
    double currAlpha,currLambda,currNSERate;
    std::string codon;
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    unsigned codonIndex,positionalRFPCount;

    std::vector <unsigned> positions = gene.geneData.getPositionCodonID();
    std::vector <unsigned> rfpCounts = gene.geneData.getRFPCount(0);

    std::vector<double> local_lgamma_currentAlpha(getGroupListSize(), -1000);
    std::vector<double> local_log_currentLambda(getGroupListSize(), 1000);
    std::vector<double> local_prob_successful(getGroupListSize(),1000);


    unsigned alphaCategory = parameter->getMutationCategory(k);
    unsigned lambdaCategory = parameter->getSelectionCategory(k);
    unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);
    unsigned mixtureElement = parameter->getMixtureAssignment(geneIndex);
    
    unsigned Y = parameter->getTotalRFPCount();
    double U = getPartitionFunction(mixtureElement, false)/Y;

    double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
    double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);
    double logPhi = std::log(phiValue); 
    double logPhi_proposed = std::log(phiValue_proposed);
    
    double currSigma = 0;
    

    for (unsigned positionIndex = 0; positionIndex < positions.size();positionIndex++)
    {
        positionalRFPCount = rfpCounts[positionIndex];
        codonIndex = positions[positionIndex];
        codon = gene.geneData.indexToCodon(codonIndex);
        
        currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
        currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        
        if (local_lgamma_currentAlpha[codonIndex] < -5)
        {
            local_lgamma_currentAlpha[codonIndex] = std::lgamma(currAlpha);
        }
        if (local_log_currentLambda[codonIndex] > 500)
        {
            local_log_currentLambda[codonIndex] = std::log(currLambda) + std::log(U);
        }

        double currLgammaRFPAlpha = std::lgamma(currAlpha+positionalRFPCount);
        logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount, phiValue, std::exp(currSigma), 
                                local_lgamma_currentAlpha[codonIndex],local_log_currentLambda[codonIndex], logPhi, currLgammaRFPAlpha);
        logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount, phiValue_proposed, std::exp(currSigma), 
                                local_lgamma_currentAlpha[codonIndex],local_log_currentLambda[codonIndex], logPhi_proposed, currLgammaRFPAlpha);
        if (local_prob_successful[codonIndex] > 500)
        {
            local_prob_successful[codonIndex] = elongationProbabilityLog(currAlpha, currLambda,1/currNSERate);
            if (prob_successful[codonIndex] > 0.0)
            {
                local_prob_successful[codonIndex] = 0.0;
            }
        }
        currSigma = currSigma + local_prob_successful[codonIndex];
        if (currSigma > 0)
        {
            my_print("WARNING: sigma for gene % is greater than 1.\n",geneIndex);
        }
    }


    //Double check math here
    double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(synthesisRateCategory, false);
    double logPhiProbability = Parameter::densityLogNorm(phiValue, (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5), stdDevSynthesisRate, true);
    double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5), stdDevSynthesisRate, true);
	
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

	logProbabilityRatio[0] = (proposedLogPosterior - currentLogPosterior) - (std::log(phiValue) - std::log(phiValue_proposed));//Is recalulcated in MCMC
	logProbabilityRatio[1] = currentLogPosterior - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogPosterior - std::log(phiValue);
	logProbabilityRatio[3] = currentLogPosterior;
	logProbabilityRatio[4] = proposedLogPosterior;
	logProbabilityRatio[5] = logLikelihood;
	logProbabilityRatio[6] = logLikelihood_proposed;
}



void PANSEModel::fillMatrices(Genome& genome)
{
    unsigned n = getNumMixtureElements();
    //These vectors should be visible to all threads. If one thread tries to overwrite another one, should only result in it replacing the same value
    for (unsigned j = 0; j < n; j++)
    {
        std::vector<double> tmp = std::vector<double> (getGroupListSize(), -1000);
        std::vector<double> tmp_2 = std::vector<double> (getGroupListSize(),1000);
        lgamma_currentAlpha.push_back(tmp);
        log_currentLambda.push_back(tmp_2);
    }
    lgamma_rfp_alpha.resize(50);
    for (unsigned k = 0; k < 50;k++)
    {
        lgamma_rfp_alpha[k].resize(n);
        for (unsigned j = 0; j < n; j++)
        {
            std::vector<double> tmp = std::vector<double> (getGroupListSize(), -1000);
            lgamma_rfp_alpha[k][j]=tmp;
        }
    }
    prob_successful.resize(getGroupListSize(),1000);

}

void PANSEModel::clearMatrices()
{
    lgamma_currentAlpha.clear();
    log_currentLambda.clear();
    lgamma_rfp_alpha.clear();
    prob_successful.clear();
}

void PANSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures)
{
    
}


void PANSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures,std::string param)
{
    std::vector<std::string> groups = parameter -> getGroupList();
   
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    double propAlpha, propLambda, propNSERate;
    
    double currAlpha, currLambda, currNSERate;
    double currAdjustmentTerm = 0;
    double propAdjustmentTerm = 0;
    Gene *gene;
    //unsigned index = SequenceSummary::codonToIndex(grouping);
    std::vector<double> codonSigmas;
    unsigned n = getNumMixtureElements();
    unsigned Y = genome.getSumRFP();

    fillMatrices(genome);

#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambda,currNSERate,propAlpha,propLambda,propNSERate) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        double prop_prob_successful = 1000;
        gene = &genome.getGene(i);
        unsigned positionalRFPCount;
        unsigned codonIndex;
        std::string codon;
        double currLgammaRFPAlpha;
        
        unsigned mixtureElement = parameter->getMixtureAssignment(i);
        
        double U = getPartitionFunction(mixtureElement, false)/Y;

        // how is the mixture element defined. Which categories make it up
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        

        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
        double log_phi = std::log(phiValue);
        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <unsigned> rfpCounts = gene->geneData.getRFPCount(0);
    
        double propSigma = 0;
        double currSigma = 0;
       
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            positionalRFPCount = rfpCounts[positionIndex];
            codonIndex = positions[positionIndex];
            codon = gene->geneData.indexToCodon(codonIndex);
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);


            if (lgamma_currentAlpha[alphaCategory][codonIndex] < -5)
            {
                lgamma_currentAlpha[alphaCategory][codonIndex] = std::lgamma(currAlpha);
            }
           
            if (log_currentLambda[lambdaCategory][codonIndex] > 500)
            {
                log_currentLambda[lambdaCategory][codonIndex] = std::log(currLambda)+ std::log(U);
            }
            
            if (positionalRFPCount < 50)
            {
                if (lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] < -5)
                {
                    lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] = std::lgamma(currAlpha + positionalRFPCount);
                }
                currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
            }
            else
            {
                currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
            }
           


            if(codon == grouping)
            {
                if (param == "Elongation")
                {
                    propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true);
                    propLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, true);
                    logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambda * U, positionalRFPCount,
                                    phiValue,std::exp(propSigma),std::lgamma(propAlpha),std::log(propLambda) + std::log(U),log_phi,std::lgamma(propAlpha+positionalRFPCount));
                    if (prop_prob_successful > 500)
                    {
                        prop_prob_successful = elongationProbabilityLog(propAlpha, propLambda,1/currNSERate);
                        if (prop_prob_successful > 0.0)
                        {
                            prop_prob_successful = 0.0;
                        }
               
                    }
                }
                else
                {
                    propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
                    logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount,
                                    phiValue,std::exp(propSigma),std::lgamma(currAlpha),std::log(currLambda) + std::log(U),log_phi,std::lgamma(currAlpha+positionalRFPCount));
                    if (prop_prob_successful > 500)
                    {
                        prop_prob_successful = elongationProbabilityLog(currAlpha, currLambda,1/propNSERate);
                        if (prop_prob_successful > 0.0)
                        {
                            prop_prob_successful = 0.0;
                        }
                    }
                }
                propSigma = propSigma + prop_prob_successful;
           }
           else
           {
                logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount,
                                  phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],std::log(currLambda) + std::log(U),log_phi,currLgammaRFPAlpha);
                if (prob_successful[codonIndex] > 500)
                {
                    prob_successful[codonIndex] = elongationProbabilityLog(currAlpha, currLambda,1/currNSERate);
                    if (prob_successful[codonIndex] > 0.0)
                    {
                        prob_successful[codonIndex] = 0.0;
                    }
                }
                propSigma = propSigma + prob_successful[codonIndex];  
            }
            


            logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda*U, positionalRFPCount,
                                    phiValue, std::exp(currSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[lambdaCategory][codonIndex],log_phi,currLgammaRFPAlpha);  
            
            if (prob_successful[codonIndex] > 500)
            {
               
                prob_successful[codonIndex] = elongationProbabilityLog(currAlpha, currLambda,1/currNSERate);
                if (prob_successful[codonIndex] > 0.0)
                {
                    prob_successful[codonIndex] = 0.0;
                }
            }

            currSigma = currSigma + prob_successful[codonIndex];
        }
    }
    
    for (unsigned j = 0; j < n; j++)
    {
        unsigned alphaCategory = parameter->getMutationCategory(j);
        unsigned lambdaCategory = parameter->getSelectionCategory(j);
        if (param == "Elongation")
        {
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, false);
            propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, true);
            propLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, true);
            if (std::isnan(logLikelihood_proposed))
            {
                my_print("WARNING: proposed logLikelihood for % is NaN\n",grouping);
                my_print("\tProposed alpha: % \n",getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, true));
                my_print("\tProposed lambda: %\n",getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, true));
            }
            currAdjustmentTerm += std::log(currAlpha) + std::log(currLambda);
            propAdjustmentTerm += std::log(propAlpha) + std::log(propLambda);
        }
        else
        {
            currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, false);
            propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, true);
            if (std::isnan(logLikelihood_proposed))
            {
                my_print("WARNING: proposed logLikelihood for % is NaN\n",grouping);
                my_print("\tProposed NSE Rate: %\n",getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, true));

            }
            currAdjustmentTerm += std::log(currNSERate);
            propAdjustmentTerm += std::log(propNSERate);
        }
    }
    
    //Should never accept parameters that give NaN, so just check proposed parameters
    
    logAcceptanceRatioForAllMixtures[0] = logLikelihood_proposed - logLikelihood - (currAdjustmentTerm - propAdjustmentTerm);
	logAcceptanceRatioForAllMixtures[3] = logLikelihood;
	logAcceptanceRatioForAllMixtures[4] = logLikelihood_proposed;
	logAcceptanceRatioForAllMixtures[1] = logLikelihood;
	logAcceptanceRatioForAllMixtures[2] = logLikelihood_proposed;

    clearMatrices();
}


void PANSEModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{

    double lpr = 0.0; // this variable is only needed because OpenMP doesn't allow variables in reduction clause to be reference

    unsigned selectionCategory = getNumSynthesisRateCategories();
    std::vector<double> currentStdDevSynthesisRate(selectionCategory, 0.0);
    std::vector<double> currentMphi(selectionCategory, 0.0);
    std::vector<double> proposedStdDevSynthesisRate(selectionCategory, 0.0);
    std::vector<double> proposedMphi(selectionCategory, 0.0);


    for (unsigned i = 0u; i < selectionCategory; i++)
    {
        currentStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, false);
        currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) / 2);
        proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
        proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) / 2);
        // take the Jacobian into account for the non-linear transformation from logN to N distribution
        lpr -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
     }


    if (withPhi)
    {
        // one for each noiseOffset, and one for stdDevSynthesisRate
        logProbabilityRatio.resize(getNumPhiGroupings() + 2);
    }
    else
        logProbabilityRatio.resize(2);



#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        unsigned mixture = getMixtureAssignment(i);
        mixture = getSynthesisRateCategory(mixture);
        double phi = getSynthesisRate(i, mixture, false);
        lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedStdDevSynthesisRate[mixture], true) -
            Parameter::densityLogNorm(phi, currentMphi[mixture], currentStdDevSynthesisRate[mixture], true);
    }

    logProbabilityRatio[0] = lpr;

    Gene *gene;
    double currAlpha, currLambda, currNSERate;
    double logLikelihood = 0;
    double logLikelihood_proposed = 0;
    unsigned n = getNumMixtureElements();
    lpr = 0.0;

    unsigned Y = parameter -> getTotalRFPCount();
    fillMatrices(genome);
    
    
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambda,currNSERate) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        //std::vector<double> prop_prob_successful(getGroupListSize(),1000);
        unsigned positionalRFPCount;
        unsigned codonIndex;
        std::string codon;
        double currLgammaRFPAlpha;
        
        gene = &genome.getGene(i);
        unsigned mixtureElement = parameter->getMixtureAssignment(i);

        // how is the mixture element defined. Which categories make it up
        double currU = getPartitionFunction(mixtureElement, false)/Y;
        double propU = getPartitionFunction(mixtureElement, true)/Y;

        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        // get non codon specific values, calculate likelihood conditional on these
        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
        double log_phi = std::log(phiValue);
        
        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <unsigned> rfpCounts = gene->geneData.getRFPCount(/*RFPCountColumn*/ 0);
  
        double currSigma = 0;
        double propSigma = 0;

        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            positionalRFPCount = rfpCounts[positionIndex];
            codonIndex = positions[positionIndex];
            codon = gene->geneData.indexToCodon(codonIndex);
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);          
            

            if (lgamma_currentAlpha[alphaCategory][codonIndex] < -5)
            {
                lgamma_currentAlpha[alphaCategory][codonIndex] = std::lgamma(currAlpha);
            }
            
            if (log_currentLambda[lambdaCategory][codonIndex] > 500)
            {
                log_currentLambda[lambdaCategory][codonIndex] = std::log(currLambda)+ std::log(currU);
            }

             //Check if lgamma(alpha + rfp_count) already counted for this codon 
            if (positionalRFPCount < 50)
            {
                if (lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] < -5)
                {
                    lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] = std::lgamma(currAlpha + positionalRFPCount);
                }
                currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
            }
            else
            {
                currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
            }


            logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, (currLambda*propU), positionalRFPCount,
                                   phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],(std::log(currLambda) + std::log(propU)),log_phi,currLgammaRFPAlpha);
            logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda*currU, positionalRFPCount,
                                   phiValue,std::exp(currSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[lambdaCategory][codonIndex],log_phi,currLgammaRFPAlpha);   

            if (prob_successful[codonIndex] > 500)
            {
                prob_successful[codonIndex] = elongationProbabilityLog(currAlpha, currLambda,1/currNSERate);
                if (prob_successful[codonIndex] > 0.0)
                {
                    prob_successful[codonIndex] = 0.0;
                }
            }

            currSigma = currSigma + prob_successful[codonIndex];
            propSigma = propSigma + prob_successful[codonIndex];
        }
          
    }
   

    //double currZ = getPartitionFunction(0, false);
    //double propZ = getPartitionFunction(0, true);
    for (unsigned j = 0; j < n; j++)
    {
        lpr -= (std::log(getPartitionFunction(j, false)) - std::log(getPartitionFunction(j, true)));
    }
    lpr += logLikelihood_proposed - logLikelihood; 
    logProbabilityRatio[1] = lpr;
    clearMatrices();
    

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


void PANSEModel::initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
    parameter->initAllTraces(samples, num_genes, estimateSynthesisRate);
}


void PANSEModel::writeRestartFile(std::string filename)
{
    return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double PANSEModel::getCategoryProbability(unsigned i)
{
    return parameter->getCategoryProbability(i);
}


unsigned PANSEModel::getMutationCategory(unsigned mixture)
{
    return parameter->getMutationCategory(mixture);
}


unsigned PANSEModel::getSelectionCategory(unsigned mixture)
{
    return parameter->getSelectionCategory(mixture);
}


unsigned PANSEModel::getSynthesisRateCategory(unsigned mixture)
{
    return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> PANSEModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
    return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned PANSEModel::getGroupListSize()
{
    return parameter->getGroupListSize();
}


std::string PANSEModel::getGrouping(unsigned index)
{
    return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double PANSEModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
    return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double PANSEModel::getCurrentStdDevSynthesisRateProposalWidth()
{
    return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}



void PANSEModel::updateStdDevSynthesisRate()
{
    parameter->updateStdDevSynthesisRate();
}


//---------------------------------------------------//
//----------- partitionFunction Functions -----------//
//---------------------------------------------------//


double PANSEModel::getPartitionFunction(unsigned selectionCategory, bool proposed)
{
    return parameter->getPartitionFunction(selectionCategory, proposed);
}


double PANSEModel::getCurrentPartitionFunctionProposalWidth()
{
    return parameter->getCurrentPartitionFunctionProposalWidth();
}



void PANSEModel::updatePartitionFunction()
{
    parameter->updatePartitionFunction();
}



//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double PANSEModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
    return parameter->getSynthesisRate(index, mixture, proposed);
}


void PANSEModel::updateSynthesisRate(unsigned i, unsigned k)
{
    parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned PANSEModel::getLastIteration()
{
    return parameter->getLastIteration();
}


void PANSEModel::setLastIteration(unsigned iteration)
{
    parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void PANSEModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
    parameter->updateStdDevSynthesisRateTrace(sample);
}


void PANSEModel::updatePartitionFunctionTrace(unsigned sample)
{
    parameter->updatePartitionFunctionTrace(sample);
}


void PANSEModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
    parameter->updateSynthesisRateTrace(sample, i);
}


void PANSEModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
    parameter->updateMixtureAssignmentTrace(sample, i);
}


void PANSEModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
    parameter->updateMixtureProbabilitiesTrace(sample);
}


void PANSEModel::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
    parameter->updateCodonSpecificParameterTrace(sample, codon);
}


void PANSEModel::updateHyperParameterTraces(unsigned sample)
{
    updateStdDevSynthesisRateTrace(sample);
    updatePartitionFunctionTrace(sample);
    if (withPhi)
    {
        updateNoiseOffsetTrace(sample);
        updateObservedSynthesisNoiseTrace(sample);
    }
}


void PANSEModel::updateTracesWithInitialValues(Genome & genome)
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


void PANSEModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PANSEModel::adaptPartitionFunctionProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptPartitionFunctionProposalWidth(adaptiveWidth, adapt);
}


void PANSEModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PANSEModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
    parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void PANSEModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
    adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
    adaptPartitionFunctionProposalWidth(adaptiveWidth, adapt);
    if (withPhi)
        adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}



//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void PANSEModel::proposeCodonSpecificParameter()
{
    parameter->proposeCodonSpecificParameter();
}


void PANSEModel::proposeHyperParameters()
{
    parameter->proposeStdDevSynthesisRate();
    parameter->proposePartitionFunction();
    if (withPhi)
    {
        parameter->proposeNoiseOffset();
    }
}


void PANSEModel::proposeSynthesisRateLevels()
{
    parameter->proposeSynthesisRateLevels();
}


unsigned PANSEModel::getNumPhiGroupings()
{
    return parameter->getNumObservedPhiSets();
}


unsigned PANSEModel::getMixtureAssignment(unsigned index)
{
    return parameter->getMixtureAssignment(index);
}


unsigned PANSEModel::getNumMixtureElements()
{
    return parameter->getNumMixtureElements();
}


unsigned PANSEModel::getNumSynthesisRateCategories()
{
    return parameter->getNumSynthesisRateCategories();
}


void PANSEModel::setNumPhiGroupings(unsigned value)
{
    parameter->setNumObservedPhiSets(value);
}


void PANSEModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
    parameter->setMixtureAssignment(i, catOfGene);
}


void PANSEModel::setCategoryProbability(unsigned mixture, double value)
{
    parameter->setCategoryProbability(mixture, value);
}


void PANSEModel::updateCodonSpecificParameter(std::string codon)
{
    my_print("This should not be called in PANSEModel.\n");
    //parameter->updateCodonSpecificParameter(codon);
}


void PANSEModel::updateCodonSpecificParameter(std::string codon,std::string param)
{
    parameter->updateCodonSpecificParameter(codon,param);
}


void PANSEModel::completeUpdateCodonSpecificParameter()
{
    parameter->completeUpdateCodonSpecificParameter();
}

//Noise offset functions

double PANSEModel::getNoiseOffset(unsigned index, bool proposed)
{
    return parameter->getNoiseOffset(index, proposed);
}


double PANSEModel::getObservedSynthesisNoise(unsigned index)
{
    return parameter->getObservedSynthesisNoise(index);
}


double PANSEModel::getCurrentNoiseOffsetProposalWidth(unsigned index)
{
    return parameter->getCurrentNoiseOffsetProposalWidth(index);
}


void PANSEModel::updateNoiseOffset(unsigned index)
{
    parameter->updateNoiseOffset(index);
}


void PANSEModel::updateNoiseOffsetTrace(unsigned sample)
{
    parameter->updateNoiseOffsetTraces(sample);
}


void PANSEModel::updateObservedSynthesisNoiseTrace(unsigned sample)
{
    parameter->updateObservedSynthesisNoiseTraces(sample);
}


void PANSEModel::adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}



void PANSEModel::updateGibbsSampledHyperParameters(Genome &genome)
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


void PANSEModel::updateAllHyperParameter()
{
    updateStdDevSynthesisRate();
    updatePartitionFunction();
    if (withPhi)
    {
        for (unsigned i =0; i < parameter->getNumObservedPhiSets(); i++)
        {
            updateNoiseOffset(i);
        }
    }
}


void PANSEModel::updateHyperParameter(unsigned hp)
{
    // NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
    if (hp == 0)
    {
        updateStdDevSynthesisRate();
    }
    else if (hp == 1)
    {
        updatePartitionFunction();
    }       
    else if (hp > 1 and withPhi)
    {   
        //subtract off 2 because the first two parameters withh be the updateStdDevSynthesisRate
        updateNoiseOffset(hp - 2);
    }
}

void PANSEModel::simulateGenome(Genome &genome)
{
    unsigned Y = genome.getSumRFP();
    double Z = 0;
    std::vector<std::vector<double>> wait_times;
    wait_times.resize(genome.getGenomeSize());

    for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
    {

        unsigned mixtureElement = getMixtureAssignment(geneIndex);
        Gene gene = genome.getGene(geneIndex);
        double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
        SequenceSummary sequence = gene.geneData;
        Gene tmpGene = gene;
        std::vector <unsigned> positions = sequence.getPositionCodonID();
        wait_times[geneIndex].resize(positions.size());
        std::vector <unsigned> rfpCount;
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaCategory = parameter->getSelectionCategory(mixtureElement);

        //double U = getPartitionFunction(mixtureElement, false)/Y;
        double sigma =  1.0;
        double v;
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            unsigned codonIndex = positions[positionIndex];
            std::string codon = sequence.indexToCodon(codonIndex);
            if (positionIndex == 0 && codon != "ATG")
            {
                my_print("Gene % does not start with ATG\n",gene.getId());
            }
            if (codon == "TAG" || codon == "TGA" || codon == "TAA")
            {
                my_print("Stop codon being used during simulations\n");
            }
            double alpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            double lambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            double NSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
            v = 1.0 / NSERate;
#ifndef STANDALONE
            RNGScope scope;
            NumericVector wt(1);
            wt = rgamma(1, alpha,1/(lambda));
            wait_times[geneIndex][positionIndex] = wt[0];
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
#else
            std::gamma_distribution<double> GDistribution(alpha, 1.0/(lambda));
            wait_times[geneIndex][positionIndex] = GDistribution(Parameter::generator);
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
#endif
            Z += phi * wait_times[geneIndex][positionIndex] * sigma;
        }
    }
    double U = Z/Y;

    for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
    {

        unsigned mixtureElement = getMixtureAssignment(geneIndex);
        Gene gene = genome.getGene(geneIndex);
        double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
        SequenceSummary sequence = gene.geneData;
        Gene tmpGene = gene;
        std::vector <unsigned> positions = sequence.getPositionCodonID();
        wait_times[geneIndex].resize(positions.size());
        std::vector <unsigned> rfpCount;
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaCategory = parameter->getSelectionCategory(mixtureElement);

        double sigma =  1.0;
        double v;
        
        
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            unsigned codonIndex = positions[positionIndex];
            std::string codon = sequence.indexToCodon(codonIndex);
            if (positionIndex == 0 && codon != "ATG")
            {
                my_print("Gene % does not start with ATG\n",gene.getId());
            }
            if (codon == "TAG" || codon == "TGA" || codon == "TAA")
            {
                my_print("Stop codon being used during simulations\n");
            }
            //double alpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            //double lambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            double NSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
            v = 1.0 / NSERate;
#ifndef STANDALONE
            RNGScope scope;
            NumericVector ribo_count(1);
            ribo_count = rpois(1, wait_times[geneIndex][positionIndex] * (1.0/U) * phi * sigma);
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
            rfpCount.push_back(ribo_count[0]);
#else
            std::poisson_distribution<unsigned> PDistribution(phi * wait_times[geneIndex][positionIndex] * (1.0/U) * sigma);
            unsigned simulatedValue = PDistribution(Parameter::generator);
            sigma *= (v/(tmp + v));
            rfpCount.push_back(simulatedValue);
#endif
        }

        tmpGene.geneData.setRFPCount(rfpCount, RFPCountColumn);
        genome.addGene(tmpGene, true);
    }
}


void PANSEModel::printHyperParameters()
{
    for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
    {
        my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, parameter->getStdDevSynthesisRate(i));
        my_print("partition function posterior estimate for selection category %: %\n", i, parameter->getPartitionFunction(i,false));
    }
    my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());
}


/* getParameter (RCPP EXPOSED)
 * Arguments: None
 *
 * Returns the PANSEParameter of the model.
 */
PANSEParameter* PANSEModel::getParameter()
{
    return parameter;
}


void PANSEModel::setParameter(PANSEParameter &_parameter)
{
    parameter = &_parameter;
}


double PANSEModel::calculateAllPriors()
{
    return 0.0; //TODO(Cedric): implement me, see ROCModel
}


double PANSEModel::getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal)
{
    return parameter->getParameterForCategory(category, param, codon, proposal);
}

//Continued fractions helper function for upper incomplete gamma
double PANSEModel::UpperIncompleteGammaHelper(double s, double x)
{
    double rv;
    int i;

    rv = 10000.0 / x;
    for(i = 10000; i > 0; i--){
        if(i % 2 == 0) rv = (double) ((int) (i / 2)) / (x + rv);
        if(i % 2 != 0) rv = ((double) (((int) (i / 2) + 1) - s)) / (1 + rv);
    }

    return x + rv;
}

// double PANSEModel::UpperIncompleteGammaHelper(double s, double x)
// {
//     double summation = 1.0;
//     double num = 1.0;
//     unsigned count = 1;
//     for (unsigned i = 0; i < 10; i++)
//     {
//         num = num * (s - count);
//         summation = summation + num/std::pow(x,count);
//         count += 1;
//     }
//     return summation;
// }



//Upper incomplete gamma function
double PANSEModel::UpperIncompleteGamma(double s, double x)
{
    double d, rv;

    rv = pow(x, s) * exp(0 - x);

    d = UpperIncompleteGammaHelper(s, x);

    return rv/d;

}

//log of upper incomplete gamma function
double PANSEModel::UpperIncompleteGammaLog(double s, double x)
{
    double rv, d;

    rv = s * std::log(x) - x;
    d = std::log(UpperIncompleteGammaHelper(s, x));

    return rv - d;

}


//Calculation of the probability of elongation at current codon
double PANSEModel::elongationProbability(double currAlpha, double currLambda, double currNSE)
{
    return std::pow(currLambda * currNSE, currAlpha) * std::exp(currLambda * currNSE) * UpperIncompleteGamma(1-currAlpha, currLambda * currNSE);
}

//Log probability of elongation at current codon
double PANSEModel::elongationProbabilityLog(double currAlpha, double currLambda, double currNSE)
{ 
    double x = (currLambda * currNSE) + currAlpha * (std::log(currLambda) + std::log(currNSE)) + UpperIncompleteGammaLog(1- currAlpha, currLambda * currNSE);
    return x;
}






double PANSEModel::elongationUntilIndexApproximation1Probability(double alpha, double lambda, double v, double current)
{
	return (current+(alpha/(lambda * v)));
}

double PANSEModel::elongationUntilIndexApproximation2Probability(double alpha, double lambda, double v, bool proposed)
{
    if (proposed)
    {
        propSigmaCalculationSummationFor2 += (alpha/(lambda * v)) *
                (-1*elongationUntilIndexApproximation1Probability(alpha, lambda, v, proposed))
                + (alpha/(lambda * lambda * v * v));
        return 1 + propSigmaCalculationSummationFor2;
    }
    else
    {
        currSigmaCalculationSummationFor2 += (alpha/(lambda * v)) *
                (-1*elongationUntilIndexApproximation1Probability(alpha, lambda, v, proposed))
                + (alpha/(lambda * lambda * v * v));
        return 1 + currSigmaCalculationSummationFor2;
    }
}

double PANSEModel::elongationUntilIndexApproximation1ProbabilityLog(double alpha, double lambda, double v, double current)
{
    return current -= (alpha/(lambda * v));   	   
}
double PANSEModel::elongationUntilIndexApproximation2ProbabilityLog(double alpha, double lambda, double v, double current)
{
	return current + (-(alpha/(lambda * v)) + (alpha/(lambda * lambda * v * v))
                      + ((alpha/(lambda * v)) * (alpha/(lambda * v))) / 2);
}


//TODO: Add sigma to parameter object to keep it from being calulcated and initialize as need using MCMC adjust functions to reflect this
double psi2phi(double psi, double sigma){
    return sigma * psi;
}
double phi2psi(double phi, double sigma){
    return phi / sigma;
}
