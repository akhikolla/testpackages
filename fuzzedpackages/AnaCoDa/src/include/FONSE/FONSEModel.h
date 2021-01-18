#ifndef FONSEMODEL_H
#define FONSEMODEL_H

#include "../base/Model.h"
#include "FONSEParameter.h"

class FONSEModel : public Model
{
	private:
		FONSEParameter *parameter;

		double calculateLogLikelihoodRatioPerAA(Gene& gene, std::string grouping, double *mutation, double *selection, double phiValue,double a1_value);
		double calculateMutationPrior(std::string grouping, bool proposed = false);



	public:
		//Constructors & Destructors:
		FONSEModel(bool _withPhi = false, bool _fix_sEpsilon = false);
		virtual ~FONSEModel();

		std::string type = "FONSE";
		

		//Likelihood Ratio Functions:
		virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio);
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures); 
		virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> &logProbabilityRatio);



		//Initialization and Restart Functions:
		virtual void initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate = true);
		virtual void writeRestartFile(std::string filename);



		//Category Functions:
		virtual double getCategoryProbability(unsigned i);
		virtual unsigned getMutationCategory(unsigned mixture);
		virtual unsigned getSelectionCategory(unsigned mixture);
		virtual unsigned getSynthesisRateCategory(unsigned mixture);
		virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k);



		//Group List Functions:
		virtual unsigned getGroupListSize(); //TODO: make not hardcoded?
		virtual std::string getGrouping(unsigned index);



		//stdDevSynthesisRate Functions:
		virtual double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false);
		virtual double getCurrentStdDevSynthesisRateProposalWidth();
		virtual void updateStdDevSynthesisRate();



		//Synthesis Rate Functions:
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false);
		virtual void updateSynthesisRate(unsigned i, unsigned k);



		//Iteration Functions:
		virtual unsigned getLastIteration();
		virtual void setLastIteration(unsigned iteration);



		//Trace Functions:
		virtual void updateStdDevSynthesisRateTrace(unsigned sample);
		virtual void updateInitiationCostParameterTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i);
		virtual void updateMixtureProbabilitiesTrace(unsigned sample);
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);
		virtual void updateHyperParameterTraces(unsigned sample);
		virtual void updateTracesWithInitialValues(Genome &genome);



		//Adaptive Width Functions:
		virtual void adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptInitiationCostProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt = true);
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt = true);



		//Other Functions:
		virtual void proposeCodonSpecificParameter();
		virtual void proposeHyperParameters();
		virtual void proposeSynthesisRateLevels();

		virtual double getInitiationCost(bool proposed);
		virtual double getCurrentInitiationCostProposalWidth();
		virtual void updateInitiationCost();

		virtual unsigned getNumPhiGroupings();
		virtual unsigned getMixtureAssignment(unsigned index);
		virtual unsigned getNumMixtureElements();
		virtual unsigned getNumSynthesisRateCategories();

		virtual void setNumPhiGroupings(unsigned value);
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene);
		virtual void setCategoryProbability(unsigned mixture, double value);

		virtual void updateCodonSpecificParameter(std::string grouping);
		virtual void completeUpdateCodonSpecificParameter();

		//virtual void updateGibbsSampledHyperParameters(Genome &genome);
		virtual void updateAllHyperParameter();
		virtual void updateHyperParameter(unsigned hp);

		virtual void simulateGenome(Genome &genome);
		virtual void printHyperParameters();
		FONSEParameter* getParameter();
		void setParameter(FONSEParameter &_parameter);
		virtual double calculateAllPriors();
		void calculateLogCodonProbabilityVector(unsigned numCodons, unsigned position, unsigned minIndexValue,
					double* mutation, double* selection, double phi, double a1_value, std::vector <double> &codonProb);
		void calculateCodonProbabilityVector(unsigned numCodons, unsigned position, double* mutation, double* selection,
					double phi, double a1_value, double codonProb[]);
		virtual void getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal,
					double* returnValue);

		virtual double getNoiseOffset(unsigned index, bool proposed = false);
		virtual double getObservedSynthesisNoise(unsigned index) ;
		virtual double getCurrentNoiseOffsetProposalWidth(unsigned index);
		virtual void updateNoiseOffset(unsigned index);
		virtual void updateNoiseOffsetTrace(unsigned sample);
		virtual void updateObservedSynthesisNoiseTrace(unsigned sample);
		virtual void adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void updateGibbsSampledHyperParameters(Genome &genome);

		//R Section:
#ifndef STANDALONE
		std::vector<double> CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi, double a1_value, unsigned position);
#endif

	protected:

		
};

#endif // FONSEMODEL_H
