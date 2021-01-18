#ifndef MODEL_H
#define MODEL_H


#include "../Genome.h"
#include "Parameter.h"

class Model
{
    private:

		double calculatePriorForCodonSpecificParam(Parameter *parameter, std::string grouping, unsigned paramType,
					bool proposed = false);


    public:
		//Constructors & Destructors:
        explicit Model();
		// TODO: Rule of Three dictates we may need a copy assignment operator as well (operator=)
        virtual ~Model();
        Model& operator=(const Model& rhs);

        std::string type;

        //Likelihood Ratio Functions:
        virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k,
					double* logProbabilityRatio) = 0;
        virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
         			std::vector<double> &logAcceptanceRatioForAllMixtures) = 0;
		virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration,
					std::vector <double> &logProbabilityRatio) = 0;

		virtual double calculateAllPriors() = 0;



		//Initialization and Restart Functions:
		virtual void initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate = true) = 0;
		virtual void writeRestartFile(std::string filename) = 0;



		//Category Functions:
		virtual double getCategoryProbability(unsigned i) = 0;
		virtual unsigned getMutationCategory(unsigned mixture) = 0;
		virtual unsigned getSelectionCategory(unsigned mixture) = 0;
		virtual unsigned getSynthesisRateCategory(unsigned mixture) = 0;
		virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k) = 0;



		//Group List Functions:
		virtual unsigned getGroupListSize() = 0;
		virtual std::string getGrouping(unsigned index) = 0;



		//stdDevSynthesisRate Functions:
		virtual double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false) = 0;
		//propose function?
		//set function?
		virtual double getCurrentStdDevSynthesisRateProposalWidth() = 0;
		virtual void updateStdDevSynthesisRate() = 0;



		//Synthesis Rate Functions:
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false) = 0;
		virtual void updateSynthesisRate(unsigned i, unsigned k) = 0;



		//Iteration Functions:
		virtual unsigned getLastIteration() = 0;
		virtual void setLastIteration(unsigned iteration) = 0;


		//Trace Functions:
		virtual void updateStdDevSynthesisRateTrace(unsigned sample) = 0;
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i) = 0;
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i) = 0;
		virtual void updateMixtureProbabilitiesTrace(unsigned sample) = 0;
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping) = 0;
		virtual void updateHyperParameterTraces(unsigned sample) = 0;
		virtual void updateTracesWithInitialValues(Genome &genome) = 0;

		//Adaptive Width Functions:
		virtual void adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt) = 0;
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt) = 0;
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration,
					bool adapt) = 0;
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt) = 0;

		//noise functions:
		virtual double getNoiseOffset(unsigned index, bool proposed = false) = 0;
		virtual double getObservedSynthesisNoise(unsigned index) = 0;
		virtual double getCurrentNoiseOffsetProposalWidth(unsigned index) = 0;
		virtual void updateNoiseOffset(unsigned index) = 0;
		virtual void updateNoiseOffsetTrace(unsigned sample) = 0;
		virtual void updateObservedSynthesisNoiseTrace(unsigned sample) = 0;
		virtual void adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt = true) = 0;
		virtual void updateGibbsSampledHyperParameters(Genome &genome) = 0;
		//Other Functions:
		virtual void proposeCodonSpecificParameter() = 0;
		virtual void proposeHyperParameters() = 0;
		virtual void proposeSynthesisRateLevels() = 0;

		virtual unsigned getNumPhiGroupings() = 0;
		virtual unsigned getMixtureAssignment(unsigned index) = 0;
		virtual unsigned getNumMixtureElements() = 0;
		virtual unsigned getNumSynthesisRateCategories() = 0;

		virtual void setNumPhiGroupings(unsigned value) = 0;
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene) = 0;
		virtual void setCategoryProbability(unsigned mixture, double value) = 0;

		virtual void updateCodonSpecificParameter(std::string grouping) = 0;
		virtual void completeUpdateCodonSpecificParameter() = 0;
		//virtual void updateGibbsSampledHyperParameters(Genome &genome) = 0;
		virtual void updateAllHyperParameter() = 0;
		virtual void updateHyperParameter(unsigned hp) = 0;

		virtual void simulateGenome(Genome &genome) =0;
		virtual void printHyperParameters() = 0;

	protected:
		bool withPhi;
		bool fix_sEpsilon;
		
};

#endif // MODEL_H
