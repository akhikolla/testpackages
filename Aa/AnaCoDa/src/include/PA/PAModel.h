#ifndef PAMODEL_H
#define PAMODEL_H


#include "../base/Model.h"
#include "PAParameter.h"

#include <sstream>

class PAModel: public Model
{
	private:
		PAParameter *parameter;
		unsigned RFPCountColumn;



		double calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
				unsigned currRFPValue, unsigned currNumCodonsInMRNA, double phiValue);
		virtual void calculateZ(std::string grouping,Genome& genome,std::vector<double> &Z);

	public:
		//Constructors & Destructors:
		explicit PAModel(unsigned RFPCountColumn = 0u, bool _withPhi = false, bool _fix_sEpsilon = false);
		virtual ~PAModel();

		std::string type = "PA";



		//Likelihood Ratio Functions:
		virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k,
				double* logProbabilityRatio); // Depends on RFPCountColumn
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
				std::vector<double> &logAcceptanceRatioForAllMixtures); // Depends on RFPCountColumn
		virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration,
				std::vector <double> &logProbabilityRatio);


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
		virtual unsigned getGroupListSize();
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
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i);
		virtual void updateMixtureProbabilitiesTrace(unsigned sample);
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);
		virtual void updateHyperParameterTraces(unsigned sample);
		virtual void updateTracesWithInitialValues(Genome &genome);


		//Adaptive Width Functions:
		virtual void adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt = true);
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt = true);


		//Other Functions:
		virtual void proposeCodonSpecificParameter();
		virtual void proposeHyperParameters();
		virtual void proposeSynthesisRateLevels();

		virtual unsigned getNumPhiGroupings();
		virtual unsigned getMixtureAssignment(unsigned index);
		virtual unsigned getNumMixtureElements();
		virtual unsigned getNumSynthesisRateCategories();

		virtual void setNumPhiGroupings(unsigned value);
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene);
		virtual void setCategoryProbability(unsigned mixture, double value);

		virtual void updateCodonSpecificParameter(std::string aa);
		virtual void completeUpdateCodonSpecificParameter();
		//virtual void updateGibbsSampledHyperParameters(Genome &genome);
		virtual void updateAllHyperParameter();
		virtual void updateHyperParameter(unsigned hp);


		virtual void simulateGenome(Genome &genome); // Depends on RFPCountColumn
		virtual void printHyperParameters();
		//virtual void printCodonSpecificParameters();
		PAParameter* getParameter();
		void setParameter(PAParameter &_parameter);
		virtual double calculateAllPriors();
		virtual double getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal);

	    virtual double getNoiseOffset(unsigned index, bool proposed = false);
		virtual double getObservedSynthesisNoise(unsigned index) ;
		virtual double getCurrentNoiseOffsetProposalWidth(unsigned index);
		virtual void updateNoiseOffset(unsigned index);
		virtual void updateNoiseOffsetTrace(unsigned sample);
		virtual void updateObservedSynthesisNoiseTrace(unsigned sample);
		virtual void adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void updateGibbsSampledHyperParameters(Genome &genome);

	protected:
		

};

#endif // PAModel_H
