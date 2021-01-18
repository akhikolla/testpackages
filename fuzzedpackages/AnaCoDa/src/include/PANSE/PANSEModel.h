#ifndef PANSEMODEL_H
#define PANSEMODEL_H


#include "../base/Model.h"
#include "PANSEParameter.h"
#include "../SequenceSummary.h"
#include <sstream>

class PANSEModel: public Model
{
	private:
		PANSEParameter *parameter;
		unsigned RFPCountColumn;
		double currSigmaCalculationSummationFor1, currSigmaCalculationSummationFor2;
		double propSigmaCalculationSummationFor1, propSigmaCalculationSummationFor2;
		double calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
				unsigned currRFPObserved, double phiValue, double prevSigma, double lgamma_currAlpha, double log_currLambdaPrime, double log_phi,double lgamma_rfp_alpha);
		std::vector<std::vector<double>> lgamma_currentAlpha;
   		std::vector<std::vector<double>> log_currentLambda;
        std::vector<std::vector<std::vector<double>>> lgamma_rfp_alpha;
      
        std::vector<double> prob_successful;

        virtual void calculateZ(std::string grouping,Genome& genome,std::vector<double> &Z,std::string param);

        virtual void fillMatrices(Genome& genome);
        virtual void clearMatrices();

	public:
		//Constructors & Destructors:
		explicit PANSEModel(unsigned RFPCountColumn = 0u, bool _withPhi = false, bool _fix_sEpsilon = false);
		virtual ~PANSEModel();

		std::string type = "PANSE";

		//Likelihood Ratio Functions:
		virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k,
				double* logProbabilityRatio); // Depends on RFPCountColumn
		//purely a placeholder
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
				std::vector<double> &logAcceptanceRatioForAllMixtures);
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
				std::vector<double> &logAcceptanceRatioForAllMixtures,std::string param="Elongation"); // Depends on RFPCountColumn
		

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


        //Partition Function Functions:
        virtual double getPartitionFunction(unsigned mixture, bool proposed = false);
        virtual double getCurrentPartitionFunctionProposalWidth();
        virtual void updatePartitionFunction();


		//Synthesis Rate Functions:
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false);
		virtual void updateSynthesisRate(unsigned i, unsigned k);


		//Iteration Functions:
		virtual unsigned getLastIteration();
		virtual void setLastIteration(unsigned iteration);


		//Trace Functions:
		virtual void updateStdDevSynthesisRateTrace(unsigned sample);
        virtual void updatePartitionFunctionTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i);
		virtual void updateMixtureProbabilitiesTrace(unsigned sample);
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);
		virtual void updateHyperParameterTraces(unsigned sample);
		virtual void updateTracesWithInitialValues(Genome &genome);


		//Adaptive Width Functions:
		virtual void adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
        virtual void adaptPartitionFunctionProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt = true);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt = true);
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt = true);


		//Other Functions:
		virtual void proposeCodonSpecificParameter();
		virtual void proposeHyperParameters();
		virtual void proposeSynthesisRateLevels();
		virtual void completeUpdateCodonSpecificParameter();

		virtual unsigned getNumPhiGroupings();
		virtual unsigned getMixtureAssignment(unsigned index);
		virtual unsigned getNumMixtureElements();
		virtual unsigned getNumSynthesisRateCategories();

		virtual void setNumPhiGroupings(unsigned value);
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene);
		virtual void setCategoryProbability(unsigned mixture, double value);

		virtual void updateCodonSpecificParameter(std::string codon);
		virtual void updateCodonSpecificParameter(std::string codon,std::string param="Elongation");
		//virtual void updateGibbsSampledHyperParameters(Genome &genome);
		virtual void updateAllHyperParameter();
		virtual void updateHyperParameter(unsigned hp);

		virtual void simulateGenome(Genome &genome); // Depends on RFPCountColumn
		virtual void printHyperParameters();
		PANSEParameter* getParameter();

		void setParameter(PANSEParameter &_parameter);
		virtual double calculateAllPriors();
		virtual double getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal);

		double UpperIncompleteGammaHelper(double s, double x);
		double UpperIncompleteGamma(double s, double x);
		double UpperIncompleteGammaLog(double s, double x);
      
        double elongationProbability(double currAlpha, double currLambda, double currNSE);
        double elongationProbabilityLog(double currAlpha, double currLambda, double currNSE);
        
        double elongationUntilIndexApproximation1Probability(double alpha, double lambda, double v, double current);
        double elongationUntilIndexApproximation2Probability(double alpha, double lambda, double v, bool proposed);
        double elongationUntilIndexApproximation1ProbabilityLog(double alpha, double lambda, double v, double current);
        double elongationUntilIndexApproximation2ProbabilityLog(double alpha, double lambda, double v, double current);

        //Psi-Phi Conversion Functions
        double psi2phi(double psi, double sigma);
        double phi2psi(double phi, double sigma);


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

#endif // PANSEMODEL_H
