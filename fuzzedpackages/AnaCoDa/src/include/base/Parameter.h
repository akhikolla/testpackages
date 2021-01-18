#ifndef PARAMETER_H
#define PARAMETER_H


#include "../Genome.h"
#include "../CovarianceMatrix.h"
#include "Trace.h"


#include <vector>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <ctime>
#include <sstream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

/* Note 1) -- on getSelectionCategory and getSynthesisRateCategory
 * These two functions are technically the same for readability.
 * Selection and Synthesis Rate are directly related even if they are not known
 * and thus are represented by the same variable. By splitting this
 * into Selection and Synthesis, avoids confusing the two, however.
*/

class Parameter {
	private:

		//STATICS - Sorting Functions:
		void quickSortPair(double a[], int b[], int first, int last);
		static int pivotPair(double a[], int b[], int first, int last);

		

		std::vector<double> codonSpecificPrior;

		bool fix_stdDevSynthesis = false;
	public:

		static const std::string allUnique;
		static const std::string selectionShared;
		static const std::string mutationShared;

		static const unsigned dM;
		static const unsigned dEta;
		static const unsigned dOmega;
		static const unsigned alp;
		static const unsigned lmPri;
		static const unsigned nse;

#ifdef STANDALONE
		static std::default_random_engine generator; // static to make sure that the same generator is used during the runtime.
#endif


		//Constructors & Destructors:
		Parameter();
		Parameter(unsigned maxGrouping);
		Parameter& operator=(const Parameter& rhs);
		virtual ~Parameter();

		std::vector <std::string> CSPToUpdate;

		//Initialization and Restart Functions: TODO: test
		void initParameterSet(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures,
			std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> mixtureDefinitionMatrix,
			bool splitSer = true, std::string _mutationSelectionState = "allUnique"); //Mostly tested; TODO caveats
		void initBaseValuesFromFile(std::string filename);
		void writeBasicRestartFile(std::string filename);
		void initCategoryDefinitions(std::string mutationSelectionState,
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix);
		void InitializeSynthesisRate(Genome& genome, double sd_phi);
		void InitializeSynthesisRate(double sd_phi);
		void InitializeSynthesisRate(std::vector<double> expression);
		std::vector<double> readPhiValues(std::string filename); //General function, possibly move


		//Prior functions: TODO: test
		double getCodonSpecificPriorStdDev(unsigned paramType);


		//Mixture Definition Matrix and Category Functions: Mostly tested, see comments
		void setNumMutationSelectionValues(std::string mutationSelectionState,
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix); //TODO: test
		void printMixtureDefinitionMatrix(); //Untested
		double getCategoryProbability(unsigned mixtureElement);
		void setCategoryProbability(unsigned mixtureElement, double value);
		unsigned getNumMutationCategories(); //TODO caveat
		unsigned getNumSelectionCategories(); //TODO caveat
		unsigned getNumSynthesisRateCategories(); //TODO caveat
		unsigned getMutationCategory(unsigned mixtureElement);
		unsigned getSelectionCategory(unsigned mixtureElement); //see Note 1) at top of file.
		unsigned getSynthesisRateCategory(unsigned mixtureElement); //see Note 1) at top of file.
		std::vector<unsigned> getMixtureElementsOfMutationCategory(unsigned category); //TODO caveat
		std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned category); //TODO caveat
		std::string getMutationSelectionState();
		unsigned getNumAcceptForCspForIndex(unsigned i); //Only for unit testing.


		//Group List Functions: All tested
		void setGroupList(std::vector<std::string> gl);
		std::string getGrouping(unsigned index);
		std::vector<std::string> getGroupList();
		unsigned getGroupListSize();


		//stdDevSynthesisRate Functions: Mostly tested, see comments.
		double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false);
		virtual void proposeStdDevSynthesisRate(); //TODO: test
		void fixStdDevSynthesis();
		void setStdDevSynthesisRate(double stdDevSynthesisRate, unsigned selectionCategory);
		double getCurrentStdDevSynthesisRateProposalWidth();
		unsigned getNumAcceptForStdDevSynthesisRate(); //Only for unit testing.
		void updateStdDevSynthesisRate(); //TODO: test
		double getStdCspForIndex(unsigned i); //Only for unit testing.


		//Synthesis Rate Functions: Mostly tested, see comments
		double getSynthesisRate(unsigned geneIndex, unsigned mixtureElement, bool proposed = false);
		double getCurrentSynthesisRateProposalWidth(unsigned expressionCategory, unsigned geneIndex);
		double getSynthesisRateProposalWidth(unsigned geneIndex, unsigned mixtureElement);
		void proposeSynthesisRateLevels(); //TODO: test
		void setSynthesisRate(double phi, unsigned geneIndex, unsigned mixtureElement);
		void updateSynthesisRate(unsigned geneIndex); //TODO: test
		void updateSynthesisRate(unsigned geneIndex, unsigned mixtureElement); //TODO: test
		unsigned getNumAcceptForSynthesisRate(unsigned expressionCategory, unsigned geneIndex); //Only for unit testing


		//Noise Functions...updating AnaCoDa to allow all models (instead of just ROC) to use empirical gene expression values to inform estimation of \phi

		double getObservedSynthesisNoise(unsigned index);
		void setObservedSynthesisNoise(unsigned index, double se);

		//noiseOffset Functions:
		double getNoiseOffset(unsigned index, bool proposed = false);
		double getCurrentNoiseOffsetProposalWidth(unsigned index);
		void proposeNoiseOffset();
		void setNoiseOffset(unsigned index, double _NoiseOffset);
		void updateNoiseOffset(unsigned index);

		void setNumObservedPhiSets(unsigned _phiGroupings);

		// noise Functions:
		void setInitialValuesForSepsilon(std::vector<double> seps);

		//Posterior, Variance, and Estimates Functions for noise:
		double getNoiseOffsetPosteriorMean(unsigned index, unsigned samples);
		double getNoiseOffsetVariance(unsigned index, unsigned samples, bool unbiased = true);

		//Adaptive Width Functions:
		void adaptNoiseOffsetProposalWidth(unsigned adaptationWidth, bool adapt);

		//Iteration Functions: All tested
		unsigned getLastIteration();
		void setLastIteration(unsigned iteration);


		//Trace Functions: TODO: test
		Trace& getTraceObject();
		void setTraceObject(Trace _trace);
		void updateObservedSynthesisNoiseTraces(unsigned sample);
		void updateNoiseOffsetTraces(unsigned sample);
		void updateStdDevSynthesisRateTrace(unsigned sample);
		void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		void updateMixtureProbabilitiesTrace(unsigned samples);


		//Adaptive Width Functions: TODO: test
		void adaptStdDevSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt);
		void adaptSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt);


		//Posterior, Variance, and Estimates Functions: TODO: test
		double getStdDevSynthesisRatePosteriorMean(unsigned samples, unsigned mixture);
		double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, bool log_scale=false);

		double getCodonSpecificPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, bool withoutReference = true, bool byGene = false, bool log_scale = false);
		double getStdDevSynthesisRateVariance(unsigned samples, unsigned mixture, bool unbiased);
		double getSynthesisRateVariance(unsigned samples, unsigned geneIndex,
			bool unbiased = true, bool log_scale = false);
		double getCodonSpecificVariance(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, bool unbiased, bool withoutReference = true, bool log_scale = false);
	        std::vector<double> getCodonSpecificQuantile(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, std::vector<double> probs, bool withoutReference, bool log_scale = false);
		std::vector<double> getExpressionQuantile(unsigned samples, unsigned geneIndex,
			std::vector<double> probs, bool log_scale = false);
		std::vector<double> calculateQuantile(std::vector<float> &parameterTrace, unsigned samples, std::vector<double> probs, bool log_scale=false);
		unsigned getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex);
		std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);


		//Other Functions: Mostly tested, see comments
		unsigned getNumParam();
		unsigned getNumMixtureElements();
		unsigned getNumObservedPhiSets();
		void setMixtureAssignment(unsigned gene, unsigned value);
		unsigned getMixtureAssignment(unsigned gene);
		virtual std::vector <std::vector <double> > calculateSelectionCoefficients(unsigned sample); //TODO: test


		//Static Functions: TODO: test
		static double calculateSCUO(Gene& gene);
		static void drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b),
			double* randomNumbers);
		static void drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumber);
		static double randNorm(double mean, double sd);
		static double randLogNorm(double m, double s);
		static double randExp(double r);
		static double randGamma(double shape, double rate);
		static void randDirichlet(std::vector <double> &input, unsigned numElements, std::vector <double> &output);
		static double randUnif(double minVal, double maxVal);
		static unsigned randMultinom(double *probabilities, unsigned mixtureElements);
		static unsigned randMultinom(std::vector <double> &probabilities, unsigned mixtureElements);
		static double densityNorm(double x, double mean, double sd, bool log = false);
		static double densityLogNorm(double x, double mean, double sd, bool log = false);
		//double getMixtureAssignmentPosteriorMean(unsigned samples, unsigned geneIndex);
		// TODO: implement variance function, fix Mean function (won't work with 3 groups)





		//R Section:

#ifndef STANDALONE

		//Initialization and Restart Functions:
		void initializeSynthesisRateByGenome(Genome& genome, double sd_phi);
		void initializeSynthesisRateByRandom(double sd_phi);
		void initializeSynthesisRateByList(std::vector<double> expression);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);



		//Mixture Definition Matrix and Category Functions:
		unsigned getMutationCategoryForMixture(unsigned mixtureElement);
		unsigned getSelectionCategoryForMixture(unsigned mixtureElement);
		unsigned getSynthesisRateCategoryForMixture(unsigned mixtureElement);
		std::vector<std::vector<unsigned>> getCategories();
		void setCategories(std::vector<std::vector<unsigned>> _categories);
		void setCategoriesForTrace();
		void setNumMutationCategories(unsigned _numMutationCategories);
		void setNumSelectionCategories(unsigned _numSelectionCategories);



		//Synthesis Rate Functions:
		std::vector<std::vector<double>> getSynthesisRateR();
		std::vector<double> getCurrentSynthesisRateForMixture(unsigned mixture);



		//Posterior, Variance, and Estimates Functions:
		double getSynthesisRatePosteriorMeanForGene(unsigned samples, unsigned geneIndex, bool log_scale);
		double getSynthesisRateVarianceForGene(unsigned samples, unsigned geneIndex, bool unbiased, bool log_scale);
		unsigned getEstimatedMixtureAssignmentForGene(unsigned samples, unsigned geneIndex);

		std::vector<double> getEstimatedMixtureAssignmentProbabilitiesForGene(unsigned samples, unsigned geneIndex);

		double getCodonSpecificPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon,
			unsigned paramType, bool withoutReference, bool log_scale = false);
		double getCodonSpecificVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon,
			unsigned paramType, bool unbiased, bool withoutReference, bool log_scale = false);
        	std::vector<double> getCodonSpecificQuantileForCodon(unsigned mixtureElement, unsigned samples,
        		std::string &codon, unsigned paramType, std::vector<double> probs, bool withoutReference, bool log_scale = false);
		std::vector<double> getExpressionQuantileForGene(unsigned samples,
			unsigned geneIndex, std::vector<double> probs, bool log_scale);



		//Other Functions:
		SEXP calculateSelectionCoefficientsR(unsigned sample);
		std::vector<unsigned> getMixtureAssignmentR();
		void setMixtureAssignmentR(std::vector<unsigned> _mixtureAssignment);
		unsigned getMixtureAssignmentForGeneR(unsigned geneIndex);
		void setMixtureAssignmentForGene(unsigned geneIndex, unsigned value);
		void setNumMixtureElements(unsigned _numMixtures);

#endif

	protected:
		Trace traces;

		unsigned adaptiveStepPrev;
		unsigned adaptiveStepCurr;
		
		std::vector<CovarianceMatrix> covarianceMatrix;
		std::vector<mixtureDefinition> categories;
		std::vector<double> categoryProbabilities;
		std::vector<std::vector<unsigned>> mutationIsInMixture;
		std::vector<std::vector<unsigned>> selectionIsInMixture;
		unsigned numMutationCategories; //TODO Probably needs to be renamed
		unsigned numSelectionCategories; //TODO Probably needs to be renamed
		std::vector<unsigned> numAcceptForCodonSpecificParameters;
		std::string mutationSelectionState; //TODO: Probably needs to be renamed

        //<Alpha or Lambda or Mutation or Selection < Mixture < Codon >>>
		std::vector<std::vector<std::vector<double>>> proposedCodonSpecificParameter;
		std::vector<std::vector<std::vector<double>>> currentCodonSpecificParameter;

		std::vector<unsigned> mixtureAssignment;
		std::vector<std::string> groupList;
		unsigned maxGrouping;


		std::vector<double> stdDevSynthesisRate_proposed;
		std::vector<double> stdDevSynthesisRate;
		double bias_stdDevSynthesisRate; //NOTE: Currently, this value is always set to 0.0
		double std_stdDevSynthesisRate;
		unsigned numAcceptForStdDevSynthesisRate;
		std::vector<double> std_csp;


		std::vector <double> observedSynthesisNoise;

		std::vector <double> noiseOffset_proposed;
		std::vector <double> noiseOffset; //A_Phi
		std::vector <double> std_NoiseOffset;
		std::vector <double> numAcceptForNoiseOffset;
		
        //Unknown indexing hoping (mixture) then gene
		std::vector<std::vector<double>> proposedSynthesisRateLevel;
		std::vector<std::vector<double>> currentSynthesisRateLevel;
		std::vector<std::vector<unsigned>> numAcceptForSynthesisRate;

		unsigned lastIteration;

		unsigned int numParam;
		unsigned numMixtures;
		unsigned obsPhiSets;

		double bias_phi; //NOTE: Currently, this value is always set to 0.0
		std::vector<std::vector<double>> std_phi;

};

#endif // PARAMETER_H
