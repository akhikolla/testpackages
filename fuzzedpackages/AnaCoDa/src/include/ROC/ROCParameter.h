#ifndef ROCPARAMETER_H
#define ROCPARAMETER_H


#include "../base/Trace.h"
#include "../base/Parameter.h"


#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <array>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class ROCParameter : public Parameter
{
	private:

		// std::vector <double> observedSynthesisNoise;

		// std::vector <double> noiseOffset_proposed;
		// std::vector <double> noiseOffset; //A_Phi
		// std::vector <double> std_NoiseOffset;
		// std::vector <double> numAcceptForNoiseOffset;
		std::vector<std::vector<double>> mutationVariance;
		std::vector<std::vector<double>> selectionVariance;
		double bias_csp;

		//double mutation_prior_sd;
		std::vector<std::vector<double>> mutation_prior_mean;
		std::vector<std::vector<double>> mutation_prior_sd;

		bool fix_dM=false;
		bool fix_dEta=false;
    	bool propose_by_prior=false;
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B);
	public:





		//Constructors & Destructors:
		explicit ROCParameter();
		explicit ROCParameter(std::string filename);
		ROCParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures,
					std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix,
					bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		ROCParameter& operator=(const ROCParameter& rhs);
		virtual ~ROCParameter();



		//Initialization, Restart, Index Checking:
		void initROCParameterSet();
		void initROCValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writeROCRestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes,bool estimateSynthesisRate = true);
		void initMutationCategories(std::vector<std::string> files, unsigned numCategories,bool fix = false);
		void initSelectionCategories(std::vector<std::string> files, unsigned numCategories,bool fix = false);



		//Trace Functions:
		// void updateObservedSynthesisNoiseTraces(unsigned sample);
		// void updateNoiseOffsetTraces(unsigned sample);
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);


		//Covariance Functions:
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);


		//observedSynthesisNoise Functions:
		// double getObservedSynthesisNoise(unsigned index);
		// void setObservedSynthesisNoise(unsigned index, double se);



		// //noiseOffset Functions:
		// double getNoiseOffset(unsigned index, bool proposed = false);
		// double getCurrentNoiseOffsetProposalWidth(unsigned index);
		// void proposeNoiseOffset();
		// void setNoiseOffset(unsigned index, double _NoiseOffset);
		// void updateNoiseOffset(unsigned index);

		// // noise Functions:
		// void setInitialValuesForSepsilon(std::vector<double> seps);


		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		void proposeCodonSpecificParameter();
        void setProposeByPrior(bool _propose_by_prior);
		void updateCodonSpecificParameter(std::string grouping);
		void completeUpdateCodonSpecificParameter();


		//Prior Functions:
		std::vector<std::vector<double>> getMutationPriorMean();
		std::vector<std::vector<double>> getMutationPriorStandardDeviation();
		std::vector<double> getMutationPriorMeanForCategory(unsigned category);
		std::vector<double> getMutationPriorStandardDeviationForCategory(unsigned category);
		void getMutationPriorMeanForCategoryForGroup(unsigned category, std::string aa, double *returnSet);
		void getMutationPriorStandardDeviationForCategoryForGroup(unsigned category, std::string aa, double *returnSet);
		void setMutationPriorMean(std::vector<std::vector<double>> _mutation_prior_mean);
		void setMutationPriorStandardDeviation(std::vector<std::vector<double>> _mutation_prior_sd);



		//Posterior, Variance, and Estimates Functions:
		//double getNoiseOffsetPosteriorMean(unsigned index, unsigned samples);
		//double getNoiseOffsetVariance(unsigned index, unsigned samples, bool unbiased = true);

		//Adaptive Width Functions:
		//void adaptNoiseOffsetProposalWidth(unsigned adaptationWidth, bool adapt);

		//Other Functions:
		//void setNumObservedPhiSets(unsigned _phiGroupings);
		void getParameterForCategory(unsigned category, unsigned parameter, std::string aa, bool proposal, double *returnValue);

		void fixDM();
		void fixDEta();
		bool isDMFixed();
		bool isDEtaFixed();


		//R Section:

#ifndef STANDALONE

		//Constructors & Destructors:
		ROCParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix,
					bool splitSer = true);
		ROCParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
					bool splitSer = true, std::string _mutationSelectionState = "allUnique");



		//Initialization, Restart, Index Checking:
		void initCovarianceMatrix(SEXP matrix, std::string aa);
		void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa);
		void initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa);
		void setMutationPriorMeanR(std::vector<double> _mutation_prior_mean);
		void setMutationPriorStandardDeviationR(std::vector<double> _mutation_prior_sd);

		//CSP Functions:
		std::vector<std::vector<double>> getProposedMutationParameter();
		std::vector<std::vector<double>> getCurrentMutationParameter();
		std::vector<std::vector<double>> getProposedSelectionParameter();
		std::vector<std::vector<double>> getCurrentSelectionParameter();


		void setProposedMutationParameter(std::vector<std::vector<double>> _proposedMutationParameter);
		void setCurrentMutationParameter(std::vector<std::vector<double>> _currentMutationParameter);
		void setProposedSelectionParameter(std::vector<std::vector<double>> _proposedSelectionParameter);
		void setCurrentSelectionParameter(std::vector<std::vector<double>> _currentSelectionParameter);

#endif //STANDALONE

	protected:
};

#endif // ROCPARAMETER_H
