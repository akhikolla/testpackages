#ifndef FONSEPARAMETER_H
#define FONSEPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "../base/Parameter.h"
#include "../base/Trace.h"

class FONSEParameter : public Parameter
{
	private:

		double bias_csp;
		double mutation_prior_sd;
		double a1;
		double a1_proposed;
		double std_a1;
        unsigned numAcceptForA1;

        bool fix_dM=false;
		bool fix_dOmega=false;
        bool fix_a1 = false;

		std::vector <double> propose(std::vector <double> currentParam, double(*proposal)(double a, double b), double A, std::vector <double> B);


	public:

		//Constructors & Destructors:
		FONSEParameter();
		explicit FONSEParameter(std::string filename);
		FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector <unsigned> geneAssignment,
				   std::vector <std::vector <unsigned> > thetaKmatrix, bool splitSer = true,
				   std::string _mutationSelectionState = "allUnique",double _a1 = 4);
		FONSEParameter& operator=(const FONSEParameter& rhs);
		FONSEParameter(const FONSEParameter &other); //TODO: No longer needed?
		virtual ~FONSEParameter();



		//Initialization, Restart, Index Checking:
		void initFONSEParameterSet(double _a1 = 4);
		void initFONSEValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writeFONSERestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate = true);
		void initMutationCategories(std::vector<std::string> files, unsigned numCategories, bool fix = false);
		void initSelectionCategories(std::vector<std::string> files, unsigned numCategories, bool fix = false);



		//Trace Functions:
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);
		void updateInitiationCostParameterTrace(unsigned sample);

		//Covariance Functions:
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);

		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping);
		void completeUpdateCodonSpecificParameter();


		//Prior Functions:
		double getMutationPriorStandardDeviation();
		void setMutationPriorStandardDeviation(double _mutation_prior_sd);

		double getInitiationCost(bool proposed);

		//Other functions:
		void getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal, double *returnSet);
		void proposeHyperParameters();
		void updateInitiationCost();
		double getCurrentInitiationCostProposalWidth();
		void adaptInitiationCostProposalWidth(unsigned adaptationWidth, bool adapt);

		void fixDM();
		void fixDOmega();
		bool isDMFixed();
		bool isDOmegaFixed();
		void fixedInitiationCost();


		//R Section:

#ifndef STANDALONE
		//Constructors & Destructors:
		FONSEParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix,
		 				bool splitSer = true, double _a1 = 4);
		FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
						bool splitSer = true, std::string _mutationSelectionState = "allUnique", double _a1 = 4);



		//Initialization, Restart, Index Checking:
		void initCovarianceMatrix(SEXP matrix, std::string aa);
		void initMutation(std::vector <double> mutationValues, unsigned mixtureElement, std::string aa);
		void initSelection(std::vector <double> selectionValues, unsigned mixtureElement, std::string aa);
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories,
						std::string paramType);


		//CSP Functions:
		std::vector< std::vector <double> > getCurrentMutationParameter();
		std::vector< std::vector <double> > getCurrentSelectionParameter();
		void setCurrentMutationParameter(std::vector<std::vector<double>> _currentMutationParameter);
		void setCurrentSelectionParameter(std::vector<std::vector<double>> _currentSelectionParameter);



#endif //STANDALONE

	protected:
};
#endif // FONSEPARAMETER_H
