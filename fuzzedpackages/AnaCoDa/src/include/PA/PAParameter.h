#ifndef PAParameter_H
#define PAParameter_H


#include "../base/Trace.h"
#include "../base/Parameter.h"


#include <vector>
#include <random>
#include <string>
#include <iostream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class PAParameter: public Parameter {
	private:

		std::vector<std::vector<double>> lambdaValues; //Currently not used.
		double bias_csp;
        bool div_flag;

	public:
		static const unsigned dalpha;
		static const unsigned dlambdaprime;

		//Constructors & Destructors:
		explicit PAParameter();
		PAParameter(std::string filename);
		PAParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true,
				std::string _mutationSelectionState = "allUnique");
		PAParameter& operator=(const PAParameter& rhs);
		virtual ~PAParameter();


		//Initialization, Restart, Index Checking:
		void initPAParameterSet();
		void initRFPValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writePARestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate = true);
		void initAlpha(double alphaValue, unsigned mixtureElement, std::string codon); //R?
		void initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon); //R?
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories,
				unsigned paramType); //TODO: function needs to be renamed


		//Trace Functions:
		void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);


		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned index);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping);
		void completeUpdateCodonSpecificParameter();


		//Adaptive Width Functions:
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt); //may make virtual


		//Other functions:
		double getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal);





		//R Section:

#ifndef STANDALONE

		//Constructors & Destructors:
		PAParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment,
			std::vector<unsigned> _matrix, bool splitSer = true);
		PAParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
			bool splitSer = true, std::string _mutationSelectionState = "allUnique");



		//Initialization, Restart, Index Checking:
		void initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);

		//CSP Functions:
		std::vector<std::vector<double>> getProposedAlphaParameter();
		std::vector<std::vector<double>> getProposedLambdaPrimeParameter();
		std::vector<std::vector<double>> getCurrentAlphaParameter();
		std::vector<std::vector<double>> getCurrentLambdaPrimeParameter();
		void setProposedAlphaParameter(std::vector<std::vector<double>> alpha);
		void setProposedLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime);
		void setCurrentAlphaParameter(std::vector<std::vector<double>> alpha);
		void setCurrentLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime);



		//Posterior, Variance, and Estimates Functions:
		double getAlphaPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getLambdaPrimePosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);

		double getAlphaVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
		double getLambdaPrimeVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);



		//Other Functions:
		double getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal);

#endif //STANDALONE

	protected:

};

#endif // PAParameter_H
