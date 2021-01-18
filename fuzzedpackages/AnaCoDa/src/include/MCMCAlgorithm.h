#ifndef MCMCALGORITHM_H
#define MCMCALGORITHM_H


#include "ROC/ROCModel.h"
#include "PA/PAModel.h"
#include "PANSE/PANSEModel.h"
#include "FONSE/FONSEModel.h"
#include "SequenceSummary.h"


#include <vector>
#include <cstdlib>
#include <sstream>
#include <chrono>
#include <fstream>
#include <stdlib.h> //can be removed later

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class MCMCAlgorithm
{
	private:
		unsigned samples;
		unsigned thinning;
		unsigned adaptiveWidth;
		unsigned lastConvergenceTest;
		int stepsToAdapt;


		bool estimateSynthesisRate;
		bool estimateCodonSpecificParameter;
		bool estimateHyperParameter;
		bool estimateMixtureAssignment;
		bool writeRestartFile;


		std::vector<double> posteriorTrace;
		std::vector<double> likelihoodTrace;
		std::vector<double> tmp;


		std::string file;
		unsigned fileWriteInterval;
		bool multipleFiles;


		//Acceptance Rejection Functions:
		double acceptRejectSynthesisRateLevelForAllGenes(Genome& genome, Model& model, int iteration);
		void acceptRejectCodonSpecificParameter(Genome& genome, Model& model, int iteration);
		void acceptRejectCodonSpecificParameter(Genome& genome, PANSEModel& model, int iteration);
		void acceptRejectHyperParameter(Genome &genome, Model& model, unsigned iteration);

	public:

		//Constructors & Destructors:
		explicit MCMCAlgorithm();
		MCMCAlgorithm(unsigned samples, unsigned thinning, unsigned _adaptiveWidth = 100,
					  bool _estimateSynthesisRate = true, bool _estimateCodonSpecificParameter = true,
					  bool _estimateHyperParameter = true);
		virtual ~MCMCAlgorithm();
	

		//MCMC Functions:
		void run(Genome& genome, Model& model, unsigned numCores = 1u, unsigned divergenceIterations = 0u); //TODO: UNTESTED
		void run_PANSE(Genome& genome, PANSEModel& model, unsigned numCores = 1u, unsigned divergenceIterations = 0u);
		void varyInitialConditions(Genome& genome, Model& model, unsigned divergenceIterations); //TODO: UNTESTED
		double calculateGewekeScore(unsigned current_iteration); //TODO: UNTESTED

		bool isEstimateSynthesisRate();
		bool isEstimateCodonSpecificParameter();
		bool isEstimateHyperParameter();
		bool isEstimateMixtureAssignment();

		void setEstimateSynthesisRate(bool in);
		void setEstimateCodonSpecificParameter(bool in);
		void setEstimateHyperParameter(bool in);
		void setEstimateMixtureAssignment(bool in);

		void setRestartFileSettings(std::string filename, unsigned interval, bool multiple); //TODO: UNTESTED
		void setStepsToAdapt(unsigned steps);
		int getStepsToAdapt();

		std::vector<double> getLogPosteriorTrace();
		std::vector<double> getLogLikelihoodTrace();
		double getLogPosteriorMean(unsigned samples); //TODO: UNTESTED

		static std::vector<double> acf(std::vector<double>& x, int nrows, int ncols, int lagmax, bool correlation, bool demean); //Currently unused. TODO: UNTESTED
		static std::vector<std::vector<double>> solveToeplitzMatrix(int lr, std::vector<double> r, std::vector<double> g); //Currently unused. TODO: UNTESTED





		//R Section:

#ifndef STANDALONE

		//Other Functions (All tested):
    	unsigned getSamples();
    	unsigned getThinning();
    	unsigned getAdaptiveWidth();
    	void setSamples(unsigned _samples);
    	void setThinning(unsigned _thinning);
    	void setAdaptiveWidth(unsigned _adaptiveWidth);
		void setLogPosteriorTrace(std::vector<double> _posteriorTrace);
		void setLogLikelihoodTrace(std::vector<double> _likelihoodTrace);
#endif //STANDALONE


	protected:
};

#endif // MCMCALGORITHM_H


/* ----------------------------- RCPP NOTE !!!!! -------------------------------------------- */
//                                                                                            //
// The functions declared and defined in the standalone block could be moved up to the C++    //
// side. These are R specific because of how RCPP deals with constructors, but having these   //
// functions in C++ would not affect the framework.                                           //
/* ----------------------------- RCPP NOTE !!!!! -------------------------------------------- */
