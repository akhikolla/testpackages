#include <R_ext/Print.h>
#include <R_ext/Error.h>

#include <sstream>

#include <Rcpp.h>

#include "problem.h"
#include "RmalschainsEvaluate.h"

#include "debug.h"
#include "ilocalsearch.h"
#include "cross.h"
#include "ea.h"
#include "ssga.h"
#include "cmaeshan.h"
#include "srandom.h"
#include "malschains.h"
#include "get_util.h"
#include <iostream>
#include <cstdio>
#include <cassert>

using namespace realea;
using namespace std;
using std::unique_ptr;

EvalBase *ev = NULL;

tFitness rFitnessFunc(const tGen *x, int n) {
	Rcpp::NumericVector par(n);

	for (int i = 0; i < n; i++) {

		par[i] = x[i];
		//	Rprintf("%Le\n", par[i]);
	}

	double res = ev->eval(par);

	//Rprintf("%Le\n", res);

	return res;
}

RcppExport SEXP RmalschainsWrapper(SEXP p_fcall, SEXP p_dim, SEXP p_lower, SEXP p_upper,
		SEXP p_rho, SEXP p_popsize, SEXP p_maxEval, SEXP p_argls,
		SEXP p_debugMA, SEXP p_istep, SEXP p_effort, SEXP p_alpha, SEXP p_targetValue, SEXP p_threshold, 
		SEXP p_optMin, SEXP p_initialpop, SEXP p_seed, SEXP p_lsOnly, SEXP p_lsParam1, SEXP p_lsParam2) {

	BEGIN_RCPP;

	unsigned int dim = Rcpp::as<unsigned int>(p_dim);
	unsigned int popsize = Rcpp::as<unsigned int>(p_popsize);
	unsigned int maxEval = Rcpp::as<unsigned int>(p_maxEval);
	unsigned int istep = Rcpp::as<unsigned int>(p_istep);

	Rcpp::NumericVector upper(p_upper);
	Rcpp::NumericVector lower(p_lower);

	unsigned int debugMA = Rcpp::as<unsigned int>(p_debugMA);

	if(debugMA == 0) {
		disable_print_info();
		disable_print_debug();

	} else if (debugMA == 1) {
		disable_print_info();
		enable_print_debug();

        } else {
		enable_print_info();
		enable_print_debug();
		set_InitVerbose();
	}

	double effort = Rcpp::as<double>(p_effort);
	double alpha = Rcpp::as<double>(p_alpha);

	double targetValue = Rcpp::as<double>(p_targetValue);
	double threshold = Rcpp::as<double>(p_threshold);
	bool optMin = Rcpp::as<bool>(p_optMin);

	//unsigned int maxEval = relMaxEval + istep;

	std::string arg_ls = Rcpp::as<std::string>(p_argls);

	unsigned long seed = Rcpp::as<unsigned long>(p_seed);

	bool lsOnly = Rcpp::as<bool>(p_lsOnly);

	Random random(new SRandom(seed));

	if (TYPEOF(p_fcall) == EXTPTRSXP) {
		ev = new EvalCompiled(p_fcall, p_rho);
	} else {
		ev = new EvalStandard(p_fcall, p_rho);
	}

	ProblemPtr prob(new Problem());
	prob->setDimension(dim);

	for(unsigned int i=0; i < dim; i++) {

		prob->setDomainValues(i, lower(i), upper(i), true);
	}

	prob->setOptimize(targetValue, threshold);
	prob->setMaxEval(maxEval);

	if(optMin)
		prob->setMinimize();
	else
		prob->setMaximize();

	//Rprintf("targetValue: %Le\n", targetValue);

	prob->setEval(rFitnessFunc);

	tChromosomeReal sol(dim);
	tFitness fitness = 0;

	DomainReal *domain = prob->getDomain();

	ILocalSearch *ls = get_LS(arg_ls, domain, &random);

        //

  if (dynamic_cast<CMAESHansen*>(ls) != NULL)
  {
        int lambda = Rcpp::as<int>(p_lsParam1);
        int mu = Rcpp::as<int>(p_lsParam2);
        ((CMAESHansen*) ls)->setPopsize(lambda);
        ((CMAESHansen*) ls)->setParentsSize(mu);
  }

	PopulationReal *pop;
	MALSChains *ma = NULL;
	SSGA *ssga = NULL;
	Hybrid *hybrid = NULL;
    EA *alg;

	if(lsOnly) {

		//the local search needs a population (at least cmaes and sw need it)
		pop = new PopulationReal(&random, popsize, popsize);
		pop->reset(domain);

	} else {

		ssga = new SSGA(&random);
		ssga->setCross(new CrossBLX(alpha));
		ssga->setMutation(new MutationBGA());
		ssga->setSelect(new SelectNAM(3));
		ssga->setReplacement(new ReplaceWorst());

		ma = new MALSChains(ssga, ls);

		if (debugMA) {
			ma->setDebug();
		}

		//There is a bug in cmaes, this line: fbestever = cmaes_Get(&params->evo, "fbestever"); 
		//sets fbestever to zero and the local search returns a fitness of zero.
		//ma->setRestart(new RestartBest());
		ma->setRestart(NULL);
		hybrid = ma;
		hybrid->setEffortRatio(effort);

		print_debug("RatioLS: %f\nIstep: %d\n", effort, istep);

		std::stringstream effortStr;
		effortStr << effort;

		std::stringstream maxEvalStr;
		maxEvalStr << maxEval;

		set_Effort(hybrid, effortStr.str());
		hybrid->setIntensity(istep);
		set_MaxEval(hybrid, maxEvalStr.str());

		if (popsize > 0) {
			print_debug("Popsize: %u\n", popsize);
			hybrid->setPopsize(popsize);
		}

		//EA alg(hybrid, prob);
		alg = new EA(hybrid, prob);

		//print_debug("sol%Le\n", sol[0]);

		ma->init();

		pop = ma->getPop();
	}

	if(p_initialpop != NULL && p_initialpop != R_NilValue) {

		Rcpp::NumericMatrix initialpop(p_initialpop);
		unsigned indsize = initialpop.ncol();

		if(indsize != dim) {

			Rf_warning("Problem with your initial population: not right number of dimensions. Not using it.\n");

		} else {

			//PopulationReal *pop = ma->getPop();

			for(int s=0; s < initialpop.nrow() && s < (int) popsize; s++) {

				tChromosomeReal ind(dim);

				for(unsigned int i=0;i<dim;i++) ind[i] = initialpop(s, i);

				tFitness fitness = prob->eval(ind);
				pop->change(s, ind, fitness);
			}

			if(lsOnly) {

				//if initialpop is given, use the best individual
				//as starting point of the local search

				unsigned pos = pop->getBest();
				tIndividualRealPtr best= pop->getInd(pos);

				tChromosomeReal bestsol= best->sol();
				copy(bestsol.begin(), bestsol.end(), sol.begin());
				fitness = best->perf();
			}

		}
	}

//        Running* running;
        unsigned nevalea = 0, nevalls = 0;
	int num_improvement_ea=0, num_improvement_ls=0;
	int num_total_ea=0, num_total_ls=0;
	tFitness improvement_alg=0, improvement_ls=0;
	double time_ms_alg=0, time_ms_ls=0, time_ms_ma=0;

	if(lsOnly) {

		//run the local search

		Running* running = new Running(prob->getFinishCriterion());
		running->setMaxEval(prob->getMaxEval());

		ls->setPopulation(pop);
		ls->setProblem(prob.get());
		ls->setRandom(&random);
		ls->setRunning(running);
		ls->setEval(prob.get());

		ILSParameters *params = ls->getInitOptions(sol);

		tFitness fitness_old = 0;

		//run the local search in steps of 100 and break if no more improvement is present
		while(nevalls < maxEval) {

			ls->apply(params, sol, fitness, 100);

			if(abs(fitness - fitness_old) < threshold) break;

			nevalls += 100;
			fitness_old = fitness;
		}
	} else {

		ma->realApply(sol, fitness);

                //running = ma->getRunning();

                nevalea = ma->getNumEvalEA();
                nevalls = ma->getNumEvalLS();


                num_improvement_ea = ma->getNumImprovementEA();
                num_improvement_ls = ma->getNumImprovementLS();
                num_total_ea = ma->getNumTotalEA();
                num_total_ls = ma->getNumTotalLS();
                improvement_alg = ma->getImprovementEA();
                improvement_ls = ma->getImprovementLS();
                time_ms_alg = ma->getTimeMsEA();
                time_ms_ls = ma->getTimeMsLS();
                time_ms_ma = ma->getTimeMsMA();

	}

  if (dynamic_cast<CMAESHansen*>(ls) != NULL)
  {
        print_debug("CMAES::Popsize/Lambda: %d\n", ((CMAESHansen*) ls)->getPopsize());
        print_debug("CMAES::ParentsSize/Mu: %d\n", ((CMAESHansen*) ls)->getParentsSize());
  }

	//print_debug("%Le\n", fitness);

        //unsigned int neval = running->numEval();

	return Rcpp::List::create(Rcpp::Named("numEvalEA") = nevalea, 
                                  Rcpp::Named("numEvalLS") = nevalls,
                                  Rcpp::Named("numImprovementEA") = num_improvement_ea,
                                  Rcpp::Named("numImprovementLS") = num_improvement_ls,
                                  Rcpp::Named("numTotalEA") = num_total_ea,
                                  Rcpp::Named("numTotalLS") = num_total_ls,
                                  Rcpp::Named("improvementEA") = improvement_alg,
                                  Rcpp::Named("improvementLS") = improvement_ls,
                                  Rcpp::Named("timeMsEA") = time_ms_alg,
                                  Rcpp::Named("timeMsLS") = time_ms_ls,
                                  Rcpp::Named("timeMsMA") = time_ms_ma,
                                  Rcpp::Named("fitness") = fitness, 
                                  Rcpp::Named("sol") = sol);

	//To suppress the initialized but not used compiler warning
	(void)alg;

	END_RCPP;
}
