#ifndef MANIFOLD_OPTIM_H
#define MANIFOLD_OPTIM_H

#include <RcppArmadillo.h>
#include "ForDebug.h"
#include <iostream>
#include "randgen.h"
#include "Manifold.h"
#include "ManifoldOptimProblem.h"
#include "ProblemAdapter.h"
#include "SolversLS.h"
#include <ctime>

#include "def.h"
#include "Util.h"
#include "ProductElement.h"
#include "ProductManifold.h"
#include "ManifoldFactory.h"
#include "VariableFactory.h"
#include "SolverFactory.h"

Rcpp::List ManifoldOptim(const arma::vec& initX, const arma::mat& initH,
	ManifoldOptimProblem& prob,
	const Rcpp::List& maniDefn,
	const Rcpp::List& maniParams,
	const Rcpp::List& solverParams,
	const Rcpp::List& derivParams,
	int method,
	bool hasHHR, bool verbose);

void ParseManiDefn(const Rcpp::List& maniDefn, Manifold*& oManifold, Element*& oElement);

#endif
