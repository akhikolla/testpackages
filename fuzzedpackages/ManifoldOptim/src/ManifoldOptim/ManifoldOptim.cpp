#include "ManifoldOptim.h"

Rcpp::List ManifoldOptim(const arma::vec& initX, const arma::mat& initH,
	ManifoldOptimProblem& prob, const Rcpp::List& maniDefn,
	const Rcpp::List& maniParams, const Rcpp::List& solverParams,
	const Rcpp::List& derivParams, const std::string& method, bool hasHHR)
{
	Manifold* domain;
	Variable* x;
	ParseManiDefn(maniDefn, domain, x);

	bool isCheckManiParams = Rcpp::as<bool>(maniParams["IsCheckParams"]);
	if (isCheckManiParams) {
		domain->CheckParams();
	}

	domain->SetHasHHR(hasHHR);
   
	// Initial iterate
	if (initX.size() == 0) {
		x->RandInManifold();
	} else {
		CopyFrom(x, initX);
	}

	// Initial H
	// For now, just leave it blank
	LinearOPE* linearOPE = nullptr;
	/*
	if (initH.n_rows > 0 && initH.n_cols > 0) {
		linearOPE = new LinearOPE(n);
		CopyFrom(linearOPE, initH);
	}
	*/

	// Set up the problem
	ProblemAdapter* probAdapter;
	ManifoldOptimProblem* probPtr = dynamic_cast<ManifoldOptimProblem*>(&prob);
	if (probPtr) {
		probAdapter = new ProblemAdapter(probPtr);
	} else {
		stop("Type of ManifoldOptimProblem could not be determined");
	}

	probAdapter->SetDomain(domain);

	// Set "epsilon" for numerical differentiation
	double epsNumericalGrad = Rcpp::as<double>(derivParams["EpsNumericalGrad"]);
	double epsNumericalHessEta = Rcpp::as<double>(derivParams["EpsNumericalHessEta"]);
	prob.SetEpsNumericalGrad(epsNumericalGrad);
	prob.SetEpsNumericalHessEta(epsNumericalHessEta);

	// Set up the solver
	PARAMSMAP solver_params;
	Rcpp::StringVector keys(solverParams.names());
	for (integer i = 0; i < solverParams.size(); i++) {
		const std::string& key = Rcpp::as<std::string>(keys(i));
		double val = Rcpp::as<double>(solverParams[key]);
		solver_params[key] = val;
	}

	// Get the solver here
	Solvers* solver = SolverFactory::GetSolver(method, probAdapter, x, linearOPE);
	solver->SetParams(solver_params);

	// TBD: We have not yet implemented the ability to provide your own line
	// search function from R
	/*
	RMEX::isstopped = mexProblem::GetFieldbyName(solverParams, 0, "IsStopped");
	if (RMEX::isstopped != nullptr) {
		solver->StopPtr = &RMEX::mexInnerStop;
	}

	RMEX::LinesearchInput = mexProblem::GetFieldbyName(solverParams, 0, "LinesearchInput");
	if (RMEX::LinesearchInput != nullptr) {
		SolversLS *solverLS = dynamic_cast<SolversLS *> (solver);
		if (solverLS != nullptr) {
			solverLS->LinesearchInput = &RMEX::mexLinesearchInput;
		}
	}
	*/

	bool isCheckSolverParams = Rcpp::as<bool>(solverParams["IsCheckParams"]);
	if (isCheckSolverParams) {
		solver->CheckParams();
	}

	solver->Run();

	const arma::vec& Xopt = ToArmaVec(solver->GetXopt());

	integer lenSeries = solver->GetlengthSeries();
	Rcpp::NumericVector funSeries(lenSeries);
	Rcpp::NumericVector gradSeries(lenSeries);
	Rcpp::NumericVector timeSeries(lenSeries);
	for (integer i = 0; i < lenSeries; i++) {
		funSeries[i] = solver->GetfunSeries()[i];
		gradSeries[i] = solver->GetgradSeries()[i];
		timeSeries[i] = solver->GettimeSeries()[i];
	}

	// ProblemAdapter keeps track of whether Hessian function was called by solver...
	// TBD: What is a nice way to tell if a numerical Hessian was used?
	Rcpp::String message;
	if (probAdapter->GetUsedHessian() && prob.GetUsedNumericalHessian()) {
	 	message = "Solver used a Hessian which was computed by numerical "
	 		"differentiation. Consider changing solvers or programming an "
	 		"analytical Hessian.";
	}

	Rcpp::List ret = Rcpp::List::create(Rcpp::Named("xopt", Xopt),
		Rcpp::Named("fval", solver->Getfinalfun()),
		Rcpp::Named("normgf", solver->Getnormgf()),
		Rcpp::Named("normgfgf0", solver->Getnormgfgf0()),
		Rcpp::Named("iter", solver->GetIter()),
		Rcpp::Named("num.obj.eval", solver->Getnf()),
		Rcpp::Named("num.grad.eval", solver->Getng()),
		Rcpp::Named("nR", solver->GetnR()),
		Rcpp::Named("nV", solver->GetnV()),
		Rcpp::Named("nVp", solver->GetnVp()),
		Rcpp::Named("nH", solver->GetnH()),
		Rcpp::Named("elapsed", solver->GetComTime()),
		Rcpp::Named("funSeries", funSeries),
		Rcpp::Named("gradSeries", gradSeries),
		Rcpp::Named("timeSeries", timeSeries),
		Rcpp::Named("message", message)
	);

	// Check Gradient and Hessian evaluated against both initial value and solution
	bool isCheckGradHess = Rcpp::as<bool>(solverParams["IsCheckGradHess"]);
	if (isCheckGradHess) {
		CopyFrom(x, initX);
		probAdapter->CheckGradHessian(x);

		CopyFrom(x, Xopt);
		probAdapter->CheckGradHessian(x);
	}

	delete solver;
	delete x;
	delete domain;
	delete probAdapter;

	// TBD: Do we need to do more clean-up for product manifold?
	/*
	for (integer i = 0; i < numoftype; i++) {
		delete manifolds[i];
		delete elements[i];
	}
	delete[] manifolds;
	delete[] elements;
	*/

	return ret;
}

void ParseManiDefn(const Rcpp::List& maniDefn, Manifold*& oManifold, Element*& oElement)
{
	integer numoftype = maniDefn.size();
	//integer powsinterval[numoftype + 1];
	integer *powsinterval = new integer[numoftype + 1];
	Manifold** manifolds = new Manifold*[numoftype];

	powsinterval[0] = 0;
	for (integer i = 0; i < numoftype; i++) {
		const Rcpp::List& mani = maniDefn[i];
		powsinterval[i+1] = powsinterval[i] + Rcpp::as<int>(mani["numofmani"]);
	}

	integer numoftotal = powsinterval[numoftype];
	Element** elements = new Element*[numoftotal];

	PARAMSMAP params;
	for (integer i = 0; i < numoftype; i++)
	{
		const Rcpp::List& mani = maniDefn[i];

		const std::string& name = Rcpp::as<std::string>(mani["name"]);
		integer n = Rcpp::as<int>(mani["n"]);
		integer m = Rcpp::as<int>(mani["m"]);
		integer p = Rcpp::as<int>(mani["p"]);
		params[static_cast<std::string> ("ParamSet")] = mani["ParamSet"];

		manifolds[i] = ManifoldFactory::GetManifold(name, n, m, p);
		manifolds[i]->SetParams(params);

		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++) {
			elements[j] = VariableFactory::GetVariable(name, n, m, p);
		}

		if (manifolds[i] == nullptr || elements[i] == nullptr) {
			throw ManifoldOptimException("Could not construct manifold, Please check your definition");
		}
	}

	if (numoftotal > 1) {
		// Return a Product Manifold and Product Manifold Element
		oManifold = new ProductManifold(manifolds, numoftype, powsinterval, numoftotal);
		oElement = new ProductElement(elements, numoftotal, powsinterval, numoftype);
	} else if (numoftotal == 1) {
		// Return the singleton manifold and element
		oManifold = manifolds[0];
		oElement = elements[0];
	} else {
		throw ManifoldOptimException("Manifold definition contained no entries");
	}

	delete[] manifolds;
	delete[] elements;
	delete[] powsinterval;
}

// TBD: This is related to line search functionality. We have not yet
// implemented it for R
/*
namespace RMEX
{
	mxArray *isstopped = nullptr;
	//This function defines the stopping criterion that may be used in the C++ solver
	bool mexInnerStop(Variable *x, Vector *gf, double f, double ngf, double ngf0)
	{
		mxArray *Xmx, *gfmx, *fmx, *ngfmx, *ngf0mx;
		mexProblem::ObtainMxArrayFromElement(Xmx, x);
		mexProblem::ObtainMxArrayFromElement(gfmx, gf);
		fmx = mxCreateDoubleScalar(f);
		ngfmx = mxCreateDoubleScalar(ngf);
		ngf0mx = mxCreateDoubleScalar(ngf0);

		mxArray *lhs[1], *rhs[6];
		rhs[0] = const_cast<mxArray *> (isstopped);
		rhs[1] = const_cast<mxArray *> (Xmx);
		rhs[2] = const_cast<mxArray *> (gfmx);
		rhs[3] = const_cast<mxArray *> (fmx);
		rhs[4] = const_cast<mxArray *> (ngfmx);
		rhs[5] = const_cast<mxArray *> (ngf0mx);
		mexCallMATLAB(1, lhs, 6, rhs, "feval");
		double result = mxGetScalar(lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(gfmx);
		mxDestroyArray(fmx);
		mxDestroyArray(ngfmx);
		mxDestroyArray(ngf0mx);
		mxDestroyArray(lhs[0]);
		return (result != 0);
	};

	mxArray *LinesearchInput = nullptr;

	// This function defines the line search algorithm that may be used in the C++ solver
	double mexLinesearchInput(Variable *x1, Vector *eta1, double initialstepsize, double initialslope)
	{
		mxArray *Xmx, *eta1mx, *tmx, *smx;
		mexProblem::ObtainMxArrayFromElement(Xmx, x1);
		mexProblem::ObtainMxArrayFromElement(eta1mx, eta1);
		tmx = mxCreateDoubleScalar(initialstepsize);
		smx = mxCreateDoubleScalar(initialslope);

		mxArray *lhs[1], *rhs[5];
		rhs[0] = const_cast<mxArray *> (LinesearchInput);
		rhs[1] = const_cast<mxArray *> (Xmx);
		rhs[2] = const_cast<mxArray *> (eta1mx);
		rhs[3] = const_cast<mxArray *> (tmx);
		rhs[4] = const_cast<mxArray *> (smx);
		mexCallMATLAB(1, lhs, 5, rhs, "feval");
		double result = mxGetScalar(lhs[0]);
		mxDestroyArray(Xmx);
	 	mxDestroyArray(eta1mx);
		mxDestroyArray(tmx);
		mxDestroyArray(smx);
		mxDestroyArray(lhs[0]);
		return result;
	}
};
*/
