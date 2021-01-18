#include "lps.h"
#include "debug.h"

namespace lps {
  // newton raphson to fit model
  std::pair<arma::colvec, unsigned> lps::newtonRaphson(bool input, 
						       arma::uvec index) const
  {
    if (!input)
      index = constants;
    // onlye one beta needed, so to reshape before iterating
    arma::colvec coef = arma::zeros<arma::mat>(p, 1);
    // use constants as index to be estimated
    unsigned iter;
    for (iter = 0; iter < param.maxIter; iter++) {
      arma::colvec grad;
      ptrLoss -> gradient(grad, coef, index);
      double gpnorm = norm(grad, 2);
      printIter(iter, gpnorm);
      if (gpnorm < param.gptol) {
	printIter(iter, gpnorm, 1);
	if (param.output)
	  Rcpp::Rcout << "*** Converged ***" << std::endl;
	break;
      }
      arma::mat hess;
      ptrLoss -> hessian(hess, coef, index);
      // std::cout << hess << std::endl;
      double scale_factor = sum(hess.diag()) / 
	static_cast<double> (index.n_rows);
      double damping = scale_factor < .1 * gpnorm ? scale_factor : .1 * gpnorm;
      hess += damping * arma::eye(index.n_rows, index.n_rows);
      arma::colvec substep = -solve(hess, grad);
      arma::colvec step = arma::zeros<arma::mat> (p, 1);
      step.elem(index) = substep;
      coef += step;
    }
    std::pair<arma::colvec, unsigned> ret;
    ret.first = coef;
    ret.second = iter;
    return ret;
  }

  // tuning
  double lps::tune(const arma::colvec& current,
		   int method = AIC)
  {
    arma::uvec nonzero = find(abs(current) > param.gptol);
    bool out = param.output;
    param.output = 0; 
    std::pair<arma::colvec, unsigned> ret;
    ret = newtonRaphson(1, nonzero);
    param.output = out;
    switch(method) {
    case AIC:
      return 2 * ptrLoss -> eval(ret.first) + 2 * nonzero.n_rows / static_cast<double> (n);
    case BIC:
      return 2 * ptrLoss -> eval(ret.first) + log(static_cast<double>(n)) * nonzero.n_rows / static_cast<double> (n);
    case GACV:
      {
	double sum = 0;
	arma::colvec fx1 = ptrLoss -> getlinear();
	// NOTES: W 
	arma::mat W = arma::zeros<arma::mat> (fx1.n_rows, fx1.n_rows);
	ptrLoss -> variance(W);
	// calculate randomized GACV for every replicate
	for (unsigned rep = 0; rep < param.numRep; rep++) {
	  arma::colvec tmpVec = epsilon.col(rep);
	  arma::colvec Y1 = ptrLoss -> getY();
	  ptrLoss -> addRand(tmpVec);
	  arma::colvec Y2 = ptrLoss -> getY();
	  std::pair<arma::colvec, unsigned> refit;
	  // turn output off
	  out = param.output;
	  param.output = 0;
	  refit = newtonRaphson(1, nonzero);
	  param.output = out;
	  arma::colvec diff = ptrLoss -> getlinear() - fx1;
	  tmpVec = -tmpVec;
	  ptrLoss -> addRand(tmpVec);
	  tmpVec = Y2 - Y1;
	  arma::colvec tmp = tmpVec.t() * diff / 
	    (tmpVec.t() * tmpVec -tmpVec.t() * W * diff);
	  sum += tmp(0);
	}

	arma::colvec meanVec;
	ptrLoss -> mean(meanVec);
	arma::colvec tmp = meanVec.t() * (meanVec - ptrLoss -> getY());
	sum *= tmp(0) / static_cast<double> (n) / static_cast<double> (param.numRep);
	return ptrLoss -> eval(ret.first) + sum; 
      }
    case BGACV:
      {
	double sum = 0;
	arma::colvec fx1 = ptrLoss -> getlinear();
	arma::mat W = arma::zeros<arma::mat> (fx1.n_rows, fx1.n_rows);
	ptrLoss -> variance(W);
	// calculate randomized GACV for every replicate
	for (unsigned rep = 0; rep < param.numRep; rep++) {
	  arma::colvec tmpVec = epsilon.col(rep);
	  arma::colvec Y1 = ptrLoss -> getY();
	  ptrLoss -> addRand(tmpVec);
	  arma::colvec Y2 = ptrLoss -> getY();
	  std::pair<arma::colvec, unsigned> refit;
	  // turn output off
	  out = param.output;
	  param.output = 0;
	  refit = newtonRaphson(1, nonzero);
	  param.output = out;
	  arma::colvec diff = ptrLoss -> getlinear() - fx1;
	  tmpVec = -tmpVec;
	  ptrLoss -> addRand(tmpVec);
	  tmpVec = Y2 - Y1;
	  arma::colvec tmp = tmpVec.t() * diff / 
	    (tmpVec.t() * tmpVec -tmpVec.t() * W * diff);
	  sum += tmp(0);
	}

	arma::colvec meanVec;
	ptrLoss -> mean(meanVec);
	arma::colvec tmp = meanVec.t() * (meanVec - ptrLoss -> getY());
	sum *= tmp(0) / static_cast<double> (n) / static_cast<double> (param.numRep);
	return ptrLoss -> eval(ret.first) + log(static_cast<double> (n)) / 2.0 * sum; 
      }
    }
    return 0; // -WALL
  }

  // step wise fitting
  void lps::stepfit()
  {
    beta.reshape(p, 1);
  }

  // fit for lasso pattern search
  void lps::gridSearch(const arma::mat& lambda)
  {
    scores.reshape(lambda.n_cols, 1);
    beta.reshape(p, lambda.n_cols);
    arma::colvec init = arma::zeros<arma::mat> (p, 1);
    init.elem(constants) = .01 * arma::ones<arma::mat> (constants.n_rows, 1);
    iterSpent.reshape(lambda.n_cols, 1);
    for (unsigned i = 0; i < lambda.n_cols; i++) {
      arma::colvec tmpBeta = init;
      arma::colvec currentLambda = lambda.col(i);
      iterSpent(i) = solveLPS(currentLambda, init, tmpBeta);
      beta.col(i) = tmpBeta;
      init = tmpBeta;
      scores(i) = tune(tmpBeta, param.tuneMethod);
    }
  }

  double lps::evalLambda (const arma::colvec& current, arma::colvec& tmpBeta)
  {
    arma::colvec init = tmpBeta;
    arma::colvec currentLambda (exp(current));
    int iters = solveLPS(currentLambda, init, tmpBeta);
    if (iters < 0 || currentLambda.max() < sqrt(param.gptol))  // threshold
      return 100;
    return tune(tmpBeta, param.tuneMethod);
  }

  int lps::bestLambda(const arma::mat& matLambda, const arma::colvec& target) const
  {
    int ret = -1;
    double minDist = 1e16;
    for (unsigned i = 0; i < matLambda.n_cols; i++) {
      double dist = 0;
      for (unsigned j = 0; j < matLambda.n_rows; j++)
	dist += pow((matLambda(j, i) - target(j)), 2);
      dist = sqrt(dist);
      if (dist < minDist) {
	ret = i;
	minDist = dist;
      }
    }
    return ret;
  }

  // Nelder Mead search over lambda
  void lps::nelderMead(const arma::colvec& lambda)
  {
    const double backStep = .2;
    unsigned dim = lambda.n_rows; // number of columns of Nelder Mead
    scores.reshape(dim + 1, 1);
    // log lambda scale
    arma::mat lambda_grid (dim, dim + 1);
    lambda_grid.col(0) = log(lambda);

    arma::colvec tmpBeta = arma::zeros<arma::mat> (p, 1);
    arma::colvec currentLambda = lambda_grid.col(0);
    // store beta values for warm start
    arma::mat storeBeta = arma::zeros<arma::mat> (p, dim + 1);
    int index = -1;
    scores(0) = evalLambda(currentLambda, tmpBeta);
    arma::uvec nonZero = arma::find(abs(tmpBeta) > 1e-6);

    while (1) {
      currentLambda -= backStep;
      scores(0) = evalLambda(currentLambda, tmpBeta);
      arma::uvec nonZero = arma::find(abs(tmpBeta) > 1e-6);
      if (nonZero.n_rows > constants.n_rows * 3)
	break;
    }
    arma::colvec startLambda = currentLambda + backStep;
    for (unsigned i = 1; i <= dim; i++) {
      arma::colvec unit = arma::zeros<arma::mat> (dim, 1);
      unit(i - 1) = backStep;
      lambda_grid.col(i) = startLambda - unit;
      currentLambda = lambda_grid.col(i);
      scores(i) = evalLambda(currentLambda, tmpBeta);
      storeBeta.col(i) = tmpBeta;
    }

    // parameters for NM
    double Alpha = 1, Beta = .5, Gamma = 2, Delta = .5; 
    unsigned iter;
    unsigned minPos = 0;
    for (iter = 0; iter < param.maxIter; iter++) {
      // start of the process for Nelder Mead
      // ordering
      arma::uvec sortOrder = sort_index(scores);
      //std::cout << lambda_grid << scores << std::endl;
      // check convergence
      arma::colvec lambda1 = lambda_grid.col(sortOrder(0));
      arma::colvec lambda2 = lambda_grid.col(sortOrder(1));
      if (std::abs(scores(sortOrder(0)) - scores(sortOrder(dim))) < sqrt(param.gptol) ||
	  norm(lambda1 - lambda2, 2) < sqrt(param.gptol)) {
	minPos = sortOrder(0);
	break;
      }
      arma::colvec centroid = arma::zeros<arma::mat> (dim, 1);
      // calculate centroid
      for (unsigned i = 0; i < dim; i++)
	centroid = centroid + lambda_grid.col(sortOrder(i));
      // compute reflect lambda and its value
      arma::colvec lambda_reflect = centroid + 
	Alpha * (centroid - lambda_grid.col(sortOrder(dim)));
      currentLambda = lambda_reflect;
      index = bestLambda(lambda_grid, currentLambda);
      tmpBeta = storeBeta.col(index);
      double score_reflect = evalLambda(currentLambda, tmpBeta);
      if (score_reflect >= scores(sortOrder(0)) && 
	  score_reflect <  scores(sortOrder(dim - 1))) {
	scores(sortOrder(dim)) = score_reflect;
	lambda_grid.col(sortOrder(dim)) = lambda_reflect;
	continue;
      } else if (score_reflect < scores(sortOrder(0))) {
	// expand
	arma::colvec lambda_expand = centroid + Gamma * (lambda_reflect - centroid);
	currentLambda = lambda_expand;
	index = bestLambda(lambda_grid, currentLambda);
	tmpBeta = storeBeta.col(index);
	double score_expand = evalLambda(currentLambda, tmpBeta);
	if (score_expand < score_reflect) {
	  scores(sortOrder(dim)) = score_expand;
	  lambda_grid.col(sortOrder(dim)) = lambda_expand;
	  continue;
	} else {
	  scores(sortOrder(dim)) = score_reflect;
	  lambda_grid.col(sortOrder(dim)) = lambda_reflect;
	  continue;
	}
      }
      if (score_reflect < scores(sortOrder(dim)) &&
	  score_reflect >= scores(sortOrder(dim - 1))) {
	// contract outside
	arma::colvec lambda_contract = centroid + Beta * (lambda_reflect - centroid);
	currentLambda = lambda_contract;
	index = bestLambda(lambda_grid, currentLambda);
	tmpBeta = storeBeta.col(index);
	double score_contract = evalLambda(currentLambda, tmpBeta);
	if (score_contract <= score_reflect) {
	  scores(sortOrder(dim)) = score_contract;
	  lambda_grid.col(sortOrder(dim)) = lambda_contract;
	  continue;
	}
      } else if (score_reflect >= scores(sortOrder(dim))) { 
	// contract inside
	arma::colvec lambda_contract = centroid + Beta * 
	  (lambda_grid.col(sortOrder(dim)) - centroid);
	currentLambda = lambda_contract;
	index = bestLambda(lambda_grid, currentLambda);
	tmpBeta = storeBeta.col(index);
	double score_contract = evalLambda(currentLambda, tmpBeta);
	if (score_contract < scores(sortOrder(dim))) {
	  scores(sortOrder(dim)) = score_contract;
	  lambda_grid.col(sortOrder(dim)) = lambda_contract;
	  continue;
	}
      }
      // shrink
      for (unsigned i = 1; i <= dim; i++) {
	lambda_grid.col(sortOrder(i)) = lambda_grid.col(sortOrder(0)) + Delta * 
	  (lambda_grid.col(sortOrder(i)) - lambda_grid.col(sortOrder(0)));
	currentLambda = lambda_grid.col(sortOrder(i));
	index = bestLambda(lambda_grid, currentLambda);
	tmpBeta = storeBeta.col(index);
	scores(sortOrder(i)) = evalLambda(currentLambda, tmpBeta);
	storeBeta.col(sortOrder(i)) = tmpBeta;
      }

    }
    iterSpent = iter * arma::ones<arma::umat> (1, 1);
    currentLambda = lambda_grid.col(minPos);
    index = bestLambda(lambda_grid, currentLambda);
    tmpBeta = storeBeta.col(index);
    scores(0) = evalLambda(currentLambda, tmpBeta);
    beta = arma::zeros<arma::mat> (p, 1);
    beta.col(0) = tmpBeta;
  }

  // procedure for lasso pattern search
  unsigned lps::solveLPS(const arma::colvec& lambda, const arma::colvec& init,
		     arma::colvec& beta_vec)
  {
    // guard agains too small lambda
    if (lambda.max() < sqrt(param.gptol)) return -1;
    double currentLoss = objectiveFunc(init, lambda);
    beta_vec = init;

    double alpha = 1;
    const double alphaIncr = 2.0;
    const double alphaDecr = .8;
    const double alphaMax = 1e8;
    const double alphaMin = 1e-3;

    bool evalFullGrad = 0;
    const unsigned reducedNewton = 500;
    unsigned two_metric_threshold = 
      init.n_rows < reducedNewton ? init.n_rows : reducedNewton;

    // loop to find the minimum
    arma::colvec beta_new = arma::zeros<arma::mat>(p, 1);
    unsigned iter;
    for (iter = 1; iter < param.maxIter; ++iter) {
      double gpnorm = firstOrderMove(beta_vec, beta_new, evalFullGrad,
				     iter, lambda, alpha);
      if (gpnorm < 0)
	return iter;

      arma::colvec beta_step_first (beta_new - beta_vec);
      double sufficientDecreaseMargin = 2 * pow(norm(beta_step_first, 2), 3);

      // find the subspace identified by the first-order step
      arma::uvec nonzeroSet, zeroSet;
      nonzeroSet = find(abs(beta_new) >  param.gptol);
      zeroSet    = find(abs(beta_new) <= param.gptol);

      // check loss for first step
      double loss_first = objectiveFunc(beta_new, lambda);

      if (param.useNewton && nonzeroSet.n_rows <= two_metric_threshold) {
	// Newton step to enhance the procedure
	arma::colvec beta_newton;
	double newtonStepLength = newtonStep(beta_vec, beta_newton, lambda, gpnorm,
					     nonzeroSet, zeroSet);
	// if permissible step length along Newton direction too short
	// don't bother trying it
	if (newtonStepLength <= 1e-2) {
	  currentLoss = loss_first;
	  beta_vec = beta_new;
	} else {
	  // try the truncated Newton step
	  arma::colvec beta_step = beta_newton - beta_vec;
	  beta_newton.elem(nonzeroSet) = beta_vec.elem(nonzeroSet) + newtonStepLength
	    * beta_step.elem(nonzeroSet);

	  double loss_newton = objectiveFunc(beta_newton, lambda);
	  double loss_threshold = currentLoss - 1e-3 * 
	    pow(newtonStepLength * norm(beta_step, 2), 3);
	  if (loss_newton <= loss_first && loss_newton <= loss_threshold) {
	    beta_vec = beta_newton;
	    currentLoss = loss_newton;
	  } else {
	    if (loss_first <= currentLoss - sufficientDecreaseMargin) {
	      beta_vec = beta_new;
	      currentLoss = loss_first;
	    }
	  }
	}
      } else { // no newton step tried
	currentLoss = loss_first;
	beta_vec = beta_new;
      }

      // update loss, since the first-order step was successful,
      // even though the accelerated step failed
      if (loss_first > currentLoss - sufficientDecreaseMargin) {
	// insufficient decrease
	alpha *= alphaIncr;
	if (alpha > alphaMax) {
	  if (param.output) {
	    Rcpp::Rcout << "alpha = " << alpha;
	    Rcpp::Rcout << " alphaMax = " << alphaMax << std::endl;
	    iter *= -1; // indicate not converge exactly
	  }
	  break;
	}
      } else { // take first-order step and derease alpha for next iteration
	alpha *= alphaDecr;
	if (alpha < alphaMin)
	  alpha = alphaMin;
      } // else take first-order step
    }
    return iter;
  }

  double lps::firstOrderMove(arma::colvec& currentBeta,
			     arma::colvec& newBeta,
			     bool& evalFullGrad,
			     unsigned iter,
			     const arma::colvec& lambda,
			     double alpha)
  {
    arma::uvec index;
    // arma::uvec complement;

    if (evalFullGrad) {
      index.reshape(p, 1);
      for (unsigned i = 0; i < p; i++)
	index(i) = i;
    } else {
      // generate the set to evaluate first order steps
      arma::colvec unif = arma::randu<arma::mat> (p - constants.n_rows, 1);
      index = merge(arma::find(abs(currentBeta) > param.gptol),
		    arma::find(unif <= param.sigma));
      index = merge(index, constants);
      // NOTES: complement is not used
    }

    // find gradient for index of current beta
    arma::colvec tmpGradient;
    ptrLoss -> gradient(tmpGradient, currentBeta, index);
    arma::colvec Gradient = arma::zeros<arma::mat> (p, 1);
    Gradient.elem(index) = tmpGradient;
    arma::colvec subGrad  = Gradient;
    // correct the gradient by the subgradient of the penalty
    arma::uvec tmp;
    for (unsigned ptr = 0; ptr < lambda.n_rows; ptr++) {
      unsigned starting = ptr * numRow;
      arma::colvec localBeta = currentBeta.rows(starting, starting + numRow - 1);
      tmp = arma::find(localBeta > param.gptol);
      if (tmp.n_rows > 0)
	subGrad.elem(tmp + starting) += lambda(ptr);
      tmp = arma::find(localBeta < -param.gptol);
      if (tmp.n_rows > 0)
	subGrad.elem(tmp + starting) -= lambda(ptr);
      tmp = find(abs(localBeta) < param.gptol);
      for (unsigned i = 0; i < tmp.n_rows; i++) {
	unsigned position = tmp(i) + starting;
	if (Gradient(position) <= 0)
	  subGrad(position) = std::min(0.0, Gradient(position) + lambda(ptr));
	else
	  subGrad(position) = std::max(0.0, Gradient(position) - lambda(ptr));
      }
    }
    subGrad.elem(constants) = Gradient.elem(constants);
    arma::uvec nonzeroSet0 = find(abs(currentBeta) >  param.gptol);
    arma::uvec zeroSet0    = find(abs(currentBeta) <= param.gptol);

    // norm of subgradient, only used to determine if convergence is achieved
    double gpnorm = arma::norm(subGrad, 2);
    if (param.output)
      print(lambda, iter, gpnorm, alpha, nonzeroSet0.n_rows);
    // check convergence
    if (gpnorm < param.gptol) {
      if (evalFullGrad) {
	if (param.output) {
	  print(lambda, iter, gpnorm, alpha, nonzeroSet0.n_rows, 1);
	  Rcpp::Rcout << "Converged **" << std::endl;
	}
	return -1;
      } else {
	evalFullGrad = 1;
      }
    } else {
      evalFullGrad = 0;
    }

    // evaluate firs-order step
    // NOTES: this can be a virtual function in loss
    newBeta = arma::zeros<arma::mat> (p, 1);
    for (unsigned ptr = 0; ptr < lambda.n_rows; ptr++) {
      unsigned starting = ptr * numRow;
      arma::colvec localBeta = currentBeta.rows(starting, starting + numRow - 1);
      arma::colvec localGrad = Gradient.rows(starting, starting + numRow - 1);
      tmp = find(localGrad - alpha * localBeta < -lambda(ptr));
      if (tmp.n_rows > 0)
	newBeta.elem(tmp + starting) = localBeta.elem(tmp) - 
	  (localGrad.elem(tmp) + lambda(ptr)) / alpha;
      tmp = find(localGrad - alpha * localBeta >  lambda(ptr));
      if (tmp.n_rows > 0)
	newBeta.elem(tmp + starting) = localBeta.elem(tmp) - 
	  (localGrad.elem(tmp) - lambda(ptr)) / alpha;
    }
    newBeta.elem(constants) = currentBeta.elem(constants) - 
      Gradient.elem(constants) / alpha;
    /*
    Rcpp::Rcout << index.n_rows << ' ';
    Rcpp::Rcout << objectiveFunc(newBeta, lambda) << ' ';
    Rcpp::Rcout << objectiveFunc(currentBeta, lambda) << std::endl;
    */
    return gpnorm;
  }

  double lps::newtonStep(arma::colvec& beta_vec,
			 arma::colvec& beta_newton,
			 const arma::colvec& lambda,
			 double gpnorm,
			 const arma::uvec& nonzeroSet,
			 const arma::uvec& zeroSet)
  {
    // calculate reduced gradient for current nonZero set
    arma::colvec tmpGrad;
    // NOTES: gradient should be already evaluated in first order...
    ptrLoss -> gradient(tmpGrad, beta_vec, nonzeroSet);
    arma::colvec reducedGrad = arma::zeros<arma::mat> (p, 1);
    reducedGrad.elem(nonzeroSet) = tmpGrad;
    // store elements in constants
    tmpGrad = reducedGrad.elem(constants);


    arma::uvec tmp;
    for (unsigned ptr = 0; ptr < lambda.n_rows; ptr++) {
      unsigned starting = ptr * numRow;
      arma::colvec localBeta = beta_vec.rows(starting, starting + numRow - 1);
      tmp = find(localBeta > param.gptol);
      if (tmp.n_rows > 0)
	reducedGrad.elem(tmp + starting) += lambda(ptr);
      tmp = find(localBeta < -param.gptol);
      if (tmp.n_rows > 0)
	reducedGrad.elem(tmp + starting) -= lambda(ptr);
    }
    reducedGrad.elem(constants) = tmpGrad;

    arma::mat partialHessian;
    ptrLoss -> hessian(partialHessian, beta_vec, nonzeroSet);
    double scale_factor = sum(partialHessian.diag()) / 
      static_cast<double> (nonzeroSet.n_rows * 4);
    double damping = scale_factor < 10 * gpnorm ? scale_factor : 10 * gpnorm;
    partialHessian += damping * arma::eye(nonzeroSet.n_rows, nonzeroSet.n_rows);

    // NOTES: left inverse of the matrix solve by arma
    //arma::colvec newtonStep = solve(partialHessian, reducedGrad.elem(nonzeroSet));
    arma::colvec newtonStep;
    newtonStep = solve(partialHessian, reducedGrad.elem(nonzeroSet));
    arma::colvec beta_step = arma::zeros<arma::mat> (p, 1);
    beta_step.elem(nonzeroSet) = newtonStep;
    beta_newton = beta_vec - beta_step;
    if (zeroSet.n_rows > 0)
      beta_newton.elem(zeroSet) = arma::zeros<arma::mat> (zeroSet.n_rows, 1);

    // find the indices for which the Newton step causes a change in sign
    // --------------------------------------------------------------
    // If there were sign changes in the full Newton step, rescale
    // the step to stop ar the first kink. Evaluate function and 
    // check for the decrease over the step from the first-order model
    // ---------------------------------------------------------------

    arma::uvec index = find(beta_newton % beta_vec < 0);
    double newtonStepLength;
    if (index.n_rows > 0) {
      arma::colvec tmpRatio = abs(beta_vec.elem(index) / beta_step.elem(index));
      newtonStepLength = tmpRatio.min() < 1 ? tmpRatio.min() : 1;
    } else {
      newtonStepLength = 1;
    }
    return newtonStepLength;
  }
}
