#include "mvbfit.h"

void collectFIT(Rcpp::List& res, lps::lps& obj)
{
  res["beta"] = obj.getCoef();
  res["iters"] = obj.getIter();
}

void collectLPS(Rcpp::List& res, lps::lps& obj) 
{
  collectFIT(res, obj);
  res["score"] = obj.getScore();
}

int getMethod (const std::string& str) {
  if (str == "AIC") 
    return AIC;
  else if (str == "BIC")
    return BIC;
  else if (str == "GACV")
    return GACV;
  else if (str == "BGACV")
    return BGACV;
  else
    return -1;
}

RcppExport SEXP get_eS(SEXP fx_,
		       SEXP K_) {
  try {
    Rcpp::NumericMatrix fxr(fx_);
    int n = fxr.nrow(), givenCol = fxr.ncol();
    unsigned K = INTEGER(K_)[0];

    arma::mat fx(fxr.begin(), n, givenCol, false);
    arma::mat eS = arma::zeros<arma::mat> (n, static_cast<int>(pow(2., static_cast<double>(K)) - 1));

    // link table used to construct eS
    std::vector<std::vector<int> > link(static_cast<int>(pow(2., static_cast<double>(K)) - 1));
    std::vector<std::vector<int> > linkTable(static_cast<int>(pow(2., static_cast<double>(K)) - 1));
    unsigned pos = 0;
    for (unsigned i = 1; i <= K; i++) {
      lps::comb obj(K, i);
      while (!obj.empty())
	link[pos++] = obj.getNext();
    }
    for (unsigned i = 0; i < linkTable.size(); i++)
      for (unsigned j = 0; j <= i; j++)
	if (lps::MVBernoulli::isSubset(link[j], link[i]))
	  linkTable[i].push_back(j);

    for (unsigned i = 0; i < eS.n_cols; i++) {
      for (unsigned j = 0; j < linkTable[i].size(); j++)
	if (linkTable[i][j] < static_cast<int>(fx.n_cols))
	  eS.col(i) = eS.col(i) + fx.col(linkTable[i][j]);
      eS.col(i) = exp(eS.col(i));
    }

    return Rcpp::wrap(eS);
  } catch ( std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

RcppExport SEXP loglike(SEXP inputX,
			SEXP inputY,
			SEXP inputBeta,
			SEXP family) {
  try {
    // inputY should be a matrix for multivariate case
    Rcpp::NumericMatrix Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    Rcpp::NumericVector betar(inputBeta);
    int n = Xr.nrow(), p = Xr.ncol(), K;
    Rcpp::StringVector fam(family);
    std::string distribution(fam[0]);
    if (distribution == "mvbernoulli")
      K = Yr.ncol();
    else
      K = 1;

    arma::mat Ya(Yr.begin(), n, K, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    lps::Loss* lossObj = lps::DistriFactory::instance().createLoss(distribution, Ya, Xa);
    arma::colvec beta_vec(betar.begin(), p * lossObj -> getNumCol(), false);

    return Rcpp::wrap(lossObj -> eval(beta_vec));
  } catch ( std::exception &ex ) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall

}

RcppExport SEXP mvbfit(SEXP inputX,
		       SEXP inputY,
		       SEXP maxOrder,
		       SEXP output,
		       SEXP printIter) {
  try {
    // inputY should be a matrix for multivariate case
    Rcpp::NumericMatrix Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    int n = Xr.nrow(), p = Xr.ncol(), K = Yr.ncol();
    std::string distribution("mvbernoulli");

    arma::mat Ya(Yr.begin(), n, K, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    
    lps::lps lpsObj(distribution, Ya, Xa);
    if (INTEGER(maxOrder)[0] != 2)
      lpsObj.setOrder(INTEGER(maxOrder)[0]);
    lpsObj.setParam(INTEGER(output)[0], INTEGER(printIter)[0]);
    lpsObj.runNewton();
    Rcpp::List res;
    res["beta"] = lpsObj.getCoef();
    res["iters"] = lpsObj.getIter();
    return res;
  } catch ( std::exception &ex ) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }

  return R_NilValue; //-Wall
}

RcppExport SEXP mvblps(SEXP inputX,
		       SEXP inputY,
		       SEXP maxOrder,
		       SEXP lambda,
		       SEXP output,
		       SEXP printIter,
		       SEXP tuneMethod,
		       SEXP searchMethod) {
  try {
    // inputY should be a matrix for multivariate case
    Rcpp::NumericMatrix Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    Rcpp::NumericMatrix Lr(lambda);
    unsigned n = Xr.nrow(), p = Xr.ncol(), K = Yr.ncol();
    std::string distribution("mvbernoulli");

    arma::mat Ya(Yr.begin(), n, K, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    arma::mat lambda_grid(Lr.begin(), Lr.nrow(), Lr.ncol(), false);
    
    lps::lps lpsObj(distribution, Ya, Xa);
    if (INTEGER(maxOrder)[0] != 2)
      lpsObj.setOrder(INTEGER(maxOrder)[0]);
    lpsObj.setParam(INTEGER(output)[0], INTEGER(printIter)[0]);
    arma::uvec index = arma::zeros<arma::uvec> (K, 1);
    for (unsigned i = 1; i < K; i++)
      index(i) = i * p;
    lpsObj.setConst(index);
    Rcpp::StringVector tune_(tuneMethod);
    std::string tune(tune_[0]);
    Rcpp::StringVector search_(searchMethod);
    std::string search(search_[0]);
    lpsObj.setTune(getMethod(tune));
    if (search == "grid") {
      lpsObj.gridSearch(lambda_grid);
    } else {
      arma::colvec currentLambda = lambda_grid.col(0);
      lpsObj.nelderMead(currentLambda);
    }
    Rcpp::List res;
    collectLPS(res, lpsObj);
    return res;
  } catch ( std::exception &ex ) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }

  return R_NilValue; //-Wall
}

RcppExport SEXP mvbme(SEXP inputX,
		      SEXP inputY,
		      SEXP inputZ,
		      SEXP maxOrder,
		      SEXP output,
		      SEXP printIter) {
  try {
    // inputY should be a matrix for multivariate case
    Rcpp::NumericMatrix Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    Rcpp::NumericVector Zr(inputZ);
    int n = Xr.nrow(), p = Xr.ncol(), K = Yr.ncol();
    std::string distribution("mvbernoulli");

    arma::mat Ya(Yr.begin(), n, K, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    arma::uvec Za = arma::zeros <arma::uvec> (n, 1);
    for (int i = 0; i < n; i++)
      Za(i) = static_cast <unsigned> (Zr(i) - 1);

    lps::gme obj (Ya, Xa, Za);
    obj.setOrder(INTEGER(maxOrder)[0]);
    obj.fit();

    Rcpp::List res;
    res["beta"] = obj.getCoef();
    res["iters"] = obj.getIter();
    res["sigma"] = obj.getSigma();
    res["b"] = obj.getB();
    return res;
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }

  return R_NilValue; //-Wall
}

RcppExport SEXP unifit(SEXP inputX,
		       SEXP inputY,
		       SEXP family,
		       SEXP output) {
  try {
    // only consider vector form
    Rcpp::NumericVector Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    int n = Xr.nrow(), p = Xr.ncol();
    Rcpp::StringVector fam(family);
    std::string distribution(fam[0]);

    arma::colvec Ya(Yr.begin(), n, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    
    lps::lps lpsObj(distribution, Ya, Xa);
    lpsObj.setParam(INTEGER(output)[0]);
    lpsObj.runNewton();
    Rcpp::List res;
    res["beta"] = lpsObj.getCoef();
    res["iters"] = lpsObj.getIter();
    return res;
  } catch ( std::exception &ex ) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

RcppExport SEXP unilps(SEXP inputX,
		       SEXP inputY,
		       SEXP lambda,
		       SEXP family,
		       SEXP output,
		       SEXP tuneMethod) {
  try {
    // only consider vector form
    Rcpp::NumericVector Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    Rcpp::NumericMatrix Lr(lambda);
    int n = Xr.nrow(), p = Xr.ncol();
    Rcpp::StringVector fam(family);
    std::string distribution(fam[0]);

    arma::colvec Ya(Yr.begin(), n, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    arma::mat lambda_grid(Lr.begin(), Lr.nrow(), Lr.ncol(), false);
    
    lps::lps lpsObj(distribution, Ya, Xa);
    lpsObj.setParam(INTEGER(output)[0]);
    arma::uvec index = arma::zeros<arma::uvec>(1, 1);
    lpsObj.setConst(index);
    Rcpp::StringVector tune_(tuneMethod);
    std::string tune(tune_[0]);
    lpsObj.setTune(getMethod(tune));
    lpsObj.gridSearch(lambda_grid);
    Rcpp::List res;
    collectLPS(res, lpsObj);
    return res;
  } catch ( std::exception &ex ) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall

}

RcppExport SEXP stepfit(SEXP inputX,
			SEXP inputY,
			SEXP maxOrder,
			SEXP direction,
			SEXP tuneMethod,
			SEXP output,
			SEXP initial) {
  try {
    // inputY should be a matrix for multivariate case
    Rcpp::NumericMatrix Yr(inputY);
    Rcpp::NumericMatrix Xr(inputX);
    int n = Xr.nrow(), p = Xr.ncol(), K = Yr.ncol();
    std::string distribution("mvbernoulli");

    arma::mat Ya(Yr.begin(), n, K, false);
    arma::mat Xa(Xr.begin(), n, p, false);
    Rcpp::StringVector tune_(tuneMethod);
    std::string tune(tune_[0]);
    Rcpp::StringVector direct_(direction);
    std::string direct(direct_[0]);

    lps::lps lpsObj(distribution, Ya, Xa);
    if (INTEGER(maxOrder)[0] != 2)
      lpsObj.setOrder(INTEGER(maxOrder)[0]);    
    lpsObj.setParam(INTEGER(output)[0]);

    arma::uvec index;
    bool predefine = !(LENGTH(initial) == 1 && REAL(initial)[0] == 0);
    if (predefine) {
      Rcpp::NumericVector ir(initial);
      arma::colvec init(ir.begin(), ir.size(), false);
      index = find(abs(init) > 1e-6);
    }

    if (direct == "forward") {
      if (!predefine) 
	index.reshape(0, 1); // null model
      lpsObj.setConst(index);
      lpsObj.runNewton();
      double current = lpsObj.tune(lpsObj.getCoef(), getMethod(tune));
      // loop to forward addition
      while(1) {
	if (index.size() == lpsObj.getDim()) //  already full model
	  break;
	arma::uvec tmp(index.size() + 1);
	unsigned minPos = 0;
	double minScore = 1e16;
	for (unsigned i = 0; i < lpsObj.getDim(); i++) {
	  arma::uvec loc = find(index == i);
	  if (loc.size() > 0) continue;
	  if (index.size() > 0)
	    tmp.rows(0, index.size() - 1) = index;
	  tmp(index.size()) = i;
	  tmp = sort(tmp);
	  lpsObj.setConst(tmp);
	  lpsObj.runNewton();
	  double score = lpsObj.tune(lpsObj.getCoef(), getMethod(tune));
	  if (score < minScore) {
	    minPos = i;
	    minScore = score;
	  }
	}
	if (minScore < current) {
	  arma::uvec tmp(index.size() + 1);
	  if(index.size() > 0)
	    tmp.rows(0, index.size() - 1) = index;
	  tmp(index.size()) = minPos;
	  index = sort(tmp);
	  current = minScore;
	} else {
	  break;
	}
      }
    } else {
      if (!predefine) 
	index.reshape(lpsObj.getDim(), 1);
      for (unsigned i = 0; i < index.size(); i++)
	index(i) = i;
      lpsObj.setConst(index);
      lpsObj.runNewton();
      double current = lpsObj.tune(lpsObj.getCoef(), getMethod(tune));
      // loop to backward elimination
      while (1) {
	arma::uvec tmp(index.size() - 1);
	unsigned minPos = 0;
	double minScore = 1e16;
	for (unsigned i = 0; i < index.size(); i++) {
	  unsigned pos = 0;
	  for (arma::uvec::iterator it = index.begin();
	       it != index.end(); it++)
	    if (*it != index(i))
	      tmp(pos++) = *it;
	  lpsObj.setConst(tmp);
	  lpsObj.runNewton();
	  double score = lpsObj.tune(lpsObj.getCoef(), getMethod(tune));
	  if (score < minScore) {
	    minPos = i;
	    minScore = score;
	  }
	}
	
	if (minScore < current) {
	  unsigned pos = 0;
	  for (unsigned i = 0; i < index.size(); i++)
	    if (i != minPos)
	      tmp(pos++) = index(i);
	  index = tmp;
	  current = minScore;
	} else {
	  break;
	}
      }
    }

    lpsObj.setConst(index);
    lpsObj.runNewton();

    Rcpp::List res;
    collectFIT(res, lpsObj);
    return res;
  } catch ( std::exception &ex ) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }

  return R_NilValue; //-Wal
}		       

