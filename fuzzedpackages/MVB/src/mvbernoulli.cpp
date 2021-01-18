#include "mvbernoulli.h"

namespace lps {
  bool MVBernoulli::isSubset(const std::vector<int>& a,
			     const std::vector<int>& b) {
    if (a.size() > b.size()) return 0;
    unsigned posA, posB;
    posA = posB = 0;
    while (posA < a.size()) {
      while (a[posA] > b[posB]) posB++;
      if (a[posA] < b[posB]) return 0;
      posA++;
    }
    return 1;
  }
  
  bool MVBernoulli::biSearch(const std::vector<int>& vec,
			     int target, unsigned low, unsigned high) {
    if (target < vec[low] ||
	target > vec[high]) return 0;
    if (high - low <= 1)
      return target == vec[low] || target == vec[high];
    unsigned mid = (low + high) / 2;
    if (target <= vec[mid]) return biSearch(vec, target, low, mid);
    return biSearch(vec, target, mid, high);
  }

  void MVBernoulli::link(const arma::mat& design,
			 const arma::colvec& beta) {
    checkDim(beta, p * numCol);
    // create a column-connected beta
    arma::mat colBeta = arma::zeros<arma::mat> (p, numCol);
    for (unsigned long j = 0; j < numCol; j++)
      colBeta.col(j) = beta.subvec(j * p, j * p + p - 1);
    fx = design * colBeta;
    // add random component if there is one
    if (numLevel) {
      for (unsigned i = 0; i < n; i++)
	for (unsigned j = 0; j < numCol; j++)
	  fx(i, j) += bi( Z(i), j);
    }   
    eS = arma::zeros<arma::mat> (n, static_cast<int>(pow(2., static_cast<double>(K)) - 1));

    for (unsigned i = 0; i < eS.n_cols; i++) {
      for (unsigned j = 0; j < linkTable[i].size(); j++)
	if (linkTable[i][j] < static_cast<int>(fx.n_cols))
	  eS.col(i) = eS.col(i) + fx.col(linkTable[i][j]);
      eS.col(i) = exp(eS.col(i));
    }
    eb = arma::ones<arma::mat>(n, 1);
    for (unsigned j = 0; j < eS.n_cols; j++)
      eb += eS.col(j);
  }

  double MVBernoulli::eval(const arma::colvec& beta)
  {
    link (X, beta);
    arma::colvec ret = arma::zeros<arma::mat> (n, 1);
    for (unsigned i = 0; i < fx.n_cols; i++)
      ret += fx.col(i) % augY.col(i);
    return -arma::sum(arma::sum(ret - log(eb) )) / static_cast<double> (n);
  }

  void MVBernoulli::gradient(arma::colvec& ret,
			     const arma::colvec& beta,
			     const arma::uvec& index)
  {
    link(X, beta);
    // use a vector to record the number of elements in index
    arma::uvec counts = arma::zeros<arma::umat> (index.n_rows, 1);
    ret.reshape(index.n_rows, 1);
    for (unsigned col = 0; col < numCol; col++) {
      arma::uvec tmp = find(index / p == col);
      if (tmp.n_rows == 0) continue;
      arma::uvec tmp0(tmp.n_rows, 1);
      for (unsigned i = 0; i < tmp.n_rows; i++)
	tmp0(i) = index(tmp(i)) % p;
      arma::mat subX = subMatrix(X, tmp0);
      arma::colvec tmpVec = arma::zeros<arma::mat> (n, 1);
      for (unsigned i = 0; i < invLink[col].size(); i++) 
	tmpVec += eS.col(invLink[col][i]); // add the corresponding colum of eS
      tmpVec /= eb; 
      ret.elem(tmp) = trans(subX) * (tmpVec - augY.col(col));
    }
    ret /= static_cast<double> (n);
  }

  void MVBernoulli::hessian(arma::mat& ret,
			    const arma::colvec& beta,
			    const arma::uvec& index)
  {
    ret.reshape(index.n_rows, index.n_rows);
    for (unsigned i = 0; i < index.n_rows; i++)
      for (unsigned j = i; j < index.n_rows; j++) {
	unsigned pos1 = index(i) % p;
	unsigned pos2 = index(j) % p;
	unsigned col1 = index(i) / p;
	unsigned col2 = index(j) / p;
	arma::colvec tmpVec1 = arma::zeros<arma::mat> (n, 1);
	for (unsigned k = 0; k < invLink[col1].size(); k++) 
	  tmpVec1 += eS.col(invLink[col1][k]);
	arma::colvec tmpVec2 = arma::zeros<arma::mat> (n, 1);
	for (unsigned k = 0; k < invLink[col2].size(); k++)
	  tmpVec2 += eS.col(invLink[col2][k]);
	arma::colvec tmpVec = arma::zeros<arma::mat> (n, 1);
	for (unsigned k = 0; k < eS.n_cols; k++)
	  if (biSearch(linkTable[k], col1, 0, linkTable[k].size() - 1) &&
	      biSearch(linkTable[k], col2, 0, linkTable[k].size() - 1) ) 
	    tmpVec += eS.col(k);
	arma::colvec res = exp(log(tmpVec) - log(eb));
	res -= exp(log(tmpVec1) + log(tmpVec2) - 2.0 * log(eb));
	res = res % X.col(pos1) % X.col(pos2);
	ret(i, j) = arma::sum(res) / static_cast<double> (n);
	if (i != j)
	  ret(j, i) = ret(i, j);
      } // for (j)
  }

  // calculate mean vector
  void MVBernoulli::mean(arma::colvec& ret) {
    ret.reshape(n * numCol, 1);
    for (unsigned col = 0; col < numCol; col++) {
      arma::colvec tmpVec = arma::zeros<arma::mat> (n, 1);
      for (unsigned i = 0; i < invLink[col].size(); ++i)
	tmpVec += eS.col(invLink[col][i]);
      tmpVec /= eb;
      for (unsigned i = 0; i < n; ++i)
	ret(i * numCol + col) = tmpVec(i);
    }
  }

  // calculate W (variance) for tuning
  void MVBernoulli::variance(arma::mat& ret)
  {
    ret.reshape(numCol * n, numCol * n);
    for (unsigned col1 = 0; col1 < numCol; col1++) 
      for (unsigned col2 = col1; col2 < numCol; col2++) {
	arma::colvec tmpVec1 = arma::zeros<arma::mat> (n, 1);
	for (unsigned k = 0; k < invLink[col1].size(); k++) 
	  tmpVec1 += eS.col(invLink[col1][k]);
	arma::colvec tmpVec2 = arma::zeros<arma::mat> (n, 1);
	for (unsigned k = 0; k < invLink[col2].size(); k++)
	  tmpVec2 += eS.col(invLink[col2][k]);
	arma::colvec tmpVec = arma::zeros<arma::mat> (n, 1);
	for (unsigned k = 0; k < eS.n_cols; k++)
	  if (biSearch(linkTable[k], col1, 0, linkTable[k].size() - 1) &&
	      biSearch(linkTable[k], col2, 0, linkTable[k].size() - 1) )
	    tmpVec += eS.col(k);
	arma::colvec res = exp(log(tmpVec) - log(eb));
	res -= exp(log(tmpVec1) + log(tmpVec2) - 2.0 * log(eb));
	for (unsigned k = 0; k < n; k++) {
	  unsigned row = numCol * k + col1;
	  unsigned col = numCol * k + col2;
	  ret(row, col) = res(k);
	  if (row != col)
	    ret(col, row) = res(k);
	}
      }
  }
  
}
