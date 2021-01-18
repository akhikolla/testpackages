#ifndef MVBERNOULLI_H_
#define MVBERNOULLI_H_

#include <vector>
#include "comb.h"
#include "loss.h"

namespace lps {
  // the class to implement multivariate Bernoulli distribution
  class MVBernoulli : public Loss {
  private:
    // order of the outcomes to be considered
    const unsigned K;
    unsigned order;
    // number of columns for beta estimate
    unsigned long numCol;
    // augmented outcome matrix
    arma::mat Y;
    arma::mat augY;
    // function for n choose k
    // resuts to be used in eval, gradient and hessian
    arma::mat fx, eS;
    arma::colvec eb;
    // link table
    std::vector<std::vector<int> > linkTable;
    // inverse of link table
    std::vector<std::vector<int> > invLink;
    // randome effects and number of levels
    arma::mat bi;
    unsigned numLevel;
    arma::uvec Z;
    
    // tool function for n choose k
    static long nchoosek(int n, int k)
    {
      if (k >= n) return 1;
      k = k < n - k ? k : n - k;
      long ret = 1;
      for (int i = 0; i < k; i++)
	ret = ret * (n - i) / (i + 1);
      return ret;
    }
    // calculate number of columns needed for beta
    static unsigned long addUp(int maxOrder, int n_cols)
    {
      unsigned long ret = 0;
      for (int i = 1; i <= maxOrder; i++)
	ret += nchoosek(n_cols, i);
      return ret;
    }
    // bisearch and subset
    static bool biSearch(const std::vector<int>&,
			 int, unsigned, unsigned);

    void augmentY() {
      // augment response Y
      augY = arma::ones<arma::mat>(n, numCol);
      unsigned pos = 0;
      for (unsigned i = 1; i <= order; i++) {
	lps::comb obj(K, i);
	while(!obj.empty()) {
	  std::vector<int>& combs = obj.getNext();
	  for (std::vector<int>::const_iterator iter = combs.begin();
	       iter != combs.end(); iter++) 
	    augY.col(pos) = augY.col(pos) % Y.col(*iter);
	  pos++;
	}
      }
    }

    void setOrder(const unsigned order_)
    {
      order = order_;
      numCol = addUp(order, K);
      augmentY();
      // construct link table
      std::vector<std::vector<int> > link(static_cast<int>(pow(2., static_cast<double>(K))) - 1);
      unsigned pos = 0;
      for (unsigned i = 1; i <= K; i++) {
	lps::comb obj(K, i);
	while(!obj.empty())
	  link[pos++] = obj.getNext();
      }
      linkTable.clear();
      linkTable.resize(link.size());
      for (unsigned i = 0; i < linkTable.size(); i++)
	for (unsigned j = 0; j <= i; j++)
	  if (isSubset(link[j], link[i]))
	    linkTable[i].push_back(j);

      // construct the inverse of link table
      invLink.clear();
      invLink.resize(numCol);
      for (unsigned i = 0; i < numCol; i++) 
	for (unsigned j = i; j < linkTable.size(); j++) 
	  if (biSearch(linkTable[j], i, 0, linkTable[j].size() - 1))
	    invLink[i].push_back(j);

    }

  public:
    explicit MVBernoulli (const arma::mat& inputY,
			  const arma::mat& inputX)
      : Loss(inputX), K(inputY.n_cols), Y(inputY) 
    {
      numLevel = 0;
      setOrder(2);
    }

    virtual void toSetOrder(const unsigned input = 2) {
      setOrder(input);
    }

    static bool isSubset(const std::vector<int>&,
			 const std::vector<int>&);
 

    void setLevel(unsigned in, const arma::uvec& inputZ) {
      numLevel = in;
      Z = inputZ;
    }
    void setbi(const arma::mat& inputB) {
      bi = inputB;
    }
    arma::mat getAug() { return augY;} ;
    inline void link (const arma::mat&, const arma::colvec&);
    virtual double eval(const arma::colvec&);
    virtual void gradient(arma::colvec&,
			  const arma::colvec&,
			  const arma::uvec&);
    virtual void hessian(arma::mat&,
			 const arma::colvec&,
			 const arma::uvec&);
    virtual void mean(arma::colvec&);
    virtual void variance(arma::mat&);
    virtual void addRand(const arma::colvec& err) {
      for (unsigned i = 0; i < K; i++) {
	arma::colvec tmpVec = err.rows(i * n, (i + 1) * n - 1);
	Y.col(i) += tmpVec;
      }
      augmentY();
    }
    virtual arma::colvec getlinear() const {
      arma::colvec ret = arma::zeros<arma::mat> (fx.n_rows * fx.n_cols, 1);
      for (unsigned row = 0; row < fx.n_rows; row++)
	ret.rows(row * fx.n_cols, (row + 1) * fx.n_cols - 1) = trans(fx.row(row));
      return ret;
    }
    virtual unsigned getDim() const {return p * numCol;};
    virtual unsigned getNumCol() const { return numCol; };
    virtual unsigned getK() const {return K; };
    virtual arma::colvec getY() const {
      arma::colvec ret = arma::zeros <arma::colvec>(numCol * n, 1);
      for (unsigned i = 0; i < n; i++)
	ret.rows(i * numCol, (i + 1) * numCol - 1) = trans(augY.row(i));
      return ret;
    }
    virtual ~MVBernoulli() {};
  };
}

#endif // MVBERNOULLI_H_
