#ifndef Tree_h
#define Tree_h

#include <vector>
#include <random>
//#include <iostream>
//#include <armadillo>

#include "Data.h"

typedef unsigned int uint;

class Tree {
public:

  Tree(arma::uvec&& lc,
       arma::uvec&& rc,
       arma::uvec&& svars,
       arma::uvec&& svals,
       arma::uvec&& isl,
       arma::uvec&& pr) : left_childs(lc), right_childs(rc),
       split_vars(svars), split_values(svals), isLeaf(isl), parents(pr) { };

  Tree(arma::uvec&& lc,
       arma::uvec&& rc,
       arma::uvec&& svars,
       arma::uvec&& svals,
       arma::uvec&& isl) : left_childs(lc), right_childs(rc),
       split_vars(svars), split_values(svals), isLeaf(isl) { };

  Tree();

  // Based on the data preparation step, the covariate values are integers
  //arma::uvec get_split_vars() const;
  const arma::uvec& get_split_values() const;
  const arma::uvec& get_left_childs() const;
  const arma::uvec& get_right_childs() const;
  const arma::uvec& get_isLeaf() const;
  const arma::uvec& get_parents() const;
  const arma::uvec& get_split_vars() const;

  // cut the current large tree to be a smaller tree whose terminal node set is given
  void cut(arma::uvec& nodeTerminal);
  // set the descendant nodes to zero/null
  void setzero(uint i, uint ndcount);


  // fmat: \hat f*(t,\tau), each column is one node, each row is one time point
  // Smat: \hat S*(t,\tau)
  // isLeafTemp: 1 if the node is terminal, 0 otherwise
  // can be used in splitting based ICON
  static double get_ICONTrain(const arma::uvec& isLeafTemp,
                              const arma::mat& fmat,
                              const arma::umat& Smat);


  // iconAll is the ICON value of the trees whose sizes ranging from 1 to numLeaf
  // nodeSetList gives the terminal nodes of these trees
  // The function uses fmat and Smat to modify iconAll and nodeSetList, does not return anything
  void findOptimalSizekSubtree(arma::mat& fmat,
                               arma::umat& Smat,
                               arma::vec& iconAll,
                               arma::field<arma::uvec>& nodeSetList,
                               uint numLeaf);

  // find a sequence of beta as the tuning grid for CV
  // based on the sequence of optimal size k subtrees
  // use iconAll to modify beta and sizeTree
  static void findBeta(arma::vec& iconAll,
                       arma::vec& beta,
                       arma::uvec& sizeTree);


private:
  arma::uvec left_childs;
  arma::uvec right_childs;
  arma::uvec split_vars;
  arma::uvec split_values;
  arma::uvec isLeaf;
  arma::uvec parents;
};

#endif /* Tree_h */

