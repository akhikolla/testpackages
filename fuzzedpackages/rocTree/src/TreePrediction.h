#ifndef TreePrediction_h
#define TreePrediction_h

#include "Data2.h"

class TreePrediction {
public:

  TreePrediction(const Data2* dat2,
                 const arma::uvec& vars,
                 const arma::uvec& values,
                 const arma::uvec& lcs,
                 const arma::uvec& rcs,
                 const arma::uvec& il);

  TreePrediction(const arma::umat& zy,
                 const arma::field<arma::umat>& zt,
                 const arma::uvec& vars,
                 const arma::uvec& values,
                 const arma::uvec& lcs,
                 const arma::uvec& rcs,
                 const arma::uvec& il);

  arma::vec getSurvival(const arma::umat& zt2,
                        const arma::vec& y,
                        const arma::uvec& e,
                        const arma::uvec& vars,
                        const arma::uvec& values,
                        const arma::uvec& lcs,
                        const arma::uvec& rcs,
                        const arma::uvec& il);

  static arma::vec getSurvival(const arma::umat& zt2,
                               const arma::vec& y,
                               const arma::uvec& e,
                               const arma::umat& nodeSize0,
                               const arma::uvec& nodeLabel0,
                               const arma::uvec& tnd30,
                               const arma::umat& treeMat);

  arma::uword getHazard();

  static void transformZ(const arma::mat& z,
                         arma::umat& z2,
                         const arma::mat& dat,
                         const arma::uvec& e,
                         const arma::vec& breaks,
                         const arma::uvec& disc);

  static void transformZH(const arma::mat& z,  // z on tg, col number is the length of tg
                          const arma::vec& tg, //grid of time point
                          arma::umat& z2,
                          const arma::mat& dat,
                          const arma::vec& y,
                          const arma::uvec& e,
                          const arma::vec& breaks,
                          const arma::uvec& disc);

  static arma::vec getHazard(const arma::umat& ztvec,
                             const arma::vec& tg,
                             const arma::vec& y,
                             const arma::uvec& e,
                             const arma::mat& fy2,
                             const double h,
                             const arma::umat& nodeSize0,
                             const arma::uvec& nodeLabel0,
                             const arma::uvec& tnd30,
                             const arma::umat& treeMat);


  const arma::umat& get_nodeSize() const
  {
    return nodeSize;
  };
  const arma::uvec& get_nodeLabel() const
  {
    return nodeLabel;
  };
  const arma::uvec& get_nodeMap() const
  {
    return tnd3;
  };

private:
  // We only need Z(t) at uncensored time points, t\in T
  // matrix of node size, each row is a terminal node, each column is a t

  arma::umat nodeSize;
  // vector of node label of ZY in dat2

  arma::uvec nodeLabel;

  // vector that defines a map between node number and its location in matrix nodeSize
  arma::uvec tnd3;


  // for forest, need a vector of id of I_2b in the orignal sample
};

#endif /* TreePrediction_h */
