#ifndef ForestPrediction_h
#define ForestPrediction_h


#include "Data.h"
#include "Data2.h"
#include "Tree.h"
#include "Forest.h"
#include <memory>

class ForestPrediction {
public:
   ForestPrediction(const Data2* dat2,
                    const arma::umat& ids,
                    const std::vector<std::shared_ptr<Tree> >& trees,
                    arma::uword n);

   ForestPrediction(const arma::umat& zy,
                    const arma::field<arma::umat>& zt,
                    const arma::umat& ids,
                    const std::vector<std::shared_ptr<Tree> >& trees,
                    arma::uword n);

   static void transformZ(const arma::mat& z,
                          arma::umat& z2,
                          const arma::mat& matX,
                          const arma::uvec& e,
                          const arma::vec& breaks,
                          const arma::uvec& disc);

   static void transformZH(const arma::mat& z,  // z on tg
                           const arma::vec& tg, //grid of time point
                           arma::umat& z2,
                           const arma::mat& dat, // dat is predictor matrix
                           const arma::vec& y,
                           const arma::uvec& e,
                           const arma::vec& breaks,
                           const arma::uvec& disc);

   static void transformZ0(const arma::mat& z,
                           arma::umat& z2,
                           const arma::mat& matX,
                           const arma::uvec& e,
                           const arma::vec& breaks,
                           const arma::uvec& disc);

   arma::vec getSurvival(const arma::umat& zt2,
                         const Data2* dat2,
                         const arma::umat& ids,
                         const std::vector<std::shared_ptr<Tree> >& trees);

   static arma::vec getSurvival(const arma::umat& zt2,
                                const arma::vec& y,
                                const arma::uvec& e,
                                const arma::field<arma::uvec>&& nodeSizeB0,
                                const arma::umat&& nodeLabelB0,
                                const arma::field<arma::uvec>&& tnd3B0,
                                const arma::umat&& ids,
                                const arma::field<arma::umat>&& trees);

   static arma::vec getHazard(const arma::umat& ztvec,
                              const arma::vec& tg,
                              const arma::vec& y,
                              const arma::uvec& e,
                              const arma::mat& Kmat,
                              const double h,
                              const arma::field<arma::uvec>&& nodeSizeB0,
                              const arma::umat&& nodeLabelB0,
                              const arma::field<arma::uvec>&& tnd3B0,
                              const arma::umat&& ids,
                              const arma::field<arma::umat>&& trees);

   static arma::vec getSurvival2(const arma::umat& zt2,
                                 const arma::vec& y,
                                 const arma::uvec& e,
                                 const arma::field<arma::uvec>&& nodeSizeB0,
                                 const arma::umat&& nodeLabelB0,
                                 const arma::field<arma::uvec>&& tnd3B0,
                                 const arma::umat&& ids,
                                 const arma::field<arma::umat>&& trees);

   const arma::field<arma::uvec>& get_nodeSize() const
   {
      return nodeSizeB;
   };
   const arma::umat& get_nodeLabel() const
   {
      return nodeLabelB;
   };
   const arma::field<arma::uvec>& get_nodeMap() const
   {
      return tnd3B;
   };
private:
   arma::field<arma::uvec> nodeSizeB;

   arma::umat nodeLabelB; // Each column is for one tree

   // vector that defines a map between node number and its location in matrix nodeSize
   arma::field<arma::uvec> tnd3B;

};

#endif /* ForestPrediction_h */
