#ifndef _DISCRETE_CD_H_
#define _DISCRETE_CD_H_

#include "type.h"
#include <vector>

void NewtonIterUpdate(const int& nh, const int& di, const int& rj, const int& nrows, const Eigen::MatrixXd& dm, const Eigen::MatrixXd& ym, const Eigen::MatrixXi& nzIndex,
                      Eigen::MatrixXd& logitm, Eigen::MatrixXd& beta, bool& isBetaZero, const double& hval, const double& penalty, const double& qtol);
void NewtonIterTmp(const int& nh, const int& di, const int& rj, const int& nrows, const Eigen::MatrixXd& dm, const Eigen::MatrixXd& ym, const Eigen::MatrixXi& nzIndex,
                  const Eigen::MatrixXd& logitm, const Eigen::MatrixXd& beta, const bool& isBetaZero, const double& hval, const double& penalty,
				  const double& qtol, Eigen::MatrixXd& logitm_Tmp, Eigen::MatrixXd& beta_Tmp, bool& isBetaZero_Tmp, bool& needUpdate);
void NewtonIterTmp_GD(const int& nh, const int& di, const int& rj, const int& nrows, const Eigen::MatrixXd& dm, const Eigen::MatrixXd& ym, const Eigen::MatrixXi& nzIndex,
                      const Eigen::MatrixXd& logitm, const Eigen::MatrixXd& beta, const bool& isBetaZero, const double& hval, const double& penalty,
                      const double& qtol, Eigen::MatrixXd& grad, Eigen::MatrixXd& dif, Eigen::MatrixXd& logitm_Tmp, Eigen::MatrixXd& beta_Tmp, bool& isBetaZero_Tmp, bool& needUpdate);
void InterceptUpdate(const int& nh, const Eigen::MatrixXd& ym, Eigen::MatrixXd& logitm, Eigen::MatrixXd& beta, const double& hval, const double& qtol);
void dmFetch(Eigen::MatrixXd& dmt, Eigen::MatrixXi& nzIndt, int& rCount, const Eigen::MatrixXd& dm, const int& nrow, const int& ncol, const Eigen::VectorXi& rIn, const Eigen::VectorXi& cIn,
             const Eigen::VectorXd& scaling);
void firstDMFetch(Eigen::MatrixXd& dmt, Eigen::MatrixXi& nzIndt, int& rCount, const Eigen::MatrixXd& dm, const int& nrow, const int& ncol, const Eigen::VectorXi& rIn, const Eigen::VectorXi& cIn,
				 Eigen::VectorXd& scaling);
void OneCDLoop(const int& node, Eigen::MatrixXi& G, const int& eor_nr, const Eigen::MatrixXi& eor, std::vector<int>& active, double& MAD, const Eigen::VectorXi& nobsVec,
               const Eigen::MatrixXi& ndfs, const VectorXMXd& dM, const VectorXMXd& yM, const VectorXVXi& yNZIndex, VectorXMXd& logitM, MatrixXMXd& betaM,
			   MatrixXb& IsBetaZeros, const Eigen::MatrixXd& hvals, const Eigen::MatrixXd& penalties, const double& qtol, const VectorXVXi& obsIndex,
			   const MatrixXVXi& levelIndex, const MatrixXVXd& scales, const int& maxRows, const int& maxCols);
void CDOnePoint(const int& node, Eigen::MatrixXi& G, const int& eor_nr, const Eigen::MatrixXi& eor, const double& eps, const Eigen::VectorXi& nobsVec,
               const Eigen::MatrixXi& ndfs, const VectorXMXd& dM, const VectorXMXd& yM, const VectorXVXi& yNZIndex, VectorXMXd& logitM, MatrixXMXd& betaM,
			   MatrixXb& IsBetaZeros, const Eigen::MatrixXd& hvals, const Eigen::MatrixXd& penalties, const double& qtol, const VectorXVXi& obsIndex,
			   const MatrixXVXi& levelIndex, const MatrixXVXd& scales, const int& maxRows, const int& maxCols);
void CDAlgo(int node, int dataSize, const Eigen::MatrixXi& data, const Eigen::VectorXi& nlevels, const VectorXVXi& obsIndex, const MatrixXVXi& levelIndex,
			int eor_nr, const Eigen::MatrixXi& eor, int nlam, double eps, double convLb, double qtol, Eigen::VectorXd& lambdaSeq, Eigen::VectorXd& log_like, Eigen::VectorXd& dur,
            MatrixXMXd& betaM, MatrixXMXd& betaN, Eigen::MatrixXi& estimateG, Eigen::MatrixXd& weights, double gamma, double upperbound = -1.0, int threshold = 3);
void maxLambda(int node, int dataSize, const Eigen::MatrixXi& data, const Eigen::VectorXi& nlevels, const VectorXVXi& obsIndex, const MatrixXVXi& levelIndex, MatrixXMXd& betaM, Eigen::MatrixXd& weights, double& lambda, double gamma, double upperbound = -1.0);
// void CD_learning(int node, int dataSize, const Eigen::MatrixXi& data, const Eigen::VectorXi& nlevels, const VectorXVXi& obsIndex, const MatrixXVXi& levelIndex,
                 // int eor_nr, const Eigen::MatrixXi& eor, double eps, double convLb, double qtol, double& log_like, MatrixXMXd& betaM, Eigen::MatrixXd& adaptWeights);

#endif
