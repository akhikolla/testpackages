#ifndef SEQUENTIALR_H_
#define SEQUENTIALR_H_
#include "distribution.h"

class SequentialR : virtual public distribution, Logistic, Normal, Cauchit, Student, Gumbel, Gompertz{
public:
  SequentialR();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gompertz(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gompertz(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gumbel(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gumbel(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const;
  virtual Eigen::MatrixXd inverse_derivative_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const ;

  // List GLMseq(std::string response,
  //             StringVector explanatory_complete,
  //             StringVector explanatory_proportional,
  //             std::string distribution,
  //             SEXP categories_order,
  //             DataFrame dataframe);


};

#endif
