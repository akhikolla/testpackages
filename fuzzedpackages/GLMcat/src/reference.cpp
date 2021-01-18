#include "distribution.h"
#include "reference.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;

// #include <boost/algorithm/string.hpp>


// [[Rcpp::depends(RcppEigen)]]
ReferenceF::ReferenceF(void) {
  distribution dist;
}

Eigen::VectorXd ReferenceF::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_logit( eta(j) ) / ( 1-
      std::max(1e-10, std::min(1-1e-6,cdf_logit( eta(j) )))
    );


    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::VectorXd ReferenceF::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_normal( eta(j) ) / ( 1-
      std::max(1e-10, std::min(1-1e-6,cdf_normal( eta(j) )))
    );
    norm1 += pi[j];

  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_logistic(const Eigen::VectorXd& eta2 ) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_logistic(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_logit( eta2(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,cdf_logit(eta2(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-cdf_logit(eta2(j)))) ); }

  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::MatrixXd ReferenceF::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = ReferenceF::inverse_normal(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_normal( eta(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,cdf_normal(eta(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-cdf_normal(eta(j)))) ); }
  return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    pi[j] = Cauchit::cdf_cauchit( eta(j) ) / ( 1-Cauchit::cdf_cauchit( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_cauchit(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_cauchit(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_cauchit( eta2(j) ) /
    (Cauchit::cdf_cauchit(eta2(j)) * (1-Cauchit::cdf_cauchit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_gompertz( eta(j) ) / ( 1-cdf_gompertz( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_gompertz(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_gompertz(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_gompertz( eta2(j) ) /
    (cdf_gompertz(eta2(j)) * (1-cdf_gompertz(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_gumbel( eta(j) ) / ( 1-cdf_gumbel( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_gumbel(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_gumbel(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_gumbel( eta2(j) ) /
    (cdf_gumbel(eta2(j)) * (1-cdf_gumbel(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    double num = Student::cdf_student(eta(j),freedom_degrees);
    double den = std::max(1e-10, std::min(1-1e-6, 1 - Student::cdf_student(eta(j),freedom_degrees)));
    pi[j] = (num / den);
    norm1 += pi[j];
  }
  return (pi/norm1);
}


Eigen::MatrixXd ReferenceF::inverse_derivative_student(const Eigen::VectorXd& eta2, const double& freedom_degrees) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_student(eta2, freedom_degrees);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  {
    double num = Student::pdf_student( eta2(j) , freedom_degrees);
    double den1 = Student::cdf_student(eta2(j), freedom_degrees) ;
    double den2 = 1-Student::cdf_student(eta2(j),freedom_degrees) ;
    D1(j,j) = (num / std::max(1e-10, std::min(1-1e-6, (den1 * den2)) ));
  }
  Eigen::MatrixXd D3 = pi1 * (pi1.transpose());
  Eigen::MatrixXd D2 = Eigen::MatrixXd(pi1.asDiagonal());
  Eigen::MatrixXd FINAL = D1 * ( D2 - D3 );
  return FINAL;
}


// RCPP_MODULE(referencemodule){
  // Rcpp::function("GLMref", &GLMref,
  //                List::create(_["formula"],
  //                             _["categories_order"],
  //                             _["proportional"] = CharacterVector::create(NA_STRING),
  //                             _["data"],
  //                             _["distribution"] = "logistic",
  //                             _["freedom_degrees"] = 1),
  //                             "Reference model");

  // Rcpp::function("Discrete_CM", &Discrete_CM,
  //                List::create(_["formula"] = R_NaN,
  //                             _["case_id"] = "a",
  //                             _["alternatives"] = "a",
  //                             _["reference"] = R_NaN,
  //                             _["alternative_specific"] = CharacterVector::create( NA_STRING),
  //                             _["data"] = NumericVector::create( 1, NA_REAL, R_NaN, R_PosInf, R_NegInf),
  //                             _["distribution"] = "a",
  //                             _["freedom_degrees"] = 1.0,
  //                             _["ratio"] = "reference"),
  //                             "Discrete Choice Model");

  // Rcpp::function("predict_glmcat_Response", &predict_glmcat_Response,
  //                List::create(_["model_object"] = R_NaN,
  //                             _["data"] = NumericVector::create( 1, NA_REAL, R_NaN, R_PosInf, R_NegInf)
  //                ),
  //                "predict_glmcat_Response Choice Model");

  // Rcpp::class_<ReferenceF>("ReferenceF")
  //   .constructor()
  //   .method( "inverse_logistic", &ReferenceF::inverse_logistic )
  // ;
// }
