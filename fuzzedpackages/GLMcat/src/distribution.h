#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <RcppEigen.h>
using namespace std;
using namespace Rcpp;

class distribution{
public:
  double _epsilon_0 = 1e-10;
  double _epsilon_1 = 1e-6;

  std::string concatenate(std::string x, std::string level);

  List All_pre_data_or(Formula formula,
                       DataFrame input_data,
                       CharacterVector categories_order,
                       CharacterVector proportional_effect,
                       std::string threshold = "NA",
                       std::string ratio = "non_cum");

  List All_pre_data_NEWDATA(Formula formula,
                            DataFrame NEWDATA,
                            CharacterVector categories_order,
                            CharacterVector proportional_effect,
                            int N_cats
  );

  List select_data_nested(Formula formula,
                          String individuals,
                          String Alternatives,
                          CharacterVector ref_cat,
                          CharacterVector var_alt_specific,
                          DataFrame input_data
                            //   ,
                            // String ratio
  );

  distribution();
};

class Logistic : virtual public distribution{
public:
  virtual Eigen::VectorXd in_open_corner(const Eigen::VectorXd& p) const;
  virtual double cdf_logit(const double& value) const;
  virtual double pdf_logit(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Logistic();
};

class Normal : virtual public distribution{
public:
  virtual double cdf_normal(const double& value) const;
  virtual double pdf_normal(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Normal();
};

class Cauchit : virtual public distribution{
public:
  virtual double cdf_cauchit(const double& value) const;
  virtual double pdf_cauchit(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Cauchit();
};

class Student :  virtual public distribution{
public:
  virtual double cdf_student(const double& value, const double& freedom_degrees) const;
  virtual double pdf_student(const double& value, const double& freedom_degrees) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Student();
};

class Gumbel :  virtual public distribution{
public:
  virtual double cdf_gumbel(const double& value) const;
  virtual double pdf_gumbel(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Gumbel();
};

class Gompertz : virtual public distribution{
public:
  virtual double cdf_gompertz(const double& value) const;
  virtual double pdf_gompertz(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Gompertz();
};


#endif
