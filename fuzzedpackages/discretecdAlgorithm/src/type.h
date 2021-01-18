#ifndef _MYTYPE_H_
#define _MYTYPE_H_

// #include "Eigen/Dense"
#include <RcppEigen.h>
//
//typedef enum {
//    BUGGY_KINDERMAN_RAMAGE,
//    AHRENS_DIETER,
//    BOX_MULLER,
//    USER_NORM,
//    INVERSION,
//    KINDERMAN_RAMAGE
//} N01type;

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
typedef Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic> MatrixXMXd;
typedef Eigen::Matrix<Eigen::MatrixXi, Eigen::Dynamic, Eigen::Dynamic> MatrixXMXi;
typedef Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic, 1> VectorXMXd;
typedef Eigen::Matrix<Eigen::MatrixXi, Eigen::Dynamic, 1> VectorXMXi;

typedef Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> MatrixXVXd;
typedef Eigen::Matrix<Eigen::VectorXi, Eigen::Dynamic, Eigen::Dynamic> MatrixXVXi;
typedef Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, 1> VectorXVXd;
typedef Eigen::Matrix<Eigen::VectorXi, Eigen::Dynamic, 1> VectorXVXi;

#endif
