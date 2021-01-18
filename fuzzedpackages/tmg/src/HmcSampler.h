/* 
 * File:   HmcSampler.h
 * Author: aripakman
 *
 * Created on July 4, 2012, 10:44 AM
 */

#ifndef HMCSAMPLER_H
#define	HMCSAMPLER_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <random>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

struct LinearConstraint{
  VectorXd f;
  double  g;
};

struct QuadraticConstraint{
    MatrixXd A;
    VectorXd B;
    double  C;    
};


class HmcSampler   {
public:
    
    HmcSampler(const int & d, const int & seed);

    void setInitialValue(const VectorXd & initial);
    void addLinearConstraint(const VectorXd & f, const double & g);
    void addQuadraticConstraint(const MatrixXd & A, const VectorXd & B, const double & C);
    MatrixXd sampleNext(bool returnTrace = false);
    
private:
    int dim;
    VectorXd lastSample;    
    static const double min_t; 
    vector<LinearConstraint> linearConstraints;
    vector<QuadraticConstraint> quadraticConstraints;
    
    knuth_b eng1; //to sample momenta 
    uniform_real_distribution<> ud; 
    normal_distribution<> nd; 

    void _getNextLinearHitTime(const VectorXd & a, const VectorXd & b,  double & t, int & cn );
    void _getNextQuadraticHitTime(const VectorXd & a, const VectorXd & b, double & t, int & cn, const bool );
    double _verifyConstraints(const VectorXd &);
    void _updateTrace( VectorXd const & a,  VectorXd const & b, double const & tt, MatrixXd & tracePoints);
};

#endif	/* HMCSAMPLER_H */

