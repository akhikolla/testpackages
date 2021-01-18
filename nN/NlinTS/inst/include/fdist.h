/**
 * @authors Hmamouche Youssef
 **/

#ifndef FDIST_H
#define FDIST_H


double gammaln(double xx) ;
double beta(double x, double y);

double fi(int N, double x, double a, double b);

double inBeta(double x, double a, double b);
double getPvalue(double f, double n1, double n2) ;
double getStudent (double t, double n);

#endif // FDIST

