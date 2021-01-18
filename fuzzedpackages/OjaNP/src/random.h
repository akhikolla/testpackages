/* $Id: random.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef RANDOM_H
#define RANDOM_H

using namespace std; //df

void set_random_seed();

int random(int min,int max);
double Uniform(double a,double b);
double Normal(double m, double s);
double Exponential(double l);
double Gamma(double a);
double Chi2(int n);
double Fischer(int m,int n);

double gamma_2(int n);
double dfChi2(int n,double x);

double limitChi2(int n,double percentage);
#endif
