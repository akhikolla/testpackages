//
//
//  Created by Qingyan Xiang on 4/19/17.
//  Copyright © 2017 Qingyan Xiang. All rights reserved.
//


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define P1 10e-4
#define P2 10e-4

using namespace std;
using namespace Eigen;
using namespace Rcpp;


double mydnorm (double x, double mu, double sigmasq);


double ran1(long *idum);

/*N(0,1) generator */
double gasdev(long *idum);


/*Gamma generator  */

double gamdev(double a, double b, long *idum);// a:alpha,shape parameter; b: beta, rate parameter.

/* beta generator */
double betadev( double a, double b, long *idum );

/*truncated norm, N(0,1), truncated at a, larger than a  */
double tndev(double a, long *idum);


//inverse gaussian generator
double igasdev(double u, double l, long *idum);
