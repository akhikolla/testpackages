/* curstat.h */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <Rcpp.h>


#define SQR(x) ((x)*(x))

using namespace std;
using namespace Rcpp;


typedef struct
{
    double t;
    int freq1;
    int freq2;
}
SampleTime;

typedef struct
{
    double t;
    int delta;
}
SampleTime2;

List ComputeConfIntervals(DataFrame input, NumericVector x, double alpha, NumericVector bw);
NumericVector ComputeBW(DataFrame input, NumericVector x);
DataFrame ComputeMLE(DataFrame data);
NumericVector ComputeSMLE(DataFrame data, NumericVector x, double h);

int     CompareTime(const void *a, const void *b);
int     compare(const void *a, const void *b);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  bdf(double A, double B,  int njumps, double *jumploc, double *p, double h, double u);
double  K(double x);
double  KK(double x);
void    data_bootstrap(int N, int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[]);
void    data_bootstrap2(int N, int nB, int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[]);
double  varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u);


