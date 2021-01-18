//
//  Newton_raphson.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 24/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Newton_raphson__
#define __perms_mallows__Newton_raphson__
#include <iostream>

class Newton_raphson;
using namespace std;
class Newton_raphson{
public:
    double  Newton_raphson_method(double dAvg_val, double initialGuess,
                                 int distanceModel_val, int model, int j_index, long double*count);
    void    mle_theta_weighted_mallows_hamming(int m, double*h_avg , double * theta);

    
    Newton_raphson(int n){
        n_ = n;
        esp_ = NULL;
        facts_  =   new long double[ n_ + 1 ];
        facts_[ 0 ] = 1;
        for (int k = 1 ; k <= n_ ; k ++)
            facts_[ k ] = k * facts_[ k - 1 ];

}
    ~Newton_raphson(){
        delete [] facts_;
        if( esp_ != NULL){
            for (int k = 0 ; k < n_ + 1 ; k ++){
                delete [] esp_no_a_[k];
                delete [] esp_yes_a_[k];
                delete [] aux_esp_[ k ];
            }
            delete [] esp_no_a_;
            delete [] esp_yes_a_;
            delete [] esp_;
            delete [] aux_esp_ ;
            delete [] t_;
        }
    }
        double   likeli_wmh(double x[]);
    void init_optim_wmh(){
        if ( esp_ == NULL ){
            esp_    =   new long double[ n_ + 1 ];
            esp_no_a_ = new long double*[ n_ + 1 ];
            esp_yes_a_= new long double*[ n_ + 1 ];
            aux_esp_  = new long double *[ n_ + 1 ];//esp
            t_        = new long double  [ n_ + 1 ];//esp
            for (int k = 0 ; k <= n_ ; k ++){
                esp_no_a_ [ k ]= new long double [ n_ ];
                esp_yes_a_[ k ]= new long double [ n_ ];
                aux_esp_[ k ] = new long double[ n_ + 1 ];
                for (int i = 0 ; i < n_ ; i ++){
                    esp_no_a_ [ k ][ i ] = 0;
                    esp_yes_a_[ k ][ i ] = 0;
                    aux_esp_[ k ][ i ] = 0;
                }
            }
        }
    }
    void    dlikeli_wmh(double x[], double deriv []);
    void    dlikeli_wmh_reverse(double x[], double deriv []);
    
private:
    
    //newton raphson multivariate - weighted hamming
    int     m_, n_;
    double  * h_avg_;
    int     model_;
    long double  *esp_;
    long double  **esp_no_a_;
    long double  **esp_yes_a_;
    long double * facts_;
    long double  **aux_esp_ ;//esp
    long double  *t_        ;//esp
    


    long double * count_;
    double  UPPER_THETA;
    int     distance_id_; //0:cayley; 1: kendall mm; 2: kendall GMM
    double  j_index_ ; //1...n
    double  dist_avg_ ;
    
    
    double  rtsafe  (double x1, double x2, double xacc);
    void    funcd   (double theta,  double *ff, double *ffdev);
    void    frprmn  (double p[], int n, double ftol, int *iter, double *fret, double (Newton_raphson::*func)(double []), void (Newton_raphson::*dfunc)(double [], double []));
    double  f1dim   (double x);
    void    mnbrak  (double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (Newton_raphson::*func)(double));
    double  df1dim  (double x);
    double  dbrent  (double ax, double bx, double cx, double (Newton_raphson::*f)(double), double (Newton_raphson::*df)(double), double tol, double *xmin);
    void    dlinmin (double p[], double xi[], int n, double *fret,
                     double (Newton_raphson::*func)(double []),
                     void (Newton_raphson::*dfunc)(double [], double []));
    void    linmin  (double p [] ,double xi[],int n,double *fret,double (Newton_raphson::*func) (double []));
    double  brent   (double ax,double bx,double cx,double (Newton_raphson::*f)(double),double tol,double *xmin);
    
    double f    (double theta);//newton univariate
    double fdev (double theta);

    void    nrerror(char error_text[]);
    double *vector(long nl, long nh);
    void    free_vector(double *v, long nl, long nh);

    
    
    
};

#endif 

