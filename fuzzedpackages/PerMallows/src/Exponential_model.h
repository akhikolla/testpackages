//
//  Exponential_model.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 27/05/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//

#ifndef perms_mallows_Exponential_model_h
#define perms_mallows_Exponential_model_h

class Exponential_model {
protected:
    int n_;
public:
    virtual ~Exponential_model(){};
    
    virtual int     maximum_distance() = 0 ;
    virtual int     distance(int*s1, int*s2) = 0;
    virtual double  probability(int*s, int*s_0, double*theta) = 0;
    virtual int     perm2dist_decomp_vector(int*sigma, int*vec ) = 0;
    virtual void    dist_decomp_vector2perm(int* vec, int* sigma) = 0;

    virtual void    random_sample_at_dist(int dist, int m, int **samples) = 0;
    virtual void    multistage_sampling(int m, double*theta, int**samples) = 0;
    virtual void    distances_sampling(int m, double theta, int**samples) = 0;
    virtual void    gibbs_sampling(int m, double *theta, int model, int **samples) = 0;
    
    virtual void    estimate_consensus_approx (int model, int m, int **samples, int *sigma_0) = 0;
    virtual double  estimate_consensus_exact  (int model, int m, int **samples, int*sigma_0_ini, int *sigma_0) = 0;
    
    virtual long double get_likelihood(int m, int** samples, int model, int * sigma_0) = 0;
    virtual void    estimate_theta(int m, int*sigma_0, int**samples, int model, double*theta) = 0 ;
    virtual long double num_permus_at_distance(int d) = 0 ;

    virtual double  expectation(double theta) = 0 ;
    virtual void    expectation(double *theta, double *expect) = 0 ;

    
};

#endif
