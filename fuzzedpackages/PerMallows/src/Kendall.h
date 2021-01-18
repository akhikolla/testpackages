//
//  Kendall.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Kendall__
#define __perms_mallows__Kendall__

#include <iostream>
#include "Exponential_model.h"
#include "Generic.h"
#include <cstdlib>
class Kendall;
using namespace std;
class Kendall: public Exponential_model{
protected:
    //    int n_;
    long double**count_;
    
public:
    int     maximum_distance(){ return n_*(n_-1)/2;}
    int     distance(int*s1, int*s2) ;
    double  probability(int*s, int*s_0, double*theta) ;
    void    random_sample_at_dist(int dist, int m, int **samples) ;
    void    multistage_sampling(int m, double*theta, int**samples) ;
    void    distances_sampling(int m, double theta, int**samples) ;
    void    gibbs_sampling(int m, double *theta, int model, int **samples) ;
    void    estimate_consensus_approx  (int model, int m, int **samples, int *sigma_0){
        borda(samples, m, sigma_0);
    }
    double  estimate_consensus_exact  (int model, int m, int **samples, int*sigma_0_ini, int *sigma_0){
        //cout<<"Not implemented. See Mandhani, B., & Meila, M. (2009). Tractable Search for Learning Exponential Models of Rankings. Journal of Machine Learning Research, 5, 392–399."<<endl;
        //exit(0);
        return -1;
    }
    void    estimate_theta(int m, int*sigma_0, int**samples, int model, double*theta);
    int     perm2dist_decomp_vector(int*sigma, int*v ) ;
    void    dist_decomp_vector2perm(int* v, int* permu) ;
    long double get_likelihood(int m, int** samples, int model, int * sigma_0) ;
    long double num_permus_at_distance(int d);
    
    /***    end virtual     ****/
    
    Kendall(int n){
        n_=n;
        count_=new long double*[n_+1];
        for(int i=0;i<n_+1;i++) count_[i]=new long double[n_*(n_-1)/2+1];
        //https://oeis.org/A008302
        for(int i=0;i<n_+1;i++)
            for(int j=1;j<n_*(n_-1)/2+1;j++)
                count_[i][j]=0;
        for(int i=0;i<n_+1;i++) count_[i][0]=1;
        for(int i=1;i<n_+1;i++){
            for(int j=1;j<i*(i-1)/2+1;j++){
                if(j-i >= 0) count_[i][j] = count_[i][j-1] + count_[i-1][j] - count_[i-1][j-i];
                else count_[i][j] = count_[i][j-1] + count_[i-1][j] ;
            }
        }
    }
    ~Kendall(){
        for(int i=0;i<n_+1;i++) delete[]count_[i];
        delete[]count_;
    };
    

    double expectation(double theta);
    void expectation (double * theta, double*expect);
    // int get_distance_and_v(int*sigma, int*v );
    
    //void v_vector_to_permutation(int* v, int* permu);
    
    
    

    
protected:
    
    void borda(int **samples, int m, int*sigma_0);
    
    double calculate_psi( double*theta, double*psi);
    
    int distance_to_sample(int **samples, int m,  int *sigma);
    
    void random_permu_at_dist_d(int dist, int*sigma );
    
    //    long double num_permus_with_k_inversions(int k);
    
};
#endif /* defined(__perms_mallows__Kendall__) */
