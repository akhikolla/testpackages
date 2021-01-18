//
//  Cayley.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Cayley__
#define __perms_mallows__Cayley__

#include <iostream>
#include "Generic.h"
#include "Exponential_model.h"

class Cayley;
using namespace std;
class Cayley: public Exponential_model{
protected:
//    int n_;
    long double **stirling_matrix_;
    long double *facts_;
    int *sigma_inv_ ;
   
        
//    int get_most_prob_cycle(int ind, int**cycles, int len, int*leng_cycles);
    
    double calculate_psi( double*theta, double*psi_vector);
    
    void get_x_lower_bound_freq(int m, int** samples_inv_freq, int ini_pos, int*min_bound_x);
    
    double get_bound_likeli(int m, int ** samples_inv_freq, int ini_pos, int* x , int*sigma_0);
    
    void get_x_lower_bound(int m, int** sample, int ini_pos, int*x_min_bound);
    
    void generate_permu_with_k_cycles(int n, int k, int*sigma);
    
    bool same_cycle(int i, int j, int*sigma);
    
//    bool item_closes_cycle(int pos, int item, int* sigma);    int which_item_closes_cycle(int*sigma_inv, int pos);

    void get_max_item_in_future_cycles(int*sigma, int i, int j, int*max_i, int*max_j);
    void get_max_item_in_current_cycle(int*sigma, int i, int*max_i );
    
    long double count_permus_by_x_core(int i, int*x, long double*c, int c_len, int *k);
    //void local_search_dist_likeli(int m, int **samples, int*sigma, int model, int*dist2sample, double*like2sample);
    //void local_search_swap_dist_likeli(int m, int** samples, int* sigma, int model, int* dist2sample, double* like2sample);
    
    //void local_search_insert_dist_likeli(int m, int** samples, int* sigma, int model, int* dist2sample, double* like2sample);
    bool updateCycleInfoForLearning(int initial, int final, int* scan, int* xs, int* oldCycleLength,
                                    int* newCycleLength, int*maxItem, int* xAcumul, int*xAcumulVariation, bool downTree);
    int getLargest(int**freq, int*a, int*b){//y al final los pone a 0
        int max=0;
        for(int i=0;i<n_;i++)
            for(int j=0;j<n_;j++)
                if(max<freq[i][j]){
                    max=freq[i][j];
                    (*a)=i;
                    (*b)=j;
                }
        freq[(*a)][(*b)]=0;
        return max;
    }
    double get_theta_log_likelihood(int m, int*x_acumul, int*x_acumul_variation, double*theta_estim);
    
    void variable_neighborhood_search(int m, int **samples, int*sigma, int model, double*f_eval);
    void local_search_swap_mm ( int m, int ** samples, int *sigma_0, double * f_eval);
    void local_search_swap_gmm( int m, int ** samples, int *sigma_0, double * f_eval);
    void local_search_insert  ( int m, int ** samples, int *sigma_0, int model, double *f_eval);
    double estimate_consensus_exact_gmm(int m, int **samples, int*sigma_0_ini, int *sigma_0);
    double estimate_consensus_exact_mm (int m, int **samples, int*sigma_0_ini, int *sigma_0);
    double estimate_consensus_exact_gmm_core(int m, int pos, int** samples, int** samplesInv, int**samplesinvfreq,
                                             int* x_acum, int* current_sigma, int* current_sigma_inv, double current_likeli_bound,
                                             int* best_sigma, double* best_likeli);
    double estimate_consensus_exact_mm_core(int m, int pos, int** samples, int** samplesInv, int* x_acum,  int* current_sigma,
                                            int* current_sigma_inv, double current_dist_bound, int*best_sigma, double*best_dist);
    void estimate_consensus_approx_gmm(int m, int**samples, int**samples_inv, int*best_sol, double*best_likeli);
    void estimate_consensus_approx_mm (int m, int**samples, int**samples_inv, int*best_sol, double*best_distance);
    void estimate_consensus_approx_gmm_core(int m, int pos, int** samples, int** samplesInv, int**samplesinvfreq,
                                            int* x_acum, int* current_sigma, int* current_sigma_inv, double current_likeli_bound,
                                            int* best_sigma, double* best_likeli);



public:
    int     maximum_distance(){ return n_ - 1;}
    int     distance(int*s1, int*s2) ;
    double  probability(int*s, int*s_0, double*theta) ;
    void    random_sample_at_dist(int dist, int m, int **samples) ;
    void    multistage_sampling(int m, double*theta, int**samples) ;
    void    distances_sampling(int m, double theta, int**samples) ;
    void    gibbs_sampling(int m, double *theta, int model, int **samples);
    void    estimate_consensus_approx  (int model, int m, int **samples, int *sigma_0);
    double  estimate_consensus_exact( int model, int m, int **samples, int*sigma_0_ini, int *sigma_0){
        estimate_consensus_approx(model , m , samples , sigma_0);
        
        if ( model == MALLOWS_MODEL )
            return estimate_consensus_exact_mm(m, samples, sigma_0_ini, sigma_0);
        return estimate_consensus_exact_gmm(m, samples, sigma_0_ini, sigma_0);
    }
    long double get_likelihood(int m, int** samples, int model, int * sigma_0) ;
    void    estimate_theta(int m, int*sigma_0, int**samples, int model, double*theta);
    long double num_permus_at_distance(int d);
    int     perm2dist_decomp_vector(int*sigma, int*vec ) ;
    void    dist_decomp_vector2perm(int* vec, int* sigma) ;


    
    virtual double  expectation(double theta);
    virtual void    expectation(double *theta, double *expect) ;


    /******     virtual     *******/

    Cayley(int ns){
        n_=ns;
        sigma_inv_ = new int [ n_ ];
        //Generic gen;
        //gen.seed();
        stirling_matrix_ = new long double *[n_+1];
        facts_ = new long double[ n_ + 1 ];
        //for (int i = 0 ; i < n_+1; i ++ )
        for (int i = 0 ; i < n_+1; i ++ ){
            stirling_matrix_[ i ]  = new long double [n_+1];
            for(int j= 0;j<n_+1;j++)stirling_matrix_[ i ][ j ] =-1;
            if(i != 0) facts_[ i ] =facts_[i-1]*i;
            else facts_[0] = 1;
        }
        stirling_matrix_[ 0 ][ 0 ] = 1;
        for (int i = 0 ; i <= n_; i ++){
            stirling_matrix_[ i ][ i ] = 1;
            stirling_matrix_[ i ][ 0 ] = 0;
            if ( i > 0 ) stirling_matrix_[ i ][ 1 ] = facts_[ i - 1 ];
        }
        for (int i = 2 ; i <= n_; i ++)
            for (int j = 2 ; j < i ; j++)
                stirling_matrix_[ i ][ j ] = stirling_matrix_[ i - 1 ][ j - 1 ] + (i - 1) * stirling_matrix_[ i - 1 ][ j ];
    }

    ~Cayley(){
        for (int i = 0 ; i < n_+1; i ++ ) delete [] stirling_matrix_[ i ];
        delete [] stirling_matrix_;
        delete [] facts_;
        delete [] sigma_inv_;
    };
    
    int get_cycles(int*sigma, int*cycle_items, int*cycle_indices);
    
    void random_sample_at_dist_d(int d, int m, int**samples);
    
    void x_vector_to_permutation_backwards(int*x, int*sigma);
    
    void x_vector_to_permutation_forward(int*x, int*sigma);
    
    long double count_permus_with_cycles(int d);
    
    void print_stirling_matrix();
        
    int distance_to_sample(int** samples, int m, int* sigma);
        
    long double count_permus_by_x(int*x);

   


};
#endif /* defined(__perms_mallows__Cayley__) */



























