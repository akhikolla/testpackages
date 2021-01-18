//
//  Hamming.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Hamming__
#define __perms_mallows__Hamming__

#include <cmath>
#include <iostream>
#include <algorithm>
#include "Generic.h"
#include "Exponential_model.h"

class Hamming;
using namespace std;
class Hamming: public Exponential_model{
public:
    int     maximum_distance(){ return n_ ;}
    int     distance(int*s1, int*s2) ;
    double  probability(int*s, int*s_0, double*theta) ;
    void    random_sample_at_dist(int dist, int m, int **samples) ;
    void    multistage_sampling(int m, double*theta, int**samples) ;
    void    distances_sampling(int m, double theta, int**samples) ;
    void    gibbs_sampling(int m, double *theta, int model, int **samples) ;
    void    estimate_consensus_approx (int model, int m, int **samples, int *sigma_0){
        if (model == MALLOWS_MODEL ){
            //cout<<"No approx learning for MM under the Hamming distance. Try the exact learning. "<<endl;
            //exit(0);
        }else{
            double likeli;
            estimate_consensus_approx_gmm( m, samples, sigma_0, &likeli);
        }
    }
    double  estimate_consensus_exact  (int model, int m, int **samples, int*sigma_0_ini, int *sigma_0){
        if (model == MALLOWS_MODEL ) estimate_consensus_exact_mm( m, samples, sigma_0);
        //else {cout<<"No exact learning for WMM under the Hamming distance. Try the approx learning. "<<endl;}
        return 0;
    }
    long double get_likelihood(int m, int** samples, int model, int * sigma_0) ;
    void estimate_theta(int m, int*sigma_0, int**samples, int model, double*theta);
    long double num_permus_at_distance(int d){
        return count_derangements(d) * facts_[n_] / (facts_[d] * facts_[ n_ - d ]) ;
    }
    int     perm2dist_decomp_vector(int*sigma, int*vec ) ;
    void    dist_decomp_vector2perm(int* vec, int* sigma) ;

    /***********    end virtual    **************/
    void learning_experiments_approx_2(int m, int **freq_neg, int *time_lap, int*time_vns, int *sigma_0, int * dist_lap, int * dist_vns, double * likeli_lap, double * likeli_vns);
    void learning_experiments_approx_1(int m, double *theta, int ** freq_neg );

    /******     experiments ***/

    Hamming(int n){n_ = n; 
        Generic gen;
        esp_red_   = new long  double[ n + 1 ];
        esp_ini_   = new long  double[ n + 1 ];        
        esp_red_yes_a_   = new long  double[ n + 1 ];
        facts_     = new long double [ n + 1 ];
        facts_[ 0 ] = 1 ;
        for (int i = 1 ; i <= n ; i++)
            facts_[i] = facts_[ i - 1 ] * i;
        g_n_      = new long double*[ n + 1 ];
        aux_esp_  = new long double *[ n + 1 ];//esp
        t_        = new long double  [ n + 1 ];//esp
        t_sampling_= new long double [ n_ ];
        deran_num_ = new double [ n_ + 1 ];
        deran_num_[ 0 ] = 1;
        deran_num_[ 1 ] = 0;
        for (int i = 2 ; i <= n_ ; i ++)
            deran_num_[i] = deran_num_[i-1]*(i-1) + deran_num_[i-2]*(i-1);
        for ( int i = 0 ; i < n ; i ++ )
            t_sampling_[ i ] = 0;//(long double)exp( theta[ i ]) - 1 ;

        gen.init_factorials(n);
        for ( int i = 0 ; i < n + 1 ; i ++ ){
            g_n_ [ i ]      = new long double [ n_ + 1 ];
            aux_esp_[ i ]   = new long double [ n_ + 1 ];
            for ( int j = 0 ; j <= i ; j ++)
                g_n_[ i ][ j ] = gen.count_permus_with_at_least_k_unfixed_points(i, j);
            for ( int j = 0 ; j <= n ; j ++)
                aux_esp_[i][j]=0;
        }
        
    }
    
    ~Hamming(){
        delete [] deran_num_;
        delete [] esp_red_   ;
        delete [] esp_ini_;
        delete [] esp_red_yes_a_   ;
        delete [] facts_     ;
        delete [] t_;
        delete [] t_sampling_;
        for (int i = 0 ; i < n_ + 1 ; i++){
            delete [] aux_esp_[ i ];
            delete [] g_n_[ i ];
        }
        delete [] aux_esp_;
        delete [] g_n_       ;
    }
        
    
    
    void multistage_sampling_experiments(int m, double*theta, double*error, double*time);
    void distances_sampling_experiments(int m, double theta, double *error, double*time) ;
    void gibbs_sampling_experiments(int m, double*theta, double *error, double*time);
    
    void sample(double* theta, int m, int**samples);

    int distance_to_sample(int **samples, int m, int *sigma);
    

    void estimate_consensus_exact_mm(int m, int**samples, int* sigma_0);//linear assignment
//    double estimate_consensus_exact_gmm(int m, int**samples, int* sigma_0_ini, int* sigma_0);
//    double estimate_consensus_exact_gmm_core (int m, int**freq,  int * max_index_in_col, double * h_avg, int* sigma_0, int * sigma_0_inv, int pos, double current_likelihood, double * best_likelihood, int* best_sigma_0);
    
    void estimate_consensus_approx_gmm (int m, int**samples, int* sigma_0, double * best_likelihood);
    
    void random_derangement(int n, int*sigma);
    
    double count_derangements(int n){
        return deran_num_[n];
    }
    
    
    void generate_h (double*theta, int*h);
//    long double compute_marginal_slow ( int *h, double *theta, int marginal_order);
    long double compute_marginal_iterative ( int *h, double *theta, int marginal_order);
    long double compute_marginal(int *h , double * theta );
    
    
    void    sample_to_h_vector(int **samples, int m, int * sigma, double *h_avg);


    double  expectation(double theta);
    void    expectation(double *theta, double*h_expected );
    long double get_likelihood_from_h(int m, int model, double *theta, double * h_avg);


protected:
    
    //int n_;
    long double *esp_red_       ;//compute marginal
    long double *esp_ini_   ;
    long double *esp_red_yes_a_;//compute marginal
    long double theta_acum_not_in_A; //compute marginal
    int b_ ;                    //compute marginal
    long double * t_sampling_;  //compute marginal: 0..n-1
    long double *facts_;
    long double **g_n_;         //compute marginal
    long double  **aux_esp_ ;//esp
    long double  *t_        ;//esp
    
    
    double*deran_num_; // coutn the number of derangements of n items
    
    double  psi_whm_reverse(double*theta);
    double  psi_hm(double theta);
    double  psi_hm_reverse(double theta);
    double  psi_whm(double*theta);

    void    random_permu_at_dist_d(int dist, int*sigma );
    
    void    generate_permu_from_list(int*ran, int dist, int*sigma);
    
/*    void    variable_neighborhood_search(int m, int **freq, int*sigma_0, double*f_eval);
    void    local_search_swap_gmm  ( int m, int ** samples, int *sigma_0, double * f_eval);
    void    local_search_insert_gmm( int m, int ** samples, int *sigma_0, double * f_eval);
    void    bound_consensus (int m, int pos, int ** freq, int*sigma_0, int*sigma_0_inv, double*h_avg,int* max_index_in_col, double*h_avg_bounded,int* max_index_in_col_bounded);
*/

    void quasyIndepTest(double* theta);
    int indexOfPermu(int* sigma){//chapu
        int index=0;int*set=new int[n_];for(int i=0;i<n_;i++)set[i]=i+1;
        Generic gen;
        for(int i=0;i<n_;i++)
            index += gen.factorial(n_-1-i) * indexOf(sigma[i],set);
        delete[]set;
        return index;
    }
    int indexOf(int j, int* set){
        int cont=0;
        for(int i=0;i<n_;i++){
            if (set[i]==j) {set[i]=-1;return cont;}
            if (set[i]!=-1) cont++;
        }
        return -1;
    }

};
#endif /* defined(__perms_mallows__Hamming__) */
