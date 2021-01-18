//
//  Ulam.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 08/07/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Ulam__
#define __perms_mallows__Ulam__

//return values of generation functions
#define GEN_NEXT  0 //ok, print and continue
#define GEN_TERM  1 //ok, terminate
#define GEN_EMPTY 2 //ok, print EMPTY SET and continue
#define GEN_ERROR 3 //an error occured, print an error message and terminate
//http://www.aconnect.de/friends/editions/computer/combinatoricode_e.html#Partitions

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include "Newton_raphson.h"
#include "Ferrers_diagram.h"
#include "Exponential_model.h"

class Ulam;
using namespace std;
class Ulam: public Exponential_model{
protected:
    long double *num_permus_per_dist_;
    long double *first_index_at_dist_;
    std::vector<Ferrers_diagram*>* shapes_of_n_;
    std::vector<long double>num_permus_at_shape_acumul_;
    long double num_partitions_of_n_;
    long double * facts_ ;
    int *comp_, *inv_ ;//auxiliary for LIS
    //fast lis
    int * M;
    int * P ;

    
    //http://www.aconnect.de/friends/editions/computer/combinatoricode_e.html#Partitions
    int gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k);
    int gen_part_next(unsigned char *vector, unsigned char *k, int bound);
    
    void    calculate_probas_at_each_distance(double theta, double*proba);
    void    fill_shapes_of_n();
    void    generate_permu_with_given_LIS(int lis, int*sigma);
    int     longest_increasing_subsequence(int*sigma);
    double  psi(double theta);
    int     get_lower_bound_sum_of_LIS(int**samples, int m);
    int     set_median(int m, int **samples, int *sigma_0);
    
    
public:
    //    int n_;
    
    int     distance_to_sample(int **samples, int m, int *sigma);
    
    
    /***    ini virtual     ****/
    
    int     maximum_distance(){ return n_ - 1;}
    int     distance(int*s1, int*s2) ;
    double  probability(int*s, int*s_0, double*theta) ;
    void    random_sample_at_dist(int dist, int m, int **samples) ;
    void    distances_sampling(int m, double theta, int**samples) ;
    void    estimate_consensus_approx (int model, int m, int **samples, int *sigma_0){
        set_median(m, samples , sigma_0);
    }
    long double get_likelihood(int m, int** samples, int model, int * sigma_0) ;
    void    estimate_theta(int m, int*sigma_0, int**samples, int model, double*theta);
    long double num_permus_at_distance(int d);

    int     perm2dist_decomp_vector(int*sigma, int*vec ) {
        //checked from R cout<<"not supported"<<endl;exit(0);
        return -1;
    }
    void    dist_decomp_vector2perm(int* vec, int* sigma) {
        //checked from R cout<<"not supported"<<endl;exit(0);
    }
    void    multistage_sampling(int m, double*theta, int**samples) {
        //checked from R cout<<"not supported"<<endl;exit(0);
    }
    void    gibbs_sampling(int m, double *theta, int model, int **samples) ;
    double  estimate_consensus_exact  (int model, int m, int **samples, int*sigma_0_ini, int *sigma_0){
        //checked from R cout<<"not supported"<<endl;exit(0);
        return -1;
    }
    
    /***    end virtual     ****/
    //fast lis
    int search_LIS(int* M, int* A, int i, int L ) ;

    Ulam(int n){
        n_ = n;
        num_partitions_of_n_ = -1;
        shapes_of_n_ = new std::vector<Ferrers_diagram*>();
        first_index_at_dist_ = new long double [ n_ ];
        num_permus_per_dist_ = new long double [ n_ ];
        for (int i = 0 ; i < n_ ; i++) num_permus_per_dist_[ i ] = 0 ;
        facts_     = new long double [ n + 1 ];
        facts_[0]=1;
        for (int i = 1 ; i < n + 1 ;i++)
            facts_[i] = facts_[i-1] * i;
        comp_   = new int[n_];
        inv_    = new int[ n_];
        M       = new int[n_];
        P       = new int[n_];
    }
    ~Ulam(){
        for ( int i = 0 ; i < shapes_of_n_->size() ; i++ ) delete shapes_of_n_->at(i);
        shapes_of_n_->clear();
        delete shapes_of_n_ ;
        delete [] first_index_at_dist_;
        delete [] num_permus_per_dist_;
        delete [] facts_;
        delete [] comp_;
        delete [] inv_;
        delete [] M;
        delete [] P;
    }

    int     distance(int*sigma){  return n_ - longest_increasing_subsequence(sigma);    }

    long double num_permus_at_distance_approx(int d);

    double  expectation(double theta);
    void    expectation(double *theta, double *expect){
        //checked from R cout<<"not supported"<<endl;exit(0);
    }

};
#endif /* defined(__perms_mallows__Ulam__) */
