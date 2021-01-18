//
//  Generic.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Generic__
#define __perms_mallows__Generic__


#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include "Exponential_model.h"


#define CAYLEY_DISTANCE 0
#define KENDALL_DISTANCE 1
#define HAMMING_DISTANCE 2
#define ULAM_DISTANCE 3
#define ULAM_DISK_DISTANCE 4

#define MALLOWS_MODEL 0
#define GENERALIZED_MALLOWS_MODEL 1

#define FIRST_ITEM 1

#define BUCKET_EXP_HAM 200 //hamming sampling experiments. check every sample of this size
#define BUCKET_EXP_UL_100 40
#define BUCKET_EXP_HAM_LER 10 //hamming sampling experiments. check every sample of this size



class Generic;
using namespace std;
class Generic{
public:
    Generic(){
        facts_ = NULL;}
    
    ~Generic(){if ( facts_ != NULL) delete [] facts_;}
    
    Exponential_model* new_instance(int distance_id, int n);
    
    void elementary_symmetric_polynomial(double* theta, int n, long double*theta_exp_aux, long double **esp_aux, long double *esp);
    
    void split_elementary_symmetric_polynomial(long  double *esp, double *theta, int n, long double **esp_no_a, long double **esp_yes_a);

    void insert_at(int *sigma, int n, int move, int to, int*res);
    
   /* int**read_sample_file(int n, int m, char*path);

    void print_freq_matrix(int **samples, int m, int n);
    
    void print_long_double_matrix(long double**matrix, int n, int m);
    
    void print_double_matrix(double**matrix, int n, int m);
    
    void print_int_matrix(int**matrix, int n, int m);
    
    void print_int_vector(int*vec, int n);
    
    void print_double_vector(double*vec, int n);
        void print_float_vector(float*vec, int n);
    */    
    void generate_random_permutation(int len, int first_item_in_perm, int*sigma);
    
    void random_shuffle(int len, int*array);
    
    void compose(int n, int*s1, int*s2, int*res);
    
    void invert(int n, int*sigma, int*res);
    
    void invert_sample(int n, int m, int ** sample, int ** sample_inv);
    
    void get_permu_matrix(int n,int*sigma, int**matrix);
        
    long double factorial(int val) ;

    long double count_permus_with_at_least_k_unfixed_points(int n, int k);
    
    double count_perm_fixed_points(int k, int j);
    
    bool valid_permutation(int*sigma, int n);
    
    //void seed(void);
    
    void partition_function_init(int n);
    
    int get_number_of_partitions(int n);
    
    void partition_function_destructor(int n);
    
    void init_factorials (int n) ;
    
    void delete_factorials();
    
    double get_factorial_in_table(int i){
        if(facts_ == NULL || facts_n_ < i ){
            //run init_factorial(n). count_permus_no_fixed_points does it, 
            cout <<"load factorial error, get_factorial_in_table "<<endl;
            exit(1);
        }
        return facts_[i];
    }
    
    //void riffle_shuffle(int n, int m, int times, int *sigma, int **shuffle);
    
    void freq_matrix(int ** samples, int m, int n, int**freq);
    
    void compose_sample_right(int**samples, int*sigma, int m, int n, int ** composed);
    void compose_sample_right(int**samples, int*sigma, int m, int n);
    void ini_chrono();
    double end_chrono();
    double milliseconds();
    
private:
    
    int**partition_table;// matrix for the table of partition function (number of integer partitions of n)
    
    long double *facts_;
    
    int facts_n_;
    double chrono;
};

#endif /* defined(__perms_mallows__Generic__) */
