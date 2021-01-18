
//
//  main.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include <iostream>
#include "Cayley.h"
#include "Ulam.h"
#include "Ulam_disk.h"
#include "Hamming.h"
#include "Kendall.h"
#include "Generic.h"
#include "Exponential_model.h"


int get_int_with_msg(char*str);

void user_menu();
void loop_menu();
void op_1();
void op_2();
void op_3();
void op_4();
void op_5();
void op_6();
void op_7();
void op_8();
void op_9();
void op_10();
void op_11();
void op_12();
void op_13();
void op_14();
void op_15();
void op_16();
void op_17();
void op_18();
void op_19();

int get_int_with_msg(char*str);
int get_theta(double*theta,int  n, int dist_id);

void test_uar_end(int dist);
    void test_uar(int*sigma);
int getIndex(int**matrix, int*sigma);
bool is_equal_to(int*s1, int*s2);
void test_uar ( int**sample, int m);
void init_test_uar(int len, int n_);
void run_test();


/**********************           *************************/
int dist_id_, n_;
Exponential_model * exp_mod_;
/**********************           *************************/



/*****************************       PRJ       ************************************
int main(int argc, const char * argv[])
{ 
    Generic gen;
    gen.seed();
    loop_menu();
    return 0;
}
**********************     GENERAL      *************************

int get_int_with_msg(char*str){
    cout<<str<<endl<<">";
    int k;
    cin>>k;
    return k;
}

int get_theta(double*theta,int  n, int dist_id){
    bool    cond = false;
    int     p, num_params_gmm = n-1;
    if (dist_id == HAMMING_DISTANCE) num_params_gmm = n;
    do{
        cout<<"How many params?? (1 for the MM, "<<num_params_gmm<<" for the GMM): "<<endl<<">";
        cin>>p;
        cond = false;
        if ( p == 1 ) cond = true;
        if ( p == n-1 && ( dist_id == KENDALL_DISTANCE || dist_id == CAYLEY_DISTANCE) ) cond = true;
        if ( p == n && dist_id == HAMMING_DISTANCE ) cond = true;
    }while ( ! cond );
    for(int i=0;i<p;i++){
        cout<<"Intro theta["<<i<<"]: "<<endl<<">";
        cin>>theta[i];
    }
    if(p==1) {
        for(int i=1;i<n;i++)theta[i]=theta[0];
        return MALLOWS_MODEL;
    }
    return GENERALIZED_MALLOWS_MODEL;
}

**********************           *************************

void user_menu(){
    cout<<"\t 1.\tGenerate a random sample at distance d"<<endl;
    cout<<"\t 2.\tGenerate a random sample with the distances sampling method"<<endl;
    cout<<"\t 3.\tGenerate a random sample with the multistage sampling method"<<endl;
    cout<<"\t 4.\tGenerate a random sample with the gibbs sampling method"<<endl;
    cout<<"\t 5.\tLean approx"<<endl;
    cout<<"\t 6.\tLean exact"<<endl;
    cout<<"\t 7.\tCalculate the probability of a permutation"<<endl;
    cout<<"\t 8.\tNumber of permus with at dist d"<<endl;
    cout<<"\t 9.\tDistance decomposition"<<endl;
    cout<<"\t10.\tDistance decomposition"<<endl;
    cout<<"\t11.\tCalculate distance"<<endl;
    cout<<"\t15.\tExpected distance or distance decomposition"<<endl;
    cout<<"\t \t ------ Distance specific ------"<<endl;
    cout<<"\t12.\tCompute the Integer partitions (ULAM)"<<endl;
    cout<<"\t13.\tCycle decoposition (CAYLEY)"<<endl;
    cout<<"\t14.\tCompute marginal (HAMMING)"<<endl;
    cout<<"\t16.\tDistances sampling with data in disk - file generation step (ULAM_DISK)"<<endl;
    cout<<"\t17.\tDistances sampling with data in disk - sampling step (ULAM_DISK)"<<endl;
    cout<<"\t18.\tEstimate the paramters of a sample with data in disk - learning step (ULAM_DISK)"<<endl;
    cout<<"\t19.\tNumber of permus with at dist d with data in disk - learning step (ULAM_DISK)"<<endl;

    cout<<"\n\t 0.\t Quit. "<<endl;
}

void loop_menu(){
    Generic gen;
    char str_msg[100];
    strcpy(str_msg, "Intro n (number of items of the permutations): ");
    dist_id_ = CAYLEY_DISTANCE;
    n_       = get_int_with_msg(str_msg);
    exp_mod_ = gen.new_instance(dist_id_, n_);
    
    int op;
    do{
        user_menu();
        strcpy(str_msg, "Intro an option: ");
        op=get_int_with_msg(str_msg);
        switch (op) {
            case 1:
                op_1();
                break;
            case 2:
                op_2 ();
                break;
            case 3:
                op_3 ();
                break;
            case 4:
                op_4 ();
                break;
            case 5:
                op_5();
                break;
            case 6:
                op_6();
                break;
            case 7:
                op_7();
                break;
            case 8:
                op_8();
                break;
            case 9:
                op_9();
                break;
            case 10:
                op_10();
                break;
            case 11:
                op_11();
                break;
            case 15:
                op_15();
                break;
            case 12:
                op_12();
                break;
            case 13:
                op_13();
                break;
            case 14:
                op_14();
                break;
            case 16:
                op_16();
                break;
            case 17:
                op_17();
                break;
            case 18:
                op_18();
                break;
            case 19:
                op_19();
                break;
            default:
                break;
        }
    }while (op!=0);
    delete exp_mod_;
}



**********************     KENDALL      *************************



void op_1(){
    int d, m;
    char str_msg[100];
    strcpy(str_msg, "Intro a distance: ");
    
    Generic gen;
    do{
        d=get_int_with_msg(str_msg);
    }while(d>= exp_mod_->maximum_distance() || d <= 0);
    strcpy(str_msg, "Intro a num of permus in the sample: ");
    do{
        m=get_int_with_msg(str_msg);
    }while(m <= 0);
    int** samples=new int*[m];
    exp_mod_->random_sample_at_dist(d, m, samples);
    gen.print_int_matrix(samples, m, n_);
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
}

void op_2(){
    char str_msg[100];
    double theta;
    Generic gen;
    strcpy(str_msg, "Intro the number of permus in the sample: ");
    int m = get_int_with_msg(str_msg);
    int**samples=new int*[m];
    cout<<"Intro theta: "<<endl;
    cin>>theta;
    exp_mod_->distances_sampling(m, theta, samples);
    cout<<"Collection of permutations in the sample: "<<endl;
    gen.print_int_matrix(samples, m, n_);
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
}

void op_3(){
    char str_msg[100];
    double*theta=new double[n_];
    Generic gen;
    get_theta(theta, n_, dist_id_);
    strcpy(str_msg, "Intro the number of permus in the sample: ");
    int m = get_int_with_msg(str_msg);
    int**samples=new int*[m];
    exp_mod_->multistage_sampling(m, theta, samples);
    cout<<"Collection of permutations in the sample: "<<endl;
    gen.print_int_matrix(samples, m, n_);
    delete[]theta;
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
}


void op_4(){
    char str_msg[100];
    Generic gen;
    double*theta=new double[n_];
    int model=get_theta(theta, n_, dist_id_);
    strcpy(str_msg, "Intro the number of permus in the sample: ");
    int m = get_int_with_msg(str_msg);
    int**samples=new int*[m];
    exp_mod_->gibbs_sampling(m, theta, model, samples);
    cout<<"Collection of permutations in the sample: "<<endl;
    gen.print_int_matrix(samples, m, n_);
    delete[]theta;
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
}

void op_5(){
    char str_msg[100];
    double*theta=new double[n_], *theta_estim = new double[n_];;
    Generic gen;
    int*sigma_0=new int[n_];
    int model = get_theta(theta, n_, dist_id_);
    if ( dist_id_ == HAMMING_DISTANCE && model == MALLOWS_MODEL){
        cout<<"No approx learning for the MM under the Hamming distance. "<<endl;
        return;
    }
    strcpy(str_msg, "Intro the number of permus in the sample: ");
    int m = get_int_with_msg(str_msg);
    int**samples=new int*[m];
    exp_mod_->multistage_sampling(m, theta, samples);
    cout<<"Collection of permutations in the sample: "<<endl;
    gen.print_int_matrix(samples, m, n_);
    exp_mod_->estimate_consensus_approx(model, m, samples, sigma_0);
    exp_mod_->estimate_theta(m, sigma_0, samples , model , theta_estim);
    cout<<"Estimated consensus: "; gen.print_int_vector(sigma_0, n_);
    cout<<"Estimated theta: ";
    if ( model == GENERALIZED_MALLOWS_MODEL ) gen.print_double_vector(theta_estim, n_);
    else cout<<theta_estim[ 0 ]<<endl;
    
    delete[]theta;
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
    delete [] sigma_0;
    delete [] theta_estim;
}

void op_6(){
    //exact learning
    if ( dist_id_ != CAYLEY_DISTANCE){
        cout<<"No exact learning for this distance. "<<endl;
        return;
    }
    char str_msg[100];
    double*theta=new double[n_], *theta_estim = new double[n_];;
    Generic gen;
    int     *sigma_0 = new int[n_];
    int     *sigma_0_estim = new int[n_];
    int     model    = get_theta(theta, n_, dist_id_);
    strcpy(str_msg, "Intro the number of permus in the sample: ");
    int m = get_int_with_msg(str_msg);
    int     **samples = new int*[m];
    
    exp_mod_->multistage_sampling(m, theta, samples);
//    cout<<"Collection of permutations in the sample: "<<endl;gen.print_int_matrix(samples, m, n_);
    gen.generate_random_permutation(n_ , FIRST_ITEM , sigma_0);
    gen.compose_sample_right(samples , sigma_0, m , n_);
    exp_mod_->estimate_consensus_exact (model, m, samples, sigma_0, sigma_0_estim);
    
    exp_mod_->estimate_theta(m, sigma_0, samples , model , theta_estim);
    cout<<"Original  consensus: "; for(int i = 0 ; i < n_ ; i++) cout<<sigma_0[ i ]<<" ";
    cout<<exp_mod_->get_likelihood(m , samples , model , sigma_0)<<endl ;
    cout<<"Estimated consensus: ";  for(int i = 0 ; i < n_ ; i++) cout<<sigma_0_estim[ i ]<<" ";
    cout<<exp_mod_->get_likelihood(m , samples , model , sigma_0_estim)<<endl;
    cout<<"Estimated theta: ";
    if ( model == GENERALIZED_MALLOWS_MODEL ) gen.print_double_vector(theta_estim, n_);
    else cout<<theta_estim[ 0 ]<<endl;
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
    delete [] sigma_0;
    delete [] sigma_0_estim;
    delete [] theta_estim;
    delete [] theta;
}


void op_7(){
    int*sigma=new int[n_], *sigma_0=new int[n_];
    double *theta=new double[n_];
    for(int i=0;i<n_;i++) sigma_0[i] = i+1;
    Generic gen;
    gen.generate_random_permutation(n_, 1, sigma);
    int model;
    model = get_theta(theta, n_, dist_id_);
    cout<<"The consensus permutation is the ident. and the probability of ";
    gen.print_int_vector(sigma, n_);
    cout<<" is "<<exp_mod_->probability(sigma , sigma_0, theta)<<endl<<endl<<endl ;
    delete [] theta;
    delete [] sigma;
    delete [] sigma_0;
}

void op_8(){
    char str_msg[100];
    strcpy(str_msg, "Intro a distance: ");
    int dist;
    do{
        dist=get_int_with_msg(str_msg);
    }while(dist <0 || dist > exp_mod_->maximum_distance());
    cout<<"There are "<<exp_mod_->num_permus_at_distance(dist)<<" permutations of "<<n_<<" items at distance "<<dist<<endl;
}

void op_9(){
    if ( dist_id_ == ULAM_DISTANCE ){
        cout<<"No distance decomposition for Ulam nor Kendall distance. "<<endl;
        return;
    }
    int*sigma=new int[n_], *vec=new int[n_];
    Generic gen;
    gen.generate_random_permutation(n_, 1, sigma);
    cout<<"The distance decomp of permutation ";
    gen.print_int_vector(sigma, n_);
    exp_mod_->perm2dist_decomp_vector(sigma, vec );
    cout<<" is ";
        gen.print_int_vector(vec, n_);
    delete [] vec ;
    delete [] sigma;
}

void op_10(){
    if ( dist_id_ == ULAM_DISTANCE ){
        cout<<"No distance decomposition for Ulam nor Kendall distance. "<<endl;
        return;
    }
    int*sigma=new int[n_], *vec=new int[n_];
    Generic gen;
    bool valid;
    int max = n_ - 1; if (dist_id_ == HAMMING_DISTANCE ) max = n_;
    for (int i = 0 ; i < max ; i ++){
        do{
            cout<<"vec["<<i<<"]: "<<endl<<">"; cin>>vec[ i ];
            if (dist_id_ == KENDALL_DISTANCE )   valid = (vec[ i ] >= 0 && vec[ i ] < n_ - i );
            if (dist_id_ == CAYLEY_DISTANCE )    valid = (vec[ i ] == 0 || vec[ i ] == 1);
            if (dist_id_ == HAMMING_DISTANCE )   valid = (vec[ i ] == 0 || vec[ i ] == 1);
        }while (!valid);
    }
    exp_mod_->dist_decomp_vector2perm(vec, sigma);
    cout<<"A u.a.r. permutation consistent with the decomposition vector is ";
    gen.print_int_vector(sigma, n_);
    delete [] vec ;
    delete [] sigma;
}

void op_11(){
    Generic gen;
    int * s1 = new int[n_], *s2 = new int[n_];
    gen.generate_random_permutation(n_, 1, s1);
    gen.generate_random_permutation(n_, 1, s2);
    cout<<"The distance between ";
    gen.print_int_vector(s1, n_);
    cout<<" \t\t\t\t and ";
    gen.print_int_vector(s2, n_);
    cout<<" \t\t\t\t is "<< exp_mod_->distance(s1, s2)<<endl<<endl  ;
    delete [] s1;
    delete [] s2;
}
void op_15(){
    Generic gen;
    double  * h_avg = new double[ n_ ];
    double  * theta = new double [ n_ ];
    int     model;
    
    model = get_theta(theta, n_ , dist_id_);
    if (model == GENERALIZED_MALLOWS_MODEL ) {
        exp_mod_->expectation(theta, h_avg);
        cout <<" The expected distance decomposition vector equals ";gen.print_double_vector(h_avg , n_ );
    }
    else cout <<" The expected distance equals "<<exp_mod_->expectation(theta[ 0 ])<<endl ;
    cout<<endl;
    delete [] h_avg;
    delete [] theta;
}

**********************     DISTANCE SPECIFIC      *************************
void op_12(){
    if (dist_id_ != ULAM_DISTANCE ){
        cout<<"Ulam distance required"<<endl;
        return;
    }
    Ulam * ul = dynamic_cast<Ulam*>(exp_mod_);
    cout<<"The integer partitions of "<<n_<<" is as follows: ";
    ul->integer_partitions(n_);
}

void op_13(){
    //cycle decomp
    if (dist_id_ != CAYLEY_DISTANCE ){
        cout<<"Cayley distance required"<<endl;
        return;
    }

    Generic gen;
    int     * sigma = new int[n_];
    bool    *visited = new bool[ n_ ];
    int     item_index, cont=0;
    for (int i = 0 ; i < n_ ; i ++ )visited[ i ] =false;
    
    gen.generate_random_permutation(n_, 1, sigma);
    cout<<"The cycle decomposition of permutation "; gen.print_int_vector(sigma, n_);
    cout<<"\t\t is as follows: ";
    while(cont < n_) {
        item_index= 0;
        while ( visited[ item_index ] )item_index++;
        cout<<"( ";
        while (!visited[ item_index ] ) {
            visited[item_index] = true;
            cout<< item_index + 1 <<" ";
            item_index=sigma[item_index]-1;
            cont++;
        }
        cout<<") ";
    }
    cout<<endl<<endl;
    delete [] visited;
    delete [] sigma;
}

void op_14(){
    if (dist_id_ != HAMMING_DISTANCE ){
        cout<<"Hamming distance required"<<endl;
        return;
    }

    Generic gen;
    Hamming * ham   = dynamic_cast<Hamming*>(exp_mod_);
    double  * theta = new double [ n_ ];
    int     * h     = new int[ n_ ];

    get_theta(theta, n_, dist_id_);
    cout<<"Insert the partial vector decomposition, i.e., \n\t0: fixed point \n\t1: unfixed point \n\tother int: unknown \n";
    for (int i = 0 ; i < n_ ; i ++) {cout<<"H["<<i<<"]: "<<endl<<">"; cin>>h[ i ];}
    cout<<"The marginal probability of H = ";gen.print_int_vector(h, n_);
    cout<<"\t equals "<<ham->compute_marginal(h , theta)<<endl <<endl ;
    delete [] h;
    delete [] theta;
}

void op_16(){
    if (dist_id_ != ULAM_DISK_DISTANCE ){
        cout<<"Ulam_disk distance required"<<endl;
        return;
    }
    Generic gen;
    Ulam_disk * ul   = dynamic_cast<Ulam_disk*>(exp_mod_);
    ul->save_counts_to_file();
}

void op_17(){
    if (dist_id_ != ULAM_DISK_DISTANCE ){
        cout<<"Ulam_disk distance required"<<endl;
        return;
    }
    op_2();
}

void op_18(){
    if (dist_id_ != ULAM_DISK_DISTANCE ){
        cout<<"Ulam_disk distance required"<<endl;
        return;
    }
    char str_msg[100];
    double theta, theta_estim;
    Generic gen;
    int *sigma_0=new int[n_];
    strcpy(str_msg, "Intro theta: ");
    theta = get_int_with_msg(str_msg);
    strcpy(str_msg, "Intro the number of permus in the sample: ");
    int m = get_int_with_msg(str_msg);
    int **samples=new int*[m];
    exp_mod_->distances_sampling(m, theta, samples);
    cout<<"Collection of permutations in the sample: "<<endl;
    gen.print_int_matrix(samples, m, n_);
    exp_mod_->estimate_consensus_approx(MALLOWS_MODEL, m, samples, sigma_0);
    exp_mod_->estimate_theta(m, sigma_0, samples , MALLOWS_MODEL , &theta_estim);
    cout<<"Estimated consensus: "; gen.print_int_vector(sigma_0, n_);
    cout<<"Estimated theta: ";
    cout<<theta_estim<<endl;
    
    for(int i=0;i<m;i++)delete [] samples[i];
    delete [] samples;
    delete [] sigma_0;
}

void op_19(){
    if (dist_id_ != ULAM_DISK_DISTANCE ){
        cout<<"Ulam_disk distance required"<<endl;
        return;
    }
    op_8();
}


****************    AUXILIAR TO TEST U.A.R. ***************

int**perms;
int test_len;
int*c;
int test_n ;
void init_test_uar(int len, int n_){
    test_len=len+3;
    test_n = n_;
    perms = new int*[test_len];
    c=new int[test_len];
    for(int i=0;i<test_len;i++)c[i]=0;
    for(int i=0;i<test_len;i++)perms[i]=new int[n_];
    for(int i=0;i<test_len;i++)for(int j=0;j<n_;j++) perms[i][j]=-1;
}
bool is_equal_to(int*s1, int*s2){
    for(int i=0;i<test_n;i++)
        if(s1[i] != s2[i])return false;
    return true;
}
int getIndex(int**matrix, int*sigma){
    for(int i=0;i<test_len;i++){
        if(is_equal_to(matrix[i], sigma))return i;
        if(matrix[i][0]==-1){
            for(int j=0;j<test_n;j++)matrix[i][j]=sigma[j];
            return i;
        }
    }
    return -1;
}

void test_uar(int*sigma){
    int in = getIndex(perms, sigma);
    c[in]++;
}
void test_uar ( int**sample, int m){
    for (int i = 0 ; i < m ; i++){
        int in = getIndex(perms, sample[i]);
        c[in]++;
    }
}
void test_uar_end(int dist){
    Generic gen;
    cout<<"vector "<<endl;
    for(int i=0;i<test_len;i++){
        cout<<c[i]<<"  x  ";
        gen.print_int_vector(perms[i], test_n);
    }
    cout<<"If there are not 3 empty .. error!!!"<<endl;
}

*/




