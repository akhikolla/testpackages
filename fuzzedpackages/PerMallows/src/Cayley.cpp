//
//  Cayley.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.

#include <R.h>
#include "Cayley.h"
#include "Generic.h"
#include "Newton_raphson.h"
#include "Ulam.h"
#include <cmath>
#include <cfloat>

double  Cayley::expectation(double theta) {
        double sum=0, ex, denom;
        for(int j = 1 ; j < n_ ; j++){
            ex = exp(theta);
            denom = j+ex;
            sum += (double)j/denom;
            //sum += (double) j /(double)(exp(theta)+(double)j);
        }
        return (double)(sum );
    }
    
void    Cayley::expectation(double *theta, double *expect) {
  for (int i = 0 ; i < n_ - 1; i ++){
    int j = i + 1;
    expect[ i ] = (n_-j) / ((n_-j)+exp(theta[i]));
  }
  expect[ n_ - 1 ] = 0;
}

double Cayley::probability(int *s, int *s_0, double *theta){
    int     *x = new int[ n_ ],*comp = new int[ n_ ], *inv = new int[ n_ ];
    double  *psi    = new double[ n_ ];
    double  proba   =1;
    for (int i = 0 ; i < n_; i ++ ) inv[s_0[ i ]-1] =i+1;
    for (int i = 0 ; i < n_; i ++ ) comp[ i ] =s[inv[ i ]-1];
    perm2dist_decomp_vector(comp, x);
    calculate_psi(theta, psi);
    for(int i = 0 ; i < n_ - 1 ; i ++ ) proba *= exp(-theta[ i ]*x[ i ])/psi[ i ];
    delete []psi;
    delete []comp;
    delete []inv;
    delete []x;
    return proba;
}

double Cayley::get_theta_log_likelihood(int m, int *x_acumul, int *x_acumul_variation, double *theta_estim){
    double likeli= 0;
    double *psi;
    psi = new  double[ n_ ];
    theta_estim[n_-1] = 0;
    for(int i = 0;i < n_ - 1; i ++ ){
        int x_i = x_acumul[ i ];
        if(x_acumul_variation != NULL)x_i += x_acumul_variation[ i ];
        int j = i+1;
        if(x_i == 0)  x_i = 1;
        if(x_i == m)  x_i = m-1;
        double xav = (double) x_i/m;
        if(xav != 0){
            theta_estim[ i ] =log(n_-j)-log(xav / (1-xav));
            //if(theta_estim[ i ] < 0)theta_estim[ i ] = 0;
            psi[ i ] =1+(n_-j)*exp(-theta_estim[ i ]);
            likeli += x_i * theta_estim[ i ] + m * log(psi[ i ]);
            if ( likeli != likeli)
                likeli= 0;
        }else theta_estim[ i ] = 0;
    }
    delete [] psi;
    return (-1 * likeli);
}
void Cayley::get_x_lower_bound_freq(int m, int ** samples_inv_freq, int ini_pos, int *min_bound_x){
  //min bound
    int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
    int maxFreq= 0;
    for(int j = ini_pos ; j < n_ - 1 ; j++){
        for(int s = 0 ; s < n_ ; s++){
            freq[ s ] += samples_inv_freq[ j ][ s ];
            if(freq[ s ]> maxFreq)
                maxFreq=freq[ s ];
        }
        min_bound_x[ j ] = m - maxFreq;
        if(min_bound_x[ j ]<0)min_bound_x[ j ] = 0;
    }
    delete []freq;
}


double Cayley::get_bound_likeli(int m, int ** samples_inv_freq, int ini_pos, int* x , int*sigma_0){
  double likelihood = 0;
  int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
  int min_freq_j; //likeli upper bound
  int min_freq_tail;//likeli lower bound
  //min_freq_tail es el min_de_[j,n] que es lo mismo que ( m - max_de_[k,j) )
  int max_freq_middle= 0;
  
  double li_min, li_max,j;
  double xav, theta_estim, psi, likeli;
  int x_i;
  //likeli prev of the ones known
  
  for (int i = 0 ; i < ini_pos && i < n_-1; i++){
    x_i=x[i];
    if(x_i == 0)  x_i = 1;
    if(x_i == m)  x_i = m-1;
    j=i+1;
    xav = (double) x_i/m;
    theta_estim = log(n_-j)-log(xav / (1-xav));
    psi =1+(n_-j)*exp(-theta_estim);
    likeli = x[i] * theta_estim + m * log(psi);
    likelihood -= likeli;
  }
  //likeli of the bound
  for(int i = ini_pos ; i < n_ - 1 ; i++){
    int j = i+1;
    for(int s = 0 ; s < n_ ; s++)
      freq[ s ] += samples_inv_freq[ i ][ s ];
    
    min_freq_j = m;
    for(int s = 0 ; s < n_ ; s++){
      if( freq[ s ] > max_freq_middle) max_freq_middle = freq[ s ];// item s repeated max_freq_middle times in
      if( sigma_0[ s ] < 0 && samples_inv_freq[ i ][ s ] < min_freq_j ) min_freq_j = samples_inv_freq[ i ][ s ];
    }
    if( min_freq_j == 0 ) min_freq_j = 1;
    if( max_freq_middle == m ) max_freq_middle = m - 1;
    min_freq_tail = m - max_freq_middle;
    
    
    //the likelihood on the min(x), decreasing part of L
    xav = (double) min_freq_tail/m;
    theta_estim = log(n_-j)-log(xav / (1-xav));
    psi =1+(n_-j)*exp(-theta_estim);
    likeli = min_freq_tail * theta_estim + m * log(psi);
    li_min = -likeli;
    
    //the likelihood on the min(x), increasing part of L
    xav = (double) 1 - (double) min_freq_j/m;
    theta_estim = log(n_-j)-log(xav / (1-xav));
    psi =1+(n_-j)*exp(-theta_estim);
    likeli = min_freq_j * theta_estim + m * log(psi);
    if ( likeli != likeli)
      likeli= 0;
    li_max = -likeli;
    
    //cout<<(double) min_freq_j/m<<" "<<(double) min_freq_tail/m<<" - "<<min_freq_j<<" "<<min_freq_tail<<" - "<<li_max<<" "<<li_min<<endl;
    
    likelihood += li_min > li_max ? li_min : li_max;
  }
  delete []freq;
  return likelihood;
}


void Cayley::get_x_lower_bound(int m, int ** sample, int ini_pos, int *x_min_bound){
    int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
    int max_freq= 0;
    for(int j=ini_pos;j<n_-1;j++){
        for(int s= 0;s<m;s++){
            freq[sample[ s ][ j ]-1]++;
            if(freq[sample[ s ][ j ]-1]> max_freq)
                max_freq=freq[sample[ s ][ j ]-1];
        }
        x_min_bound[ j ] =m-max_freq;
        if(x_min_bound[ j ]<0)x_min_bound[ j ] = 0;
    }
    //  for (int i = 0 ; i < n; i ++ ) cout<<xMinBound[ i ]<<" ";cout<<" xbound, pos: "<<iniPos<<endl;
    delete []freq;
}

/* returns 2 vectors. item cycle_items[ i ] belongs to cycle cycle_indices[ i ]*/
int Cayley::get_cycles(int *sigma, int *cycle_items, int *cycle_indices){
    bool*visited = new bool[ n_ ];
    for (int i = 0 ; i < n_; i ++ )visited[ i ] =false;
    int item_index, cont=0, cycle_index= 0;
    while( cont < n_ ) {
        item_index = 0;
        while(visited[item_index])item_index++;
        while (!visited[item_index]) {
            visited[item_index] = true;
            cycle_items[cont] =item_index+1;
            cycle_indices[cont] =cycle_index;
            item_index=sigma[item_index]-1;
            cont++;
        }
        cycle_index++;
    }
    delete [] visited;
    return cycle_index;
}

void Cayley::random_sample_at_dist(int d, int m, int **samples){
    //samples must be initialized before as samples = new int*[ m ];
    for (int i = 0 ; i < m; i ++ ){
        samples[ i ] = new int[ n_ ];
        generate_permu_with_k_cycles(n_, n_ - d , samples[ i ]);
    }
}
double Cayley::calculate_psi(double *theta, double *psi_vector){
    double psi = 1;
    for (int i = 0 ; i <  n_ - 1; i ++ ){
        int j = i + 1;
        psi_vector[ i ] = 1 + ( n_ - j ) *(double) exp(-theta[ i ]);
        psi *= psi_vector[ i ];
    }
    return psi;
}


void Cayley::generate_permu_with_k_cycles(int n, int k, int * sigma){
    // n: 1..n_
    bool * long_cycle = new bool [ n_ ];
    int ran2;
    long double  ran1 ;
    
    while ( k > 1 ){
        ran1 =(long double) unif_rand();
        //ran1 = ((long double)rand() / ((double)RAND_MAX + 1 ) );
        if(ran1 < (long double)(stirling_matrix_[ n - 1 ][ k - 1 ] / stirling_matrix_[ n ][ k ])){
            long_cycle[ n - 1 ] = false;
            k --;
        }else{
            long_cycle[ n - 1 ] = true;
        }
        n --;
    }
    Generic gen; //use sigma_inv_ as auxiliary array
    //the n first items form a cycle
    gen.generate_random_permutation(n, 0, sigma_inv_);
    for (int i = 0 ; i < n-1; i ++ ) sigma[sigma_inv_[ i ]] = sigma_inv_[i+1]+1;
    sigma[sigma_inv_[n-1]] =sigma_inv_[0]+1;
    
    for (int i = n ; i < n_ ; i++){
        if ( long_cycle[ i ] ){
            //ran2 = rand() % ( i ) ;//0..n-2
            ran2 = (int) (unif_rand() * i );//[ 0,i)
            sigma[ i ] = sigma[ ran2 ];
            sigma[ ran2 ] = i + 1;
        }else{
            sigma[ i ] = i + 1;
        }
    }
    delete [] long_cycle;
}

void Cayley::dist_decomp_vector2perm(int* vec, int* sigma) {
    x_vector_to_permutation_forward(vec , sigma);
    //x_vector_to_permutation_backwards(vec , sigma);
}
void Cayley::x_vector_to_permutation_forward(int *x, int *sigma){
    int random, aux;
    for ( int i = 0 ; i < n_ ; i++) sigma[ i ] = i + 1;
    for ( int i = 0 ; i < n_ - 1 ; i++)
        if ( x [ i ] == 1 ){
            //random = (i + 1) + ( rand() % (n_ - i - 2 + 1) ); // random \in (i + 1)+[0..n-i-2] = [i+1..n-1]
            random = (i + 1) +  (int)( unif_rand() * (n_ - i - 2 + 1) );
            aux = sigma[  random ];
            sigma[ random ] = sigma[ i ];
            sigma[ i ] = aux;
        }        
}


void Cayley::x_vector_to_permutation_backwards(int *x, int *sigma){
    int     dist =  0;
    int     tables_num = 0 ;
    for (int i = 0 ; i < n_ ; i ++) dist += x[ i ];
    int     * tables_len = new int [ n_ - dist ];
    int     ** tables    = new int*[ n_ - dist ];
    for ( int i = 0 ; i < n_ - dist ; i++) {
        tables_len[ i ] = 1;
        tables[ i ] = new int[ dist + 1 ];
        for (int j = 0 ; j < dist ; j ++) tables[ i ][ j ] = 0 ;
    }
    
    x[ n_ - 1 ] = 0 ; //
    for ( int i = n_ - 1 ; i >= 0 ; i--){
        if ( x[ i ] == 0 ){ //sit the item at a new table
            tables[ tables_num ][ 0 ] = i;
            tables_num ++;
        }else{//unformly at random choose an item i among the ones that are sitting
            //and sit the new item at the right of i
            int table_num_insert = 0, table_pos_insert = 0;
            //table_pos_insert = (int) rand() % ( n_ - i - 1 );
            table_pos_insert = (int) ( unif_rand() * ( n_ - i - 1 ));
            while ( table_pos_insert >= tables_len[ table_num_insert ] ){
                table_pos_insert -= tables_len[ table_num_insert ] ;
                table_num_insert ++;
            }
            tables[ table_num_insert ][ tables_len[table_num_insert] ] = i;
            tables_len [table_num_insert ] ++;
        }
    }
    Generic gen;
    
    for (int i = 0 ; i < tables_num ; i++ ){
        gen.random_shuffle(tables_len[ i ], tables[ i ]);
        sigma[ tables[ i ][ tables_len[ i ] - 1 ]  ] = tables[ i ][ 0 ]  + FIRST_ITEM;
        for (int j = 0 ; j < tables_len[ i ] - 1 ; j++)
            sigma[ tables[ i ][ j ] ] = tables[ i ][ j + 1 ] + FIRST_ITEM;
  
        delete [] tables[ i ];
    }
    delete [] tables;
    delete [] tables_len;
    
}



long double Cayley::num_permus_at_distance(int d){
    return stirling_matrix_[n_ ][ n_ - d];
}
long double Cayley::count_permus_with_cycles(int d){
    return stirling_matrix_[n_ ][ n_ - d];
}


int Cayley::distance(int * s, int * t){
    int *comp = new int[ n_ ], *sigma_inv = new int[ n_ ];
    Generic gen;
    for(int j = 0 ; j < n_ ; j++) sigma_inv[t[ j ]-1] =j + 1;
    for(int i = 0 ; i < n_ ; i++) comp[ i ] = s [ sigma_inv [ i ] - 1 ];
    int dist = perm2dist_decomp_vector(comp, NULL);
    delete [] sigma_inv;
    delete [] comp;
    return dist;
}

int Cayley::perm2dist_decomp_vector(int*sigma, int*vec ){
    //x is a vector of length n
    //also updates the x vector if it isnot null
    if(vec!=NULL)for (int i = 0 ; i < n_; i ++ )vec[ i ] =1;
    //for (int i = 0 ; i < n; i ++ )cout<<sigma[ i ]<<" ";cout<<endl;
    int num_cycles=0, num_visited=0, item= 0;
    bool*visited = new bool[ n_ ];
    for (int i = 0 ; i < n_; i ++ )visited[ i ] =false;
    while(num_visited < n_ ){
        item = num_cycles;
        while ( visited[ item ]) item++;
        num_cycles++;
        int max_item_in_cycle= 0;
        do{
            if ( item > max_item_in_cycle ) max_item_in_cycle = item;
            visited[ item ] =true;
            num_visited++;
            item = sigma[ item ]-1;
        }while ( !visited[ item ] );
        if(vec != NULL )vec[ max_item_in_cycle ] = 0;
    }
    delete [] visited;
    return (n_ - num_cycles );
}


void Cayley::distances_sampling(int m, double theta, int **samples) {
    int     target_dist = 0;
    double  rand_val;
    long double *acumul = new long double[ n_ ];//+1
    acumul[0] =exp(-theta  * 0) * stirling_matrix_[ n_ ][ n_ ];//stirling1(n_,n_);
    for(int dista = 1 ; dista < n_ ; dista++)
        acumul[dista] =acumul[ dista - 1 ] +  exp(-theta  * dista) * stirling_matrix_[ n_ ][ n_ - dista ];//stirling1(n_,n_-dista);
    for (int i = 0 ; i < m; i ++ ){
        target_dist = 0;
        //rand_val = (double) acumul[ n_ - 1 ] * (double) rand() / RAND_MAX;
        rand_val = (double) (acumul[ n_ - 1 ] * unif_rand());
        while(acumul[target_dist] <= rand_val) target_dist++;
        //int *sigma=generate_permu_with_k_cycles(n_,(n_-target_dist));
        samples[ i ] = new int [ n_ ];
        generate_permu_with_k_cycles(n_,(n_-target_dist), samples[ i ]);
    }
    delete []acumul;
}



void Cayley::multistage_sampling(int m, double *theta, int **samples){
    double *psi = new double[n_-1];
    int *x = new int[ n_ ];
    Generic gen;
    
    calculate_psi(theta, psi);
    for(int samp = 0 ; samp < m ; samp++){
        for(int i= 0;i < n_ - 1; i ++ ){
            double probaX_j = (double) 1/psi[ i ];//(double) exp(-theta[ i ])/psi[ i ]; //
            //if(((double)rand() / RAND_MAX) < probaX_j) x[ i ] = 0;
            if(( unif_rand() ) < probaX_j) x[ i ] = 0;
            else x[ i ] =1;
        }
        x[n_-1] = 0;
        int *sigma = new int[ n_ ];
        dist_decomp_vector2perm(x, sigma);
        samples[samp] = sigma;
    }
    delete [] x;
    delete [] psi;
}


bool Cayley::same_cycle(int i, int j, int *sigma){
    int index = sigma[ i ] - 1;
    while(index != i && index != j)index = sigma[index]-1;
    if(index == j) return true;
    return false;
}

void Cayley::get_max_item_in_future_cycles(int *sigma, int i, int j, int *max_i, int *max_j){
    
    int pos=sigma[ i ]-1;
    *max_j = pos;

    while(pos != j ){
        pos  = sigma[pos]-1;
        if(pos > *max_j ) *max_j = pos;
    }
    pos=sigma[ j ]-1;
    *max_i = pos;
    while(pos != i){
        pos  = sigma[pos]-1;
        if(pos > *max_i ) *max_i = pos;
    }
}

void Cayley::get_max_item_in_current_cycle(int *sigma, int i, int *max_i ){
    int pos=sigma[ i ]-1;
    *max_i = pos;
    while(pos != i ){
        pos  = sigma[pos]-1;
        if(pos > *max_i ) *max_i = pos;
    }
}


void Cayley::gibbs_sampling(int m, double *theta, int model, int **samples) {
    int *sigma  = new int[ n_ ];
    Generic     gen;
    int burning_period_samples = n_*log(n_);
    gen.generate_random_permutation(n_, 1, sigma);
    
    for(int sample= 0;sample<m+burning_period_samples;sample++){
        int i,j, max_i=-1, max_j=-1, min;
        do{
            i = (int) (unif_rand() * n_); //rand() % n_;
            j = (int) (unif_rand() * n_); //rand() % n_;
        }while(i == j);
        bool make_swap = false;
        if(  same_cycle(i, j, sigma) )  make_swap=true;
        else{
            //double rand_double = (double)rand()/RAND_MAX;
            double rand_double = unif_rand();
            if(model == MALLOWS_MODEL){
                if(rand_double < exp(-theta[0])) make_swap = true;
            }else{
                get_max_item_in_current_cycle(sigma, i, &max_i);
                get_max_item_in_current_cycle(sigma, j, &max_j);
                min = (max_i < max_j) ? max_i : max_j;
                if(rand_double < exp(-theta[min])) make_swap = true;
            }
        }
        if(make_swap){
            int aux=sigma[ i ];
            sigma[ i ] =sigma[ j ];
            sigma[ j ] =aux;
        }
        //gen.print_int_vector(sigma, n_);
        if(sample>=burning_period_samples){
            samples[sample-burning_period_samples]  = new int[ n_ ];
            for (int i = 0 ; i < n_; i ++ )   samples[sample-burning_period_samples][ i ] =sigma[ i ];
        }
    }
    delete [] sigma;
}


long double Cayley::get_likelihood(int m, int** samples, int model, int * sigma_0){
    double  psi;
    int     dist    = 0 ;
    double  * theta = new double[ n_ ];
    long double loglikelihood;

    if(model == MALLOWS_MODEL){
        Newton_raphson newton(n_);
        double  *psi_vec = new double [n_];
        dist = distance_to_sample(samples, m, sigma_0);
        theta[ 0 ] = newton.Newton_raphson_method((double)dist/m, -10.001,CAYLEY_DISTANCE, MALLOWS_MODEL, -1, NULL);
         for (int i = 1 ; i < n_ -1; i++) theta[ i ] = theta[ 0 ];
         psi = calculate_psi(theta, psi_vec);
         delete [] psi_vec;
         loglikelihood = - theta[0] * dist - m * log (psi);
    }else{
        int *x = new int[ n_ ], *x_acumul = new int[ n_ ], *inv = new int[ n_ ], *comp = new int[ n_ ];
        for(int i = 0 ; i < n_ ; i++) x_acumul[ i ] = 0;
        for(int i = 0 ; i < n_ ; i++) inv [ sigma_0[ i ] - 1 ] = i + 1;
        for (int i = 0 ; i < m; i ++ ){
            for(int j = 0 ; j < n_ ; j ++) comp[ j ] = samples[ i ] [ inv [ j ] - 1 ];
            perm2dist_decomp_vector(comp, x);
            for(int j = 0 ; j < n_ - 1 ; j++) x_acumul[ j ] += x[ j ];
        }
        loglikelihood = get_theta_log_likelihood(m, x_acumul, NULL, theta);//theta is an array of length n
        delete [] x;
        delete [] x_acumul;
        delete [] inv;
        delete [] comp;
    }
    delete [] theta;
    return loglikelihood;
}


 void Cayley::estimate_theta(int m, int *sigma_0, int **samples, int model, double *theta){
     if(model == MALLOWS_MODEL){
        int dist = distance_to_sample(samples, m, sigma_0);
        Newton_raphson newton(n_);
        theta[0] = newton.Newton_raphson_method((double)dist/m, -10.001,CAYLEY_DISTANCE, MALLOWS_MODEL, -1, NULL);
    }else{
        int *x = new int[ n_ ], *x_acumul = new int[ n_ ], *inv = new int[ n_ ], *comp = new int[ n_ ];
        for(int i = 0 ; i < n_ ; i++) x_acumul[ i ] = 0;
        for(int i = 0 ; i < n_ ; i++) inv [ sigma_0[ i ] - 1 ] = i + 1;
        for (int i = 0 ; i < m; i ++ ){
            for(int j = 0 ; j < n_ ; j ++) comp[ j ] = samples[ i ] [ inv [ j ] - 1 ];
            perm2dist_decomp_vector(comp, x);
            for(int j = 0 ; j < n_ - 1 ; j++) x_acumul[ j ] += x[ j ];
        }
        get_theta_log_likelihood(m, x_acumul, NULL, theta);//theta is an array of length n
        //cout<<"likeli(Cayley::estimate_theta) "<<get_theta_log_likelihood(m, x_acumul, NULL, theta)<<endl;
        delete [] x;
        delete [] x_acumul;
        delete [] inv;
        delete [] comp;
    }
}

int Cayley::distance_to_sample(int **samples, int m, int *sigma){
    int distance= 0;
    int *comp = new int[ n_ ], *sigma_inv = new int[ n_ ];
    for(int j = 0 ; j < n_ ; j ++) sigma_inv[sigma[ j ] - 1 ] = j + 1;
    for(int s = 0 ; s < m ; s ++){
        for(int i = 0 ; i < n_ ; i ++) comp[ i ] = samples[ s ][ sigma_inv [ i ] - 1 ];
        distance += perm2dist_decomp_vector(comp, NULL);
    }
    delete []sigma_inv;
    delete []comp;
    return distance ;
}

long double Cayley::count_permus_by_x(int *x){
    int d = 0, c_len = 0;
    for(int i = 0 ; i < n_ ; i++) if(x[ i ] == 1 ) d++;
    c_len = n_ - d;
    long double *c = new long double[c_len];
    for ( int i = 0 ; i < c_len ; i ++) c[ i ] = 0 ;
    int *k= new int[ n_ ], k_in=c_len;
    for( int i = n_ - 1 ; i >= 0 ; i --){// first top cycle item with which can be joind¡ed in the cycle
        if ( x[ i ] == 0 )k_in --;
        k[ i ] = k_in;
    }
    long double l = count_permus_by_x_core(0, x, c, c_len, k);
    delete [] c;
    delete [] k;
    return l;
}

long double Cayley::count_permus_by_x_core(int index, int *x, long double *c, int c_len, int *k){
    Generic gen;
    long double l = 0, p = 1;
    if( index == n_ ) {
        for(int i = 0 ; i < c_len ; i ++) p *= gen.factorial(c[ i ]);
        return p ;
    }
    if (x[index] == 0) return count_permus_by_x_core(index + 1, x, c , c_len, k);
    for(int i = k[index] ; i < c_len ; i ++){
        long double *c_i = new long double[c_len];
        for(int j = 0 ; j < c_len ; j ++) c_i[ j ] = c[ j ];
        c_i[ i ] ++;
        l += count_permus_by_x_core(index+1, x, c_i, c_len, k);
        delete [] c_i;
    }
    return l;
}

double Cayley::estimate_consensus_exact_gmm(int m, int **samples, int*sigma_0_ini, int *sigma_0){
    int     ** samples_inv = new int*[ m ];
    int     ** samples_inv_freq = new int*[ n_ ];
    int     * x_acum = new int[ n_ ];
    int     * sigma_0_aux = new int[ n_ ];
    int     * sigma_0_inv_aux = new int[ n_ ];
    double  likelihood = -DBL_MAX;

    
    Generic gen;
    for (int i = 0 ; i < m; i ++ ){
        samples_inv[ i ]  = new int[ n_ ];
        gen.invert(n_, samples[ i ], samples_inv[ i ]);
    }
    for (int i = 0 ; i < n_; i ++ ){ //ini
        samples_inv_freq[ i ]  = new int[ n_ ];
        for(int j= 0;j<n_;j++)samples_inv_freq[ i ][ j ] = 0;
        x_acum[ i ] = 0;
    }
    for(int i = 0 ; i < m; i ++ )for(int j= 0;j<n_;j++)samples_inv_freq[ j ][samples_inv[ i ][ j ]-1]++;
    for (int i = 0 ; i < n_; i ++ ) {sigma_0_aux[ i ] = -1 ; sigma_0_inv_aux[ i ] = -1 ;}
    
    //bounds
    estimate_consensus_approx(GENERALIZED_MALLOWS_MODEL , m, samples, sigma_0);
    likelihood = get_likelihood(m , samples , GENERALIZED_MALLOWS_MODEL , sigma_0);
    
    if ( sigma_0_ini != NULL ){
        //obtain the likelihod of the proposed sigma_0_ini.
        double like_ini = get_likelihood(m , samples , GENERALIZED_MALLOWS_MODEL , sigma_0_ini);
        //check both solutions: the proposed sigma_0_ini vs. the approx. Discard the worse
        if (like_ini > likelihood) {
            likelihood = like_ini;
            for (int i = 0 ; i < n_ ; i++) sigma_0[ i ] = sigma_0_ini[ i ];
        }
  }
    
    
    double visited_nodes = estimate_consensus_exact_gmm_core(m, 0 , samples, samples_inv, samples_inv_freq,
                                                   x_acum, sigma_0_aux, sigma_0_inv_aux, 0, sigma_0, &likelihood);
    
    //delete
    for (int i = 0 ; i < m; i ++ ) delete []samples_inv[ i ];
    delete [] samples_inv;
    for (int i = 0 ; i < n_; i ++ ) delete []samples_inv_freq[ i ];
    delete [] samples_inv_freq;
    delete [] x_acum;
    delete [] sigma_0_aux;
    delete [] sigma_0_inv_aux;
    return visited_nodes;
}

double Cayley::estimate_consensus_exact_gmm_core(int m, int pos, int ** samples, int ** samples_inv, int **samples_inv_freq,
                                       int * x_acum, int * current_sigma, int * current_sigma_inv, double current_likeli_bound,
                                       int * best_sigma, double * best_likeli){
  //for (int i = 0 ; i < n_; i ++ ) cout<<current_sigma[ i ]; cout<<endl;
    if(pos == n_ && current_likeli_bound >= (*best_likeli)){
        for (int i = 0 ; i < n_; i ++ )best_sigma[ i ] =current_sigma[ i ];
        (*best_likeli)=current_likeli_bound;
        return 1;
    }
    double visited_nodes= 0;
    bool trace=false;
    int *x_acum_var = new int[ n_ ];
    double *theta_estim = new double[ n_ ];
    for(int it= 0;it<n_;it++){//sigmaInv(pos)=it
        if(current_sigma[it] ==-1){
            int xIncr= 0;
            current_sigma_inv[pos] =it+1;
            current_sigma[it] =pos+1;
            int *pos_swaps = new int[ m ];
            for(int s= 0;s<m;s++){
                pos_swaps[ s ] =-1;
                if(samples[ s ][it] != current_sigma[it]){
                    int x= samples[ s ][it];
                    int y= samples_inv[ s ][pos]-1;
                    samples[ s ][it] =pos+1;
                    samples[ s ][y] =x;
                    samples_inv[ s ][pos] =it+1;
                    samples_inv[ s ][x-1] =y+1;
                    pos_swaps[ s ] =y;
                    xIncr++;
                    samples_inv_freq[pos][y]--;
                    samples_inv_freq[x-1][it]--;
                    samples_inv_freq[pos][it]++;
                    samples_inv_freq[x-1][y]++;
                }
                /*if(trace){
                    for (int k = 0 ; k < n_ ; k ++ ) cout<<samples[ s ][ k ]<<" ";cout<<" sample ";
                    for (int k = 0 ; k < n_ ; k ++ ) cout<<samples_inv[ s ][ k ]<<" ";cout<<" sampInv ";
                    for (int i = 0 ; i < n_ ; i ++ ) cout<<current_sigma_inv[ i ]<<" ";cout<<" sigmaInv ";
                    for (int i = 0 ; i < n_ ; i ++ ) cout<<current_sigma[ i ]<<" ";cout<<" sigma (v3)"<<endl;;
                }*/
                //in this case the distance has incresed in one, x[maxItemInCylce] was 0 and now = 1
            }
            
            for (int i = 0 ; i < n_; i ++ )x_acum_var[ i ] =x_acum[ i ];
            x_acum_var[pos] += xIncr;
            double likeliBound = get_bound_likeli(m, samples_inv_freq, pos+1, x_acum_var, current_sigma);
            
            if(likeliBound >= (*best_likeli) )
                visited_nodes += estimate_consensus_exact_gmm_core(m, pos+1, samples, samples_inv, samples_inv_freq, x_acum_var, current_sigma, current_sigma_inv, likeliBound, best_sigma, best_likeli);
           // else {cout<<"bounded at "<<pos<<" - ";for (int i = 0 ; i < n_; i ++ )cout<<current_sigma[ i ]<<" ";cout<<endl;}
            current_sigma_inv[pos] =-1;
            current_sigma[it] =-1;
            for (int s = 0 ; s < m ; s ++){
                if(pos_swaps[ s ] != -1){
                    int y=pos_swaps[ s ];
                    int x=samples[ s ][y];
                    samples[ s ][it] =x;
                    samples[ s ][y] =pos+1;
                    samples_inv[ s ][pos] =y+1;
                    samples_inv[ s ][x-1] =it+1;
                    pos_swaps[ s ] =-1;
                    samples_inv_freq[pos][y]++;
                    samples_inv_freq[x-1][it]++;
                    samples_inv_freq[pos][it]--;
                    samples_inv_freq[x-1][y]--;
                }
            }
            delete []pos_swaps;
        }
    }
    delete [] theta_estim;
    delete [] x_acum_var;
    return visited_nodes+1;
}

double Cayley::estimate_consensus_exact_mm(int m, int **samples, int*sigma_0_ini, int *sigma_0){
    int **samples_inv = new int*[ m ];
    int *x_acum = new int[ n_ ];
    int *sigma_0_aux = new int[ n_ ];
    int *sigma_0_inv_aux = new int[ n_ ];
    Generic gen;
    for (int i = 0 ; i < m; i ++ ){
        samples_inv[ i ]  = new int[ n_ ];
        gen.invert(n_, samples[ i ], samples_inv[ i ]);
    }
    double best_distance = (n_ - 1) * m;//maximum
    for (int i = 0 ; i < n_; i ++ ) {
        sigma_0_aux[ i ] =-1;
        sigma_0_inv_aux[ i ] =-1;
        x_acum[ i ] = 0;
    }
    //bounds
    estimate_consensus_approx(MALLOWS_MODEL , m, samples, sigma_0);
    best_distance = distance_to_sample(samples , m , sigma_0);
    
    if ( sigma_0_ini != NULL ){
        //obtain the likelihod of the proposed sigma_0_ini.
        double dist_ini = distance_to_sample(samples , m , sigma_0_ini);
        //check both solutions: the proposed sigma_0_ini vs. the approx. Discard the worse
        if (dist_ini < best_distance) {
            best_distance = dist_ini;
            for (int i = 0 ; i < n_ ; i++) sigma_0[ i ] = sigma_0_ini[ i ];
        }
    }

    double visited_nodes = estimate_consensus_exact_mm_core(m, 0 , samples, samples_inv, x_acum, sigma_0_aux, sigma_0_inv_aux, 0, sigma_0, &best_distance);
    for (int i = 0 ; i < m; i ++ ) delete []samples_inv[ i ];
    delete [] samples_inv;
    delete [] x_acum;
    delete [] sigma_0_aux;
    delete [] sigma_0_inv_aux;
    return visited_nodes;
}
double Cayley::estimate_consensus_exact_mm_core(int m, int pos, int ** samples, int ** samples_inv, int * x_acum,  int * current_sigma,
                              int * current_sigma_inv, double current_dist_bound, int *best_sigma, double *best_dist){
    if(pos == n_ && current_dist_bound <= (*best_dist)){
        for (int i = 0 ; i < n_; i ++ )best_sigma[ i ] =current_sigma[ i ];
        (*best_dist)=current_dist_bound;
        return 1;
    }
    double visited_nodes= 0;
    bool trace=true, enc=false;
    int *x_acum_var = new int[ n_ ], *candVec = new int[ n_ ];
    int cand=n_;
    int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
    for(int s = 0;(!enc && s<m );s++){
        if((freq[samples_inv[ s ][pos]-1]++) > (m/2)){
            candVec[0] =samples_inv[ s ][pos]-1;
            enc=true;
            cand=1;
        }
    }
    if(!enc){
        for(int it = 0 ; it < n_ ; it ++) candVec[ it ] = it;
        cand=n_;
    }
    for(int index = 0;index < cand ; index++){//sigmaInv(pos)=it
        int it = candVec[index];
        if(current_sigma[it] ==-1){
            int x_incr= 0;
            current_sigma_inv[pos] =it+1;
            current_sigma[it] =pos+1;
            int *pos_swaps = new int[ m ];
            //for (int i = 0 ; i < m; i ++ )posSwaps[ i ] =-1;
            for(int s= 0;s<m;s++){
                pos_swaps[ s ] =-1;
                if(samples[ s ][it] != current_sigma[it]){
                    int x = samples[ s ][it];
                    int y = samples_inv[ s ][pos]-1;
                    samples[ s ][it] = pos+1;
                    samples[ s ][y] = x;
                    samples_inv[ s ][pos] = it + 1;
                    samples_inv[ s ][x-1] = y + 1;
                    pos_swaps[ s ] =y;
                    x_incr++;
                }
               /* if(trace)
                {   for (int k = 0 ; k < n_ ; k ++) cout<<samples[ s ][ k ]<<" ";cout<<" sample ";
                    for (int k = 0 ; k < n_ ; k ++) cout<<samples_inv[ s ][ k ]<<" ";cout<<" sampInv ";
                    for (int i = 0 ; i < n_; i ++ ) cout<<current_sigma_inv[ i ]<<" ";cout<<" sigmaInv ";
                    for (int i = 0 ; i < n_; i ++ ) cout<<current_sigma[ i ]<<" ";cout<<" sigma (v3 MM)"<<endl;;
                }*/
                //in this case the distance has incresed in one, x[maxItemInCylce] was 0 and now = 1
            }
            
            double distance_bound = 0;
            int *xbound = new int[ n_ ];for (int i = 0 ; i < n_; i ++ )xbound[ i ] = 0;
            get_x_lower_bound(m, samples_inv, pos+1, xbound);
            //if(trace){for (int i = 0 ; i < n_; i ++ ) cout<<xbound[ i ]<<" ";cout<<" xbound, pos: "<<pos<<endl;}
            for(int i=pos+1 ; i < n_; i ++ ) distance_bound+=xbound[ i ];
            delete []xbound;
            for (int i = 0 ; i < n_; i ++ ){
                x_acum_var[ i ] =x_acum[ i ];
                distance_bound += x_acum_var[ i ];
            }
            x_acum_var[pos] += x_incr;
            distance_bound += x_incr;
            
            if(distance_bound <= (*best_dist)) visited_nodes+=estimate_consensus_exact_mm_core(m, pos+1, samples, samples_inv, x_acum_var, current_sigma, current_sigma_inv, distance_bound, best_sigma, best_dist);
            //else {cout<<"bounded at "<<pos<<" - ";for (int i = 0 ; i < n; i ++ )cout<<current_sigma[ i ]<<" ";cout<<endl;}
            current_sigma_inv[pos] =-1;
            current_sigma[it] =-1;
            for(int s= 0;s<m;s++){
                if(pos_swaps[ s ] != -1){
                    int y = pos_swaps[ s ];
                    int x = samples[ s ][y];
                    samples[ s ][it] =x;
                    samples[ s ][y] =pos+1;
                    samples_inv[ s ][pos] =y+1;
                    samples_inv[ s ][x-1] =it+1;
                    pos_swaps[ s ] =-1;
                }
            }
            delete [] pos_swaps;
        }
    }
    
    delete [] x_acum_var;
    delete [] candVec;
    delete [] freq;
    return visited_nodes + 1 ;
}

void Cayley::estimate_consensus_approx(int model, int m, int **samples, int *sigma_0){
    int ** samples_inv = new int*[ m ];
    int ** samples_copy = new int*[ m ];
    for (int i = 0 ; i < m; i++) {
        samples_copy[ i ] = new int[ n_ ];
        samples_inv [ i ] = new int[ n_ ];
        for(int j = 0 ; j < n_ ; j ++) {samples_inv[ i ][ samples[ i ][ j ] - 1 ] = j + 1;samples_copy[ i ][ j ] = samples[ i ][ j ];}
    }
    double bl_data;
    double *best_likeli = &bl_data;
    if( model == MALLOWS_MODEL )
        estimate_consensus_approx_mm(m, samples_copy, samples_inv, sigma_0, best_likeli);
    else
        estimate_consensus_approx_gmm(m, samples_copy, samples_inv, sigma_0, best_likeli);
    variable_neighborhood_search(m, samples, sigma_0, model , best_likeli);
    for (int i = 0 ; i < m; i++) {
        delete [] samples_inv[ i ];
        delete [] samples_copy[ i ];
    }
    delete [] samples_copy;
    delete [] samples_inv;
}

void Cayley::estimate_consensus_approx_gmm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_likeli){
    // pre: sigma, sigmaInv =-1
    int *freq = new int[ n_ ];     
    int *x_acum = new int[ n_ ];   for (int i = 0; i < n_; i++) x_acum[ i ] = 0;
    int * sigma_0_inv = new int[ n_ ];
    
    for (int i = 0; i < n_; i++){ sigma_0[ i ] = -1;sigma_0_inv[ i ] = -1;}

    for (int item = 0 ; item < n_ ; item ++){
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  1"<<endl;}
        for (int i = 0; i < n_; i++) freq[ i ] = 0;
        int max_freq = 0;
        int pj = -1; // 0..n-1
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  2"<<endl;}
        for (int s = 0;  s < m; s++){
            freq[ samples_inv[ s ][item] - 1] ++;
            if ((freq[ samples_inv[ s ][item] - 1]) > max_freq){
                max_freq = freq[ samples_inv[ s ][item] - 1];
                pj = samples_inv[ s ][item] - 1;
            }
        }
        //for (int i = 0 ; i < n; i ++ ){for(int j= 0;j<n;j++)cout<<freq[ i ][ j ]<<" ";cout<<endl;}
        sigma_0_inv[item] = pj + 1;
        sigma_0 [pj] = item + 1;
        for (int s = 0; s < m; s++) {
            if (samples_copy[ s ][pj] != item + 1) {//swap
                int x = samples_copy[ s ][pj];
                int y = samples_inv[ s ][item] - 1;
                samples_copy[ s ][pj] = item + 1;
                samples_copy[ s ][y] = x;
                samples_inv [ s ][item] = pj + 1;
                samples_inv [ s ][x - 1] = y + 1;
                x_acum[ item ]++;
            }
        }
    }
    
    double *theta = new double[ n_ ];
    *best_likeli = get_theta_log_likelihood(m, x_acum, NULL, theta);
    delete [] theta;
    delete [] x_acum;
    delete [] sigma_0_inv;
    delete [] freq;
}

void Cayley::estimate_consensus_approx_mm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_distance){
    // pre: sigma, sigmaInv =-1
    int distance_increase = 0, remaining = n_;
    int **freq = new int*[ n_ ];for (int i = 0; i < n_; i++) freq[ i ] = new int[ n_ ];
    
    for (int i = 0; i < n_; i++) sigma_0[ i ] = -1;
    do{ //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  1"<<endl;}
        for (int i = 0; i < n_; i++)
            for (int j = 0; j < n_; j++)freq[ i ][ j ] = 0;
        int max_freq = 0;
        int pi = -1, pj = -1; // 0..n-1
        bool dirty_items = true ; //si hay alguno q aparezca >m/2 veces?
        
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n_; i ++ )cout<<samples_copy[ s ][ i ]<<" ";cout<<"samples  trace  2"<<endl;} 
        //for (int i = 0 ; i < n_; i ++ ) cout<<sigma_0[i]<<" ";cout<<endl;
        Generic gen;
        //gen.print_int_matrix(samples_copy, m, n_);
        for (int s = 0; dirty_items && s < m; s++)
            for (int i = 0; ( dirty_items && i < n_ ) ; i++){
                //gen.print_int_matrix(freq, n_, n_);
                if (sigma_0[ i ] == -1) (freq[ i ][samples_copy[ s ][ i ] - 1])++;
                if ((freq[ i ][samples_copy[ s ][ i ] - 1]) > max_freq){
                    max_freq = (freq[ i ][samples_copy[ s ][ i ] - 1]);
                    pi = i;
                    pj = samples_copy[ s ][ i ] - 1;
                    if (max_freq > m/2 ) dirty_items = false ;
                }
            }
        //for (int i = 0 ; i < n; i ++ ){for (int j = 0 ; j < n;j++)cout<<freq[ i ][ j ]<<" ";cout<<endl;}
        sigma_0[pi] = pj + 1;
        for (int s = 0; s < m; s++) {
            if (samples_copy[ s ][pi] != pj + 1) {//////////////swap
                //swapSamplesDo(s, pi, pj, samples, samples_inv, posSwap);
                int x = samples_copy[ s ][pi];
                int y = samples_inv[ s ][pj] - 1;
                samples_copy[ s ][pi] = pj + 1;
                samples_copy[ s ][y] = x;
                samples_inv[ s ][pj] = pi + 1;
                samples_inv[ s ][x - 1] = y + 1;
                distance_increase++;
            }
        }
        remaining--;
    }while(remaining > 0) ;
    (*best_distance) = distance_increase;
    for (int i = 0 ; i < n_ ; i ++) delete [] freq[ i ];
    delete [] freq;
}

void Cayley::variable_neighborhood_search(int m, int **samples, int *sigma, int model, double *f_eval){
    bool improve;
    do{
        double f_eval_ini = (*f_eval);
        //cout<<" distance ini"<<(*f_eval)<<endl;
        if( model == MALLOWS_MODEL ) local_search_swap_mm  (m, samples, sigma, f_eval);
        else local_search_swap_gmm  (m, samples, sigma, f_eval);
        //cout<<" distance new swap "<<(*f_eval)<<endl;
        local_search_insert(m, samples, sigma, model, f_eval);
        //cout<<" distance new ins "<<(*f_eval)<<endl;
        improve=false ;
        if( model == MALLOWS_MODEL  && (f_eval_ini) > *f_eval ) improve = true;
        if( model == GENERALIZED_MALLOWS_MODEL  && (f_eval_ini) < *f_eval ) improve = true;
    }while(improve);
}

void Cayley::local_search_swap_mm(int m, int **samples, int *sigma_0, double *f_eval){
    int ** samples_comp = new int*[ m ];
    for (int s = 0 ; s < m; s++)  samples_comp[ s ] = new int[ n_ ];
    int *sigma_0_inv = new int[ n_ ];
    int *cycle_items = new int [n_ ], *cycle_index = new int[ n_];
    int **same_cycle = new int*[ n_ ];
    for (int i = 0 ; i < n_ ; i ++){ same_cycle[ i ] = new int[ n_ ]; for (int j = 0 ; j < n_ ; j++) same_cycle[ i][ j ] = 0;}
    Generic gen;
    int index_i,index_j, distance_variation = 0;
    for(int i = 0 ; i < n_ ; i ++) sigma_0_inv[ sigma_0 [ i ] - 1 ] = i + 1;
    bool improve;
    do {//iterating use just sigma_0_inv
        int max_freq = 0;
        for (int i = 0 ; i < n_ ; i ++) for (int j = 0 ; j < n_ ; j++) same_cycle[ i][ j ] = 0;
        for (int s = 0 ; s < m; s++) {
            for(int j = 0 ; j < n_ ; j ++ ) samples_comp[ s ][ j ] = samples[ s ][ sigma_0_inv [ j ] - 1 ];
            get_cycles(samples_comp[ s ], cycle_items, cycle_index);
            //cout<<"Cycles ";gen.print_int_vector(cycle_items, n_); gen.print_int_vector(cycle_index, n_);
            for ( int i = 0 ; i < n_ ; i ++){
                for (int j = i + 1; j < n_ && (cycle_index[ i ] == cycle_index[ j ]) ; j++){
                    int max, min; // for a triangular matrix
                    if (cycle_index [ i ] > cycle_index[ j ]){ max = cycle_items[ i ] - 1; min  = cycle_items[ j ]-1;}
                    else {min  = cycle_items[ i ]-1; max  = cycle_items[ j ]-1;}
                    same_cycle[ min ][ max ] ++;
                    if (max_freq < same_cycle[ min ][ max ] ) {
                        max_freq = same_cycle[ min ][ max ];
                        index_i = min;
                        index_j = max;
                    }
                    //if (max_freq > m)//ERROR
                    //{ cout<<"Cycles ";gen.print_int_vector(cycle_items, n_); gen.print_int_vector(cycle_index, n_);gen.print_int_matrix(same_cycle, n_, n_);}
                }
            }
        }
        //there are m oermus, max_freq of them are going to decrease dist in 1. (m-max_freq) are going to incresse
        distance_variation = m - 2 * max_freq;
        
        improve = false;
        if( distance_variation < 0 ) {
            improve = true;
            int aux = sigma_0_inv[index_i];
            sigma_0_inv[index_i] = sigma_0_inv[index_j];
            sigma_0_inv[index_j] = aux;
            
            (*f_eval) += distance_variation;
            for(int i = 0 ; i < n_ ; i ++) sigma_0[ sigma_0_inv[ i ] - 1 ] = i + 1;
            
        }
        
        //cout<<"trace in local search swap, distance variation "<<distance_variation<<" mxa freq "<<max_freq<<endl;
    } while (improve );
    
    delete [] sigma_0_inv;
    delete [] cycle_index;
    delete [] cycle_items;
    for (int i = 0 ; i < n_ ; i++) delete [] same_cycle[ i ];
    delete [] same_cycle;
    for (int i = 0 ; i < m ; i++)  delete [] samples_comp[ i ];
    delete [] samples_comp;
}

void Cayley::local_search_swap_gmm(int m, int **samples, int *sigma_0, double *f_eval){
    int ** samples_comp = new int*[ m ];
    for (int s = 0 ; s < m; s++)  samples_comp[ s ] = new int[ n_ ];
    int *sigma_0_inv = new int[ n_ ];
    int *x_acum = new int [n_ ], *x_acum_var = new int[ n_], *x = new int[ n_], *best_x_acum = new int[ n_ ], *x_test = new int[ n_ ];
    double *theta = new double [ n_ ];
    Generic gen;
    int index_i,index_j;
    for(int i = 0 ; i < n_ ; i ++) {sigma_0_inv[ sigma_0 [ i ] - 1 ] = i + 1;x_acum[ i ] = 0;x_acum_var[ i ] = 0;}
    bool improve;
    for (int s = 0 ; s < m; s++){
        for(int j = 0 ; j < n_ ; j ++ )
            samples_comp[ s ][ j ] = samples[ s ][ sigma_0_inv [ j ] - 1 ];
        perm2dist_decomp_vector(samples_comp[ s ], x);
        for(int j = 0 ; j < n_ ; j ++ ) x_acum[ j ] += x[ j ];
    }

    do {//iterating use just sigma_0_inv
        double best_like = 0;
        for (int iter_i = 0; iter_i < n_ - 1 ; iter_i++) {
            for (int iter_j = iter_i + 1 ; iter_j < n_ ; iter_j++) {
                //test swap i j on samples_comp
                for(int k = 0 ; k < n_ ; k ++) {x_acum_var[ k ] = 0;}
                for (int s = 0 ; s < m ; s++) {
                    int max_i, max_j;
                    //cout<<"sample "; Generic gen; gen.print_int_vector(samples_comp[ s ], n_);
                    if (same_cycle(iter_i, iter_j, samples_comp[ s ])){
                        get_max_item_in_future_cycles(samples_comp[ s ], iter_i, iter_j, &max_i , &max_j);
                        x_acum_var[(max_i < max_j ? max_i : max_j )] --;
                        //cout<<" resta en "<<(max_i < max_j ? max_i : max_j )<<endl;
                    }else {
                        get_max_item_in_current_cycle(samples_comp[ s ], iter_i, &max_i );
                        get_max_item_in_current_cycle(samples_comp[ s ], iter_j, &max_j );
                        x_acum_var[(max_i < max_j ? max_i : max_j )] ++;
                        //cout<<" suma en "<<(max_i < max_j ? max_i : max_j )<<endl;
                    }
                }
                //////test
                double like = get_theta_log_likelihood(m, x_acum, x_acum_var, theta);
                if ( like > best_like || best_like == 0){
                    for(int i = 0 ; i < n_ ; i ++) best_x_acum[ i ] = x_acum_var[ i ];
                    best_like = like;
                    index_i = iter_i;
                    index_j = iter_j;
                }
            }
        }
        //there are m oermus, max_freq of them are going to decrease dist in 1. (m-max_freq) are going to incresse
        improve = false;
        

        if( best_like > *f_eval ) {
            improve = true;
            (*f_eval) = best_like;
            int aux = sigma_0_inv[index_i];
            sigma_0_inv[index_i] = sigma_0_inv[index_j];
            sigma_0_inv[index_j] = aux;
            for(int i = 0 ; i < n_ ; i ++) sigma_0[ sigma_0_inv[ i ] - 1 ] = i + 1;
            for (int s = 0 ; s < m ; s++) {
                int aux = samples_comp[ s ][index_i];
                samples_comp [ s ][index_i] = samples_comp[ s ][index_j];
                samples_comp [ s ][index_j] = aux;
            }
            for(int j = 0 ; j < n_ ; j ++ ) x_acum[ j ] += best_x_acum[ j ];
            //test a ver si son iguales
        }
        //cout<<"trace in local search swap, distance variation "<<distance_variation<<" mxa freq "<<max_freq<<endl;
    } while (improve );


    delete [] theta;
    delete [] best_x_acum;
    delete [] sigma_0_inv;
    delete [] x_acum;
    delete [] x_acum_var;
    delete [] x;
    delete [] x_test;
    for (int i = 0 ; i < m ; i++)  delete [] samples_comp[ i ];
    delete [] samples_comp;
}

void Cayley::local_search_insert (int m, int **samples, int *sigma_0, int model, double *f_eval){
    int *best_sol = new int[ n_ ], *next_sol = new int[ n_ ], *x_acumul = new int[ n_ ], *x = new int[ n_ ];
    double * theta = new double[ n_ ];
    double best_eval= 0;
    bool better;
    Generic gen;
    do {
        better = false;
        best_eval = 0;
        for(int i = 0 ; i < n_ ; i ++){
            for ( int j = 0 ; j < n_ ; j++){
                if(i != j){
                    gen.insert_at (sigma_0, n_, i, j, next_sol);
                    if(model == MALLOWS_MODEL){
                        int dist = distance_to_sample(samples, m , next_sol);
                        if ( dist < best_eval || best_eval == 0){
                            best_eval = (double )dist;
                            for (int i = 0 ; i < n_ ; i ++) best_sol [ i ] = next_sol[ i ];
                        }
                    }else {
                        for(int index = 0 ; index < n_ ; index ++ ) x_acumul[ index ] = 0;
                        for (int s = 0 ; s < m; s++){
                            perm2dist_decomp_vector(samples[ s ], x);
                            for(int index = 0 ; index < n_ ; index ++ ) x_acumul[index] += x[index];
                        }
                        double likeli = get_theta_log_likelihood(m,x_acumul, NULL, theta);
                        if (likeli > best_eval || best_eval == 0) {
                            best_eval = likeli;
                            for (int i = 0 ; i < n_ ; i ++) best_sol [ i ] = next_sol[ i ];
                        }
                    }
                }
            }
        }
        
        if( (model == MALLOWS_MODEL && best_eval < *f_eval ) || ( model == GENERALIZED_MALLOWS_MODEL && best_eval > *f_eval ) ){
            *f_eval = best_eval;
            for (int i = 0 ; i < n_ ; i ++) sigma_0 [ i ] = best_sol[ i ];
            better = true;
        }
    } while (better);
    delete [] next_sol;
    delete [] theta;
    delete [] x;
    delete [] x_acumul;
    delete [] best_sol;
}
