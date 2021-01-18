//
//  Generic.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include "Generic.h"
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include "Hamming.h"
#include "Cayley.h"
#include "Ulam.h"
#include "Kendall.h"
#include "Ulam_disk.h"
#include <R.h>

Exponential_model* Generic::new_instance(int distance_id, int n){
    if ( distance_id == HAMMING_DISTANCE )   return (new Hamming( n ));
    if ( distance_id == CAYLEY_DISTANCE )    return new Cayley( n );
    if ( distance_id == ULAM_DISTANCE )      return new Ulam( n );
    if ( distance_id == ULAM_DISK_DISTANCE ) return new Ulam_disk( n );
    if ( distance_id == KENDALL_DISTANCE )   return new Kendall( n );
    return NULL;
}

void Generic::elementary_symmetric_polynomial(double* theta, int n, long double*theta_exp_aux, long double **esp_aux, long double *esp){
    //esp[j][n]: j-th elementarySymmetricPolynomials of n items
    //theta_exp_aux , esp_aux: are defined outside because the learning process (NewtonRaphson) calls this func lots of return;

    for ( int i = 0 ; i < n ; i ++ ){
        for ( int j = 0 ; j <= n ; j ++) esp_aux[i][j]=0;
        theta_exp_aux[ i + 1 ] = (long double)exp( theta[ i ]) - 1 ;
    }
    for ( int j = 0 ; j <= n ; j ++) esp_aux[ n ][ j ] = 0;
    for ( int j = 1 ; j <= n ; j ++)
        for ( int k = 1 ; k<= j ; k ++)
            esp_aux[ 1 ][ j ] += theta_exp_aux[ k ];//la suma de los primeros
    for ( int i = 2 ; i <= n ; i ++)
        for ( int j = i ; j <= n ; j ++)
            esp_aux[i][j]=esp_aux[i][j-1]+esp_aux[i-1][j-1]*theta_exp_aux[j];//theta va de 0..n-1 y esp va de 1..n
    esp[ 0 ] = 1 ;
    
    for ( int i = 1 ; i < n + 1 ; i ++ ) esp[ i ] = esp_aux[ i ][ n ];
}

void Generic::split_elementary_symmetric_polynomial(long double *esp, double *theta,int n, long double **esp_no_a, long double **esp_yes_a){
    for (int k = 0 ; k <= n ; k ++){
        for (int i = 0 ; i < n ; i ++){
            esp_no_a [ k ][ i ] = 0;
            esp_yes_a[ k ][ i ] = 0;
        }
    }
    for (int i = 0; i < n ; i ++){
        esp_no_a[0][i] = 1;
        esp_yes_a[0][i] = 1;//default
        esp_yes_a[1][i] = (exp( theta[ i ]) - 1 );
    }
    for (int k = 1 ; k < n ; k ++)
        for (int i = 0 ; i < n; i++) {
            esp_no_a[k][i] = esp[k] - esp_yes_a[k][i];
            esp_yes_a[ k + 1 ][ i ] =  esp_no_a[k][i] * (exp( theta[ i ]) - 1 );
        }
    for (int i = 0; i < n ; i ++)
        esp_no_a[ n ][ i ] = esp[ n ] - esp_yes_a[ n ][i];
}

void Generic::freq_matrix(int **samples, int m, int n, int **freq){
    for (int i = 0 ; i < n ; i ++)
        for (int j = 0 ; j < n ; j ++)
            freq[ i ][ j ] = 0;
    for (int s = 0;  s < m; s++)
        for (int i = 0;  i < n; i ++)
            freq[ i ][ samples[ s ][i] - 1] ++;
}

void Generic::generate_random_permutation(int len, int first_item_in_perm, int*sigma){
    //implements knuth shuffle (Fisher–Yates shuffle)
    //other option is to code Feller coupling
    for(int i=0;i<len;i++) sigma[i]=i+first_item_in_perm;
    for(int i=0;i<len-1;i++){
        //int pos = rand() % (len-i) + i;
        int pos = (int) (unif_rand() * (len-i) + i);
        int aux= sigma[i];
        sigma[i]=sigma[pos];
        sigma[pos]=aux;
    }
}
void Generic::random_shuffle(int len, int*array){
    //implements knuth shuffle (Fisher–Yates shuffle)
    for(int i=0;i<len-1;i++){
        //int pos = rand() % (len-i) + i;
        int pos = (int) (unif_rand() * (len-i) + i);
        int aux= array[i];
        array[i]=array[pos];
        array[pos]=aux;
    }
}
void Generic::compose(int n, int*s1, int*s2, int*res){
    for(int i = 0 ; i < n ; i++) res[ i ] = s1 [ s2 [ i ] - FIRST_ITEM ];
}

void Generic::compose_sample_right(int **samples, int* sigma, int m, int n, int **composed){
    for (int s = 0 ; s < m ; s ++){
        composed[ s ] = new int[ n ];
        for(int i = 0 ; i < n ; i++) composed[ s ][ i ] = samples[ s ] [ sigma [ i ] - 1 ];
    }
}

void Generic::compose_sample_right(int **samples, int* sigma, int m, int n){
    int * aux = new int[ n ];
    for (int s = 0 ; s < m ; s ++){
        for(int i = 0 ; i < n ; i++) aux[ i ] = samples[ s ] [ sigma [ i ] - 1 ];
        for(int i = 0 ; i < n ; i++) samples[ s ] [ i ] =  aux[ i ];
    }
    delete [] aux;
}

void Generic::invert_sample(int n, int m, int **sample, int **sample_inv){
    //sample_inv must be init as sample_inv = new int*[m];
    for (int s = 0 ; s < m ; s++ ){
        sample_inv[ s ] = new int[ n ];
        invert(n, sample[ s ], sample_inv[ s ]);
    }
}

void Generic::invert(int n, int*sigma, int*res){
    for(int i = 0 ; i < n ; i ++) res[ sigma[ i ] - FIRST_ITEM ] = i + FIRST_ITEM;
    //for(int j = 0 ; j < n_ ; j ++) samples_inv[i][ samples[ i ][ j ] - 1 ] = j + 1;
}
void Generic::get_permu_matrix(int n,int*sigma, int**matrix){
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) matrix[i][j]=0;
    for(int i=0;i<n;i++) matrix[i][sigma[i] - FIRST_ITEM ] = 1;
}

long double Generic::factorial(int val) {
    if(val <= 0) return 1;
    //long  N, b, c, p; // use int for fast calculation and small range of calculation..
    long   b, c;
    long double p, N;
    N=(long double)val;
    c = (long)N - 1;
    p = 1;
    while (c > 0) {
        p = 0;
        b = c;
        while (b > 0) {
            if (b & 1) {
                p += N; // p = p + N;
            }
            // if you would like to use double choose the alternative forms instead shifts
            // the code is fast even!
            // you can use the same tips on double or 64 bit int etc.... but you must... ;-)
            //b >>= 1; // b/=2; (b = b / 2;) ( b >> 1; a.s.r. is more efficent for int or long..!)
            b/=2;
            //N <<= 1; // N += N; N = N + N; N = N * 2; (N <<=1; a.s.l. is more efficent for int or long..!)
            N += N;
        } // end of: while(b>0)
        N = p;
        c--; // c = c - 1;
    } // end of: while(c > 0)
    //printf("[%d] is the factorial! \n", p);
    return p;
}
void Generic::init_factorials (int n) {
    if (facts_ == NULL) {
        facts_n_ = n;
        facts_=new long double[n+1];
        facts_[0]=1;
        for(int i=1;i<=n;i++) facts_[i] = facts_[i-1] * i;
    }
}
long double Generic::count_permus_with_at_least_k_unfixed_points(int n, int k){
    if ( facts_ == NULL ) {
        init_factorials(n);}
    //else if (facts_n_ < n ) {cout<<"Check n in Generic::count_permus_no_fixed_points. ";exit(0);}
    long double  sum = 0 , aux = 0 ;
    int     multi = -1 ;
    for(int i = 1 ; i <= k ; i++){
        //num = num + (-1)^j * factorial(n-l) * factorial(l) / (factorial(j) * factorial(l-j) );
        aux = (long double) multi*(long double) facts_[k]*facts_[n - i]/(facts_[i]*facts_[k-i]);
        sum += aux;
        multi *= -1;
    }
    return facts_[ n ] + sum;
}
double Generic::count_perm_fixed_points(int k, int j){
    if(j<0 || j>k) return 0;
    if(k == 0 && j == 0) return 1;
    return count_perm_fixed_points(k-1, j-1) + count_perm_fixed_points(k - 1, j)*(k - 1  - j)
    + count_perm_fixed_points(k - 1, j + 1)*(j + 1);
    //(1) For every j < 0 or j > k : f(k, j) = 0.
    //(2) f(0, 0) = 1.
    //(3) For every k > 1 and k ≥ j ≥ 0, f(k, j) = f(k − 1, j − 1) + f(k − 1, j)·(k − 1  − j) + f(k − 1, j + 1)·(j + 1)
}

bool Generic::valid_permutation(int*sigma, int n){
    bool*per=new bool[n];
    for(int i=0;i<n;i++)per[i] = false;
    for(int i=0;i<n;i++)
        if(sigma[i] > 0 && sigma[i] <= n && !per[sigma[i]-1] )
            per[sigma[i]-1] = true ;
        else return false;
    delete[]per;
    return true;
}

/*void Generic::riffle_shuffle(int n, int m, int times, int *sigma, int **shuffle){
    // implements gilbert shannon reeds shuffle
    int l,  l_min = 0,  l_max;
    int r,  r_min,      r_max = n;
    int cont = 0 ,      ran_int = 0;
    int center, ini , end;
    int *aux = new int[ n ];
    double * acumul = new double[ n + 1 ];
    double ran ;init_factorials(n);
    acumul [ 0 ] = 1;
    for ( int i = 1 ; i < n +1 ; i++) acumul[ i ] = acumul[ i - 1 ] +  facts_[ n ] / ( facts_[ i ] * facts_[ n - i ]);
    
    
    for (int i = 0 ; i < m ; i++){
        shuffle[ i ] = new int[ n ];
        for (int j = 0 ; j < n ; j++) aux [ j ] = sigma[ j ];
        for (int cont_times = 0 ; cont_times < times ; cont_times++) {
            l_min = 0; r_max = n ; cont = 0 ;
            //ran = (double) rand() / (RAND_MAX ) ;
            
            ran *= pow(2, (double)n);
            center =  0;
            while (acumul[ center ] < ran) center++;
            //cout<<" . "<<center;
            l_max = center;
            r_min = center ;
            l = l_max - l_min;
            r = n - l ;
            while (l > 0 && r > 0){
                ran_int = rand() %(l + r );
                if (ran_int < l) {
                    shuffle[ i ][ cont ] = aux[ l_min ];
                    l_min ++;
                    l--;
                }else{
                    shuffle[ i ][ cont ] = aux[ r_min ];
                    r_min ++;
                    r--;
                }
                cont++;
            }
            if (r == 0 ) {
                ini = l_min;
                end = l_max;
            }else {
                ini = r_min;
                end = r_max;
            }
            for (int j = ini ; j < end ; j++) shuffle[ i ][ cont ++ ] = aux[ j ];
            //if (cont_times < times - 1 )
                for (int j = 0 ; j < n ; j++) aux [ j ] = shuffle[ i ][ j ];
            //print_int_vector(shuffle , n);
        }
    }
    delete [] aux;    
}

void Generic::seed(void) {
    int fd, buf;
    if ((fd = open("/dev/urandom", O_RDONLY)) < 0) {
        perror("/dev/urandom");
    }
    read(fd, &buf, sizeof (buf));
    close(fd);
    srand(buf);
}*/

void Generic::partition_function_init(int n){
    //F(n,k)=F(n,k−1)+F(n−k,k)
    partition_table=new int*[n+1];
    for(int i=0; i< n+1; i++)
        partition_table[i]=new int[n+1];
    for(int i=0; i< n+1; i++) partition_table[0][i]=1;
    for(int i=0; i< n+1; i++) partition_table[i][0]=0;
    for(int i=1; i< n+1; i++){
        for (int j = 1 ; j < n+1 ; j ++) {
            if (i - j < 0 ) partition_table [i][j] = partition_table[i][j-1];
            else partition_table [i][j] = partition_table[i][j-1] + partition_table[i-j][j];
        }
    }
}

void Generic::insert_at(int *sigma, int n, int move, int to, int*res){
    if ( move < to ){
        for (int i = 0 ; i < move  ; i ++) res[i] = sigma[i];
        for (int i = move ; i < to  ; i ++) res[i] = sigma[i+1];
        res[to] = sigma[move];
        for (int i = to+1 ; i < n  ; i ++) res[i] = sigma[i];
    }else{
        for (int i = 0 ; i < to ; i ++) res[i] = sigma[i];
        res[to] = sigma[move];
        for (int i = to+1 ; i <= move ; i ++) res[i] = sigma[i-1];
        for (int i = move+1 ; i < n  ; i ++) res[i] = sigma[i];
    }
}

int Generic::get_number_of_partitions(int n){
    return partition_table[n][n];
}

void Generic::partition_function_destructor(int n){
    for(int i=0; i< n+1; i++) delete [] partition_table[i];
    delete [] partition_table;
}
double Generic::milliseconds(){////returns milliseconds
    // cpu time in seconds since start of run.
    /*double secs;
     secs = (double)(clock()/ (double) CLOCKS_PER_SEC);//  / 1000.0
     return(secs);*/
    double msecs;
    msecs = 1000*(clock()/ (double) CLOCKS_PER_SEC);
    return(msecs);
}

void Generic::ini_chrono(){chrono = milliseconds();}

double Generic::end_chrono(){
    return  milliseconds() - chrono;
}
