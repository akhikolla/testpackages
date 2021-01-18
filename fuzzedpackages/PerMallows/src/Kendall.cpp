//
//  Kendall.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//
#include <R.h>

#include <cmath>
#include "Newton_raphson.h"
#include "Kendall.h"

double Kendall::probability(int *s, int *s_0, double *theta){
    int*v = new int[n_],*comp = new int[n_], *inv = new int[n_];
    for(int i=0;i < n_;i++) inv[s_0[i]-1]=i+1;
    double proba=1;
    for(int i=0;i<n_;i++) comp[i]=s[inv[i]-1];
    double*psi=new double[n_];
    perm2dist_decomp_vector(comp, v );
    calculate_psi(theta, psi);
    for(int i=0;i<n_-1;i++) proba *= exp(-theta[i]*v[i])/psi[i];
    delete[]psi;
    delete[]comp;
    delete[]inv;
    delete[]v;
    return proba;
}

double Kendall::calculate_psi(double*theta, double*psi_vector){
    double psi = 1;
    for(int i=0;i< n_ - 1;i++){
        int j=i+1;
        psi_vector[i]= (1 - exp(-theta[ i ]*(n_ - j + 1 ))) / (1 - exp(-theta[ i ]));
        psi *= psi_vector[ i ];
    }
    return psi;
}

int Kendall::perm2dist_decomp_vector(int*sigma, int*vec ){
    int dist = 0 ;
    if ( vec != NULL ) for(int i = 0 ; i < n_ - 1 ; i++)vec[ i ] = 0;
    for (int i = n_ - 2 ; i >= 0; i--)
        for (int j = i + 1 ; j < n_ ; j++)
            if(sigma[ i ] > sigma [ j ]){
                if ( vec != NULL ) vec[ i ]++;
                dist ++;
            }
    return dist;
}

int Kendall::distance(int*s1, int*s2){
    int*comp=new int[n_], *inv = new int[ n_];
    for(int i = 0 ; i < n_ ; i ++) inv[ s2[ i ] - 1 ] = i + 1;
    for (int i=0; i<n_ ; i++)   comp[i]=s1[inv[i]-1];
    int dist = perm2dist_decomp_vector(comp, NULL); 
    delete [] comp;
    delete [] inv;
    return dist;
}


long double Kendall::num_permus_at_distance(int d){
    if( d > n_*(n_-1)/2) return 0;
    return count_[ n_ ][ d ];
}

void Kendall::random_permu_at_dist_d( int dist, int*sigma  ){
    //it is done in O(max_dist)=O(n^2/2), so optimal ???
    //it starts by sigma[n-2]
    // for each i (decreasing) randomly chooses a position j>i.
    // Thus, there is no need to initialise sigma
    double* acum = new double[ n_ ];
    int*v = new int[ n_ ];
    v[ n_ - 1 ] = 0 ;
    int i;
    for(i = 0 ; ( i < n_ && dist > 0 );i ++ ) {
        int rest_max_dist = (n_ - i - 1 ) * ( n_ - i - 2 ) / 2;//con los restantes n' puedes tener distMAx de binom(n`)(2)
        if(rest_max_dist  >= dist )acum[ 0 ] = count_[ n_ - i - 1 ][ dist ];
        else acum[ 0 ] = 0 ;
        int max = (n_ - i < dist + 1 ) ? (n_ - i ) : dist + 1;
        for(int j = 1; j < max; j++ )
            if(rest_max_dist + j >= dist) acum[j] = acum[j-1] + count_[ n_ - i - 1 ] [ dist - j];
            else acum[ j ] = 0;
        //double bound = (double)rand() / (double)(RAND_MAX) * acum[ max - 1 ];
        double bound = unif_rand() * acum[ max - 1 ];
        int pos = 0 ;
        while(acum[pos] <= bound) pos++;
        //if(!mute)cout<<". choose "<<pos<<endl;
        dist -= pos;
        v[i] = pos;
    }
    for (int j = i; j < n_; j++) v[ j ] = 0; //the last n-i positions
    dist_decomp_vector2perm(v, sigma);
    
    delete [] v;
    delete [] acum;
}

void Kendall::random_sample_at_dist(int dist, int m, int **samples){
    for (int i = 0 ; i < m ; i ++){
        samples[i] = new int[n_];
        random_permu_at_dist_d(dist, samples[ i ]);
    }
}

void Kendall::dist_decomp_vector2perm(int* vec, int* sigma) {
    int val;
    int*ident = new int[n_];
    for(int i = 0 ; i < n_ ; i ++) ident[i]=i;
    for(int i = 0 ; i < n_ - 1 ; i++){
        val = vec[i];
        int index = 0;
        while( !(ident[ index ] != -1 && val == 0))
            if(ident[ index ++ ] != -1)
                val --;
        sigma[ i ] = index + FIRST_ITEM ;
        ident[ index ] = -1 ;
    }
    int index=0;
    while(ident[ index ] == -1 )index++;
    sigma[ n_ - 1 ] = index + FIRST_ITEM;
    delete [] ident;
}

void Kendall::distances_sampling(int m, double theta, int**samples) {
    int target_dist=0;
    int d_max = n_*(n_-1)/2 ;
    long double * acumul = new long double[ d_max + 1];//+1
    double rand_val;
    acumul[ 0 ] = exp( -theta  * 0 ) * count_[ n_ ][ 0 ];
    for(int dista = 1 ; dista < d_max + 1 ; dista++)
        acumul[ dista ] = acumul[ dista - 1 ] +  exp(-theta  * dista) * count_[ n_ ][ dista ];
    for( int i = 0 ; i < m ; i++ ){
        target_dist = 0;
        //rand_val = (double) acumul[ d_max ] * (double) rand() / RAND_MAX;
        rand_val = (double) acumul[ d_max ] * unif_rand();
        while ( acumul[ target_dist ] <= rand_val ) target_dist ++;
        int *sigma = new int[ n_ ];
        random_permu_at_dist_d( target_dist , sigma);
        samples[ i ] = sigma;
    }
    delete [] acumul;
}


void Kendall::multistage_sampling(int m, double*theta, int**samples){
    double *psi = new double[ n_ - 1 ];
    int *v = new int[ n_ ];
    long double **vprobs = new long double*[ n_];
    for ( int i = 0  ; i < n_ ; i ++) vprobs[i] = new long double [n_ ];
    
    for(int i = 0 ; i < n_ - 1 ; i ++)
        psi[i] = (1 - exp(( - n_ + i )*(theta[ i ])))/(1 - exp( -theta[i]));

    for(int j = 0 ; j < n_ - 1 ; j++){
        vprobs[j][0] = 1 / psi[j];
        for(int r = 1 ; r < n_ - j ; r ++) 
            vprobs[j][r] = vprobs[j][ r - 1 ] + exp( -theta[ j ] * r ) / psi[ j ];
        
    }
    for(int samp = 0 ; samp < m ; samp++){
        for( int i = 0 ; i < n_ - 1 ; i ++ ){
            int target_v = 0 ;
            //double rand_val = (double) vprobs[i][ n_ - i - 1 ] * (double) rand() / RAND_MAX;
            double rand_val = (double) vprobs[i][ n_ - i - 1 ] * unif_rand();
            while(vprobs[i][target_v] <= rand_val) target_v ++;
            v[i] = target_v;
        }        
        v[ n_ - 1 ] = 0;
        int *sigma=new int[ n_ ];
        dist_decomp_vector2perm(v, sigma);
        samples[samp]=sigma;
    }
    delete [] v;
    delete [] psi;
    for ( int i = 0  ; i < n_ ; i ++) delete [] vprobs[i];
    delete [] vprobs;
}


void Kendall::gibbs_sampling(int m, double *theta, int model, int **samples) {
    Generic     gen;
    int*sigma   = new int[ n_ ];
    int *v      = new int[ n_ ];
    int burning_period_samples = n_*log(n_);
    gen.generate_random_permutation( n_ , 1, sigma);
    perm2dist_decomp_vector( sigma , v);
    
    for(int sample = 0 ; sample < m + burning_period_samples ; sample ++){
        int i;
        //i = rand() % (n_ - 1);
        i = (int) ( unif_rand() * (n_ - 1) );
        int a = sigma[ i ];
        int b = sigma[ i + 1 ];
        bool make_swap = false;
        if( a > b )  make_swap = true;
        else{
            //double rand_double = (double)rand()/RAND_MAX;
            double rand_double = unif_rand();
            if( model == MALLOWS_MODEL ){
                if(rand_double < exp(-theta[0])) make_swap = true;
            }else{
                //\[	exp(-\theta_i*v_{i+1}(\sigma) - \theta_{i+1}*v_i(\sigma) +\theta_i * v_i(\sigma) + \theta_{i+1}* v_{i+1}(\sigma)) 	\]
                //equiv TODO
                double ratio = exp(-theta[ i ] * (v[ i + 1 ]+1)  -theta[ i + 1 ] * v [ i ]   +theta[ i ] * v[ i ] +theta[ i + 1 ] * v [ i + 1 ]);
                if(rand_double < ratio ) make_swap = true;
            }
        }
        if(make_swap){
            sigma[ i ] = b;
            sigma[ i + 1 ] = a;
            int aux = v[ i ];
            v[ i ] = v[ i + 1 ];
            v[ i + 1 ] = aux;
            if (a < b )  v [ i ] ++;
            else  v[ i + 1 ] --;
        }
        if(sample >= burning_period_samples){
            samples[sample-burning_period_samples] = new int[ n_ ];
            for(int i = 0  ; i < n_ ; i ++)   samples[ sample - burning_period_samples ][ i ] = sigma[ i ];
        }
    }
    delete [] sigma;
    delete [] v;
}

void Kendall::borda(int **samples, int m, int * sigma_0){
    Generic gen;
    double*average = new double[n_];
    for (int j = 0 ; j < n_ ; j++) average[j] = 0;
    int * sigma_0_inv = new int [n_];
    for (int s = 0 ; s < m ; s ++)
        for (int j = 0 ; j < n_ ; j++)
            average [ j ] += samples[ s ][ j ];
	//gen.print_double_vector(average, n_);
    for(int i = 0 ; i < n_ ; i ++){
        int min_pos = -1, cur_pos = -1 , min = -1;
        while (++cur_pos < n_ )
            if ( ( average[ cur_pos ] != -1 && average[ cur_pos ] < min ) || min == -1 ) {
                min = average[ cur_pos ];
                min_pos = cur_pos;
            }
        sigma_0[ min_pos ] = i + 1;
        average[ min_pos ] = -1;
    }
    for(int i = 0 ; i < n_ ; i++)sigma_0_inv[sigma_0[ i ] - 1 ] = i + 1;
    
    delete [] sigma_0_inv;
    delete [] average;
    //return lis;
}

long double Kendall::get_likelihood(int m, int** samples, int model, int * sigma_0){
    Newton_raphson newton(n_);
    long double loglikelihood = 0;
    double  * theta = new double[ n_ ];
    double  * psi_vec = new double [ n_ ];
    double  psi;
    if(model == MALLOWS_MODEL){
        int dist = distance_to_sample(samples, m, sigma_0);
        theta[0] = newton.Newton_raphson_method((double)dist/m, -1.001 , KENDALL_DISTANCE, MALLOWS_MODEL, -1, NULL);
        for (int i = 1 ; i < n_ -1; i++) theta[ i ] = theta[ 0 ];
        psi = calculate_psi(theta, psi_vec);
        loglikelihood = - theta[0] * dist - m * log (psi);
    }else{
        int * s_inv = new int [n_], *comp = new int[ n_ ], *v = new int[ n_ ], *v_avg = new int[ n_ ];
        for( int j = 0 ; j < n_ -1 ; j++) v_avg[ j ] = 0;
        for( int i = 0 ; i < n_ ; i ++) s_inv[ sigma_0[ i ] - 1 ] = i + 1;
        for( int i = 0 ; i < m ; i ++){
            for(int j = 0 ; j < n_ ; j++) comp[ j ] = samples[i] [ s_inv [ j ] - 1 ];
            perm2dist_decomp_vector(comp, v);
            for( int j = 0 ; j < n_ -1 ; j++) v_avg [j] += v[j];
        }
        for( int j = 0 ; j < n_ -1 ; j++)
            theta[ j ] = newton.Newton_raphson_method( (double)v_avg[j]/m , -1.001 , KENDALL_DISTANCE, model, j+1 , NULL);
        psi = calculate_psi(theta, psi_vec);
        for( int j = 0 ; j < n_ - 1 ; j++)
            loglikelihood += theta[ j ] * ((double)v_avg[ j ]/m) + log (psi_vec[ j ]);
        loglikelihood *= -m;
        delete [] comp;
        delete [] v;
        delete [] s_inv;
    }
    delete [] psi_vec;
    delete [] theta;
    return loglikelihood;
}

void Kendall::estimate_theta(int m, int*sigma_0, int**samples, int model, double*theta){
    Newton_raphson newton(n_);
    double *psi_vec = new double [ n_ ];
    if(model == MALLOWS_MODEL){
        int dist = distance_to_sample(samples, m, sigma_0);
        theta[0] = newton.Newton_raphson_method((double)dist/m, -10.001 , KENDALL_DISTANCE, MALLOWS_MODEL, -1, NULL);
    }else{
        int * s_inv = new int [n_], *comp = new int[ n_ ], *v = new int[ n_ ], *v_avg = new int[ n_ ];
        for( int j = 0 ; j < n_ -1 ; j++) v_avg[ j ] = 0;
        for( int i = 0 ; i < n_ ; i ++) s_inv[ sigma_0[ i ] - 1 ] = i + 1;
        for( int i = 0 ; i < m ; i ++){
            for(int j = 0 ; j < n_ ; j++) comp[ j ] = samples[i] [ s_inv [ j ] - 1 ];
            perm2dist_decomp_vector(comp, v);
            for( int j = 0 ; j < n_ -1 ; j++) v_avg [j] += v[j];
        }
        for( int j = 0 ; j < n_ -1 ; j++)
            theta[ j ] = newton.Newton_raphson_method( (double)v_avg[j]/m , -10.001 , KENDALL_DISTANCE, model, j+1 , NULL);
        delete [] comp;
        delete [] v;
        delete [] s_inv;
        delete [] v_avg;
    }
    delete [] psi_vec;
}


int Kendall::distance_to_sample(int **samples, int m,  int *sigma){
    int distance= 0;
    int *comp = new int[ n_ ], *sigma_inv = new int[ n_ ];
    for(int j = 0 ; j < n_ ; j ++) sigma_inv[sigma[ j ] - 1 ] = j + 1;
    for(int s = 0 ; s < m ; s ++){
        for(int i = 0 ; i < n_ ; i ++) comp[ i ] = samples[ s ][ sigma_inv [ i ] - 1 ];
        distance += perm2dist_decomp_vector(comp, NULL);//get_distance_and_X(comp, NULL);
    }
    delete []sigma_inv;
    delete []comp;
    return distance ;
}

double Kendall::expectation(double theta){
    double res =  0, aux = 0 ;
    for (int j = 1 ; j <= n_ ; j ++){
        aux =  exp(-j * theta);
        res += j * aux / (1 - aux);
    }
    aux = exp(-theta);
    return n_*aux /(1-aux) - res;
}

void Kendall::expectation(double *theta, double *expect){
    double aux1 = 0 , aux2=0;
    int j = 0 ;
    for (int i = 0 ; i < n_ - 1 ; i ++){
        j = i + 1 ;
        aux1 = exp(-theta[ i ] * (n_ - j  + 1 ));
        aux2 = exp(-theta[ i ]);
        expect[ i ] = aux2 / (1 - aux2 ) - (n_ - j + 1 ) *aux1 / (1 - aux1);
    }
    expect[ n_ - 1 ] = 0;
}

