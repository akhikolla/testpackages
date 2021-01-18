//
//  Hamming.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include "Hamming.h"
#include "Generic.h"
#include "Lap.h"
#include <cmath>
#include <float.h>
#include "Newton_raphson.h"
#include <R.h>


double Hamming::probability(int *s, int *s_0, double *theta){
    double  pro = 0;
    bool    MM = true ;
    for (int i = 0 ; (i < n_ -1 && MM); i++) if ( theta[i] != theta[i+1]) MM = false;
    if(MM){
        int sum = 0;
        for ( int i = 0 ; i < n_ ; i ++) if(s[i] != s_0[i]) sum ++;
        return exp( -sum*theta[0] )/psi_hm(theta[0]);
    }else{
        for (int i = 0 ; i < n_ ; i++)
            if (s[i] == s_0[i]) pro += 1;
            else pro +=  theta[i] ;
        return exp(-pro )/(double)psi_whm(theta);
    }
}

int Hamming::distance(int*s1, int*s2){
    int     sum = 0;
    for ( int i = 0 ; i < n_ ; i ++) if(s1[i] != s2[ i ]) sum ++; //if(sigma[i] != i + 1) sum ++;
    return sum;
}

void Hamming::random_derangement(int n, int *sigma){

    if( n== 2) {
        sigma[0]=2;
        sigma[1]=1;
    //}else if ( (n-1)*deran_num_[n-1] / deran_num_[n] > (double) rand() / RAND_MAX ){
    }else if ( (n-1)*deran_num_[n-1] / deran_num_[n] > unif_rand() ){
        random_derangement(n-1, sigma);
        //int ran = rand() % (n - 1 );
        int ran = (int) (unif_rand() * (n - 1 ));
        sigma[n-1] = sigma[ran];
        sigma[ran]= n;
    }else{
        int*deran = new int[ n - 2 ], *conv = new int[n-1];
        random_derangement( n - 2 , deran );
        //int ran = rand() % (n - 1 );
        int ran =(int)( unif_rand() * (n - 1 ));
        int j = 0;
        for (int i= 0; i < n-1 ; i ++)
            if ( i != ran ) conv[j++] = i+1;
        j = 0;
        for (int i = 0; i < n-1 ; i ++)
            if (i != ran )
                sigma[i] = conv[deran[j++]-1];
        sigma[ran]= n ;
        sigma[n-1] = ran+1;
        delete [] deran;
        delete [] conv;
    }
}

void Hamming::random_sample_at_dist(int dist, int m, int **samples){
    for (int i = 0 ; i < m ; i++) {
        samples[i] = new int[n_];
        random_permu_at_dist_d(dist, samples[i]);
    }
}

void Hamming::random_permu_at_dist_d(int dist, int *sigma){
    Generic gen;
    int     *ran = new int[n_];
    gen.generate_random_permutation(n_, 1, ran );
    generate_permu_from_list(ran, dist, sigma);
    delete [] ran;
}

int Hamming::perm2dist_decomp_vector(int*sigma, int*vec ){
    int dist = 0 ;
    for (int i = 0 ; i < n_ ;i++)
        if (sigma[ i ] != i + 1 ){
            dist ++;
            vec[ i ] = 1;
        }else
            vec[ i ] = 0;
    return dist;
}

void Hamming:: dist_decomp_vector2perm(int* vec, int* sigma){
    int     *ran = new int[n_];
    int     last = n_ - 1, first = 0;
    for (int i = 0 ; i < n_; i++) {
        if (vec[i] == 0 ) ran[ last-- ] = i + 1 ;
        else ran[first ++] = i + 1;
    }
    //if ( first == 1 )bool traza = true;
    generate_permu_from_list(ran, first, sigma);
    //Generic gen;cout<< "h ";gen.print_int_vector(h, n_); cout<< "p ";gen.print_int_vector(permu, n_);
    delete [] ran;
}

void Hamming::generate_permu_from_list(int*ran, int dist, int*sigma){
    //the first d items in ran will be deranged. for the rest, sigma[i]=i
    if ( dist == 0) {
        for (int i = 0 ; i < n_ ; i++) sigma[i] = i + 1 ;
        return;
    }
    int     *deran=new int[n_];
    if(dist > 1) random_derangement(dist, deran);
    for( int i = 0; i < dist; i++)
        sigma[ ran[i] - 1 ] = ran  [ deran  [ i ] - 1 ];
    for (int i = dist ; i < n_; i++)
        sigma[ran[i] - 1 ]= ran[i];
    delete [] deran;
}

void Hamming::distances_sampling(int m, double theta, int **samples){
    // limit: n = 90
    Generic     gen;
    int         d_max = n_ ;
    int         target_dist=0;
    long double      rand_val;
    long double      * fact       = new long double[ n_ + 1 ];
    long double      * hamm_count = new long double[ n_ + 1 ];
    long double      * acumul     = new long double[ d_max + 1];//+1
    
    fact[0]=1; fact[1]=1;
    for ( int i = 2 ; i <= d_max ; i++) fact[i] = i * fact[i-1];
    for ( int d = 0 ; d <= d_max ; d++) hamm_count[d] = (long double) deran_num_[d] *fact[n_] / (fact[d] * fact[n_ - d]) ;
    acumul[ 0 ] = exp( -theta  * 0 ) * hamm_count[ 0 ];
    for ( int dista = 1 ; dista <= d_max  ; dista++)
        acumul[ dista ] = acumul[ dista - 1 ] +  (long double) exp(-theta  * dista) * hamm_count[ dista ];
    
    //gen.print_double_vector(hamm_count, n_+1);
    //for (int d = 0 ; d <= d_max ; d++) cout<<acumul[d]<<" ";cout<<endl;
    for( int i = 0 ; i < m ; i++ ){
        target_dist = 0;
        //rand_val = (long double) acumul[ d_max ] * (double) rand() / RAND_MAX;
        rand_val = (long double) acumul[ d_max ] * unif_rand();
        while ( acumul[ target_dist ] <= rand_val ) target_dist ++;
        int *sigma = new int[ n_ ];
        random_permu_at_dist_d( target_dist , sigma);
        samples[ i ] = sigma;
    }
    delete [] fact;
    delete [] acumul;
    delete [] hamm_count;
}

double Hamming::psi_whm_reverse(double*theta){
    //oposite meaning of h_j (como el de Fligner)
    Generic     gen;
    long double res = 0 ;
    long double * esp    = new long double [ n_ + 1 ];
    double      * theta_inv = new double [ n_ ];
    for ( int k = 0 ; k < n_ ; k ++) theta_inv[ k ] = -theta[ k ];
    gen.elementary_symmetric_polynomial( theta_inv, n_, t_, aux_esp_, esp );
    for ( int k = 0 ; k <= n_ ; k ++)
        res += g_n_[ n_ ][k] * esp[k];
    delete [] esp;
    delete [] theta_inv;
    return res;
}

double Hamming::psi_whm(double*theta){
    Generic     gen;
    long double res = 0 , sum_theta = 0;
    long double *esp = new long double [ n_ + 1 ];
    for ( int k = 0 ; k < n_ ; k ++) sum_theta += theta[ k ];
    gen.elementary_symmetric_polynomial( theta, n_, t_, aux_esp_,esp );
    for ( int k = 0 ; k <= n_ ; k ++)
        res += facts_[ n_ - k ] * esp[k];
    
    delete [] esp;
    return (res * exp( -sum_theta ));
}

double Hamming::psi_hm_reverse(double theta){ // como el de fligner pero las epsilon son lo contrario
    long double  sum = 0, aux;
    Generic gen;
    //gen.init_factorials( n_ );
    for ( int i = 0 ; i <= n_ ; i ++ ){
        //aux = (gen.count_permus_with_at_least_k_unfixed_points(n_,i) / (long double)(gen.get_factorial_in_table(n_-i)))/gen.get_factorial_in_table(i);
        aux = (g_n_[ n_ ][ i ] / (facts_[ n_ - i ]))/facts_[ i ];
        sum += (long double)aux*(long double)pow((exp(-theta)-1),i);
    }
    return (double)facts_[ n_ ] * sum;
}

double Hamming::psi_hm(double theta){ //fligner's
    double  sum=0;
    Generic gen;
    for ( int i = 0 ; i <= n_ ; i ++ )
        sum += ((double)pow((exp(theta)-1),i)/(double)gen.factorial(i));
    sum = sum * exp( - n_ * theta) * gen.factorial( n_ );
    return sum;
}

/*
long double Hamming::compute_marginal_slow(int *h , double * theta, int  marginal_order){
    //|A| set of fixed points; |B| set of unfixed points
    Generic gen;
    int     a = 0 , b = 0 ;//a : |A|; b = |B|
    long double  res = 0;
    long double * esp_red_slow = new long double[ n_+1];
     double * theta_red_slow = new  double[ n_+1];
    double theta_acum_B = 0 , theta_acum_not_in_AB = 0;
    for (int i = 0 ; i < marginal_order ; i++ ){
        if ( h[ i ] == 1 ){
            theta_acum_B += theta[ i ];
            b++;
        }
    }
    a = marginal_order - b ;
    for (int i = marginal_order ; i < n_ ; i++ ){
        theta_red_slow[ i - marginal_order ] = theta[ i ];
        theta_acum_not_in_AB += theta[ i ];
    }
    gen.elementary_symmetric_polynomial( theta_red_slow, n_ - marginal_order ,t_, aux_esp_, esp_red_slow );
//    for (int i = 0 ; i <= n_ - marginal_order ; i ++ )cout<<esp_red_slow[ i ]<<" "; cout<<" slow"<<endl;//cout<<"a, b "<<a<<" "<<b<<endl;
    for (int k = 0 ; k <= n_ - marginal_order ; k++ )
        res +=  gen.count_permus_with_at_least_k_unfixed_points(n_ -a -k, b) * (esp_red_slow[k] );
    //cout<<"res s "<<res<<" "<<- theta_acum_not_in_AB - theta_acum_B<<endl;
    //cout<<" theta "<<- theta_acum_not_in_AB - theta_acum_B<<endl;
    return exp( - theta_acum_not_in_AB - theta_acum_B )* res;
}**/

long double Hamming::compute_marginal_iterative(int *h , double * theta, int  marginal_order){
    // must be initialized :
    //esp_ini: elementary_symmetric_polynomial( theta, n_ ,..., esp_ini_ );
    //t_sampling_[i]=exp(theta[i]-1
    int     a = 0 ;//a : |A|; b = |B| is global
    long double  result = 0;
    int     current_var = marginal_order - 1;
    int     num_vars = n_ - marginal_order;
    long double  res = 0;
    

    if (marginal_order == 1 ){//the first time it is called
        theta_acum_not_in_A = 0;
        b_ = 0 ;
        for (int i = 0 ; i < n_ ; i++ ){
            theta_acum_not_in_A += (long double) theta[ i ];
            esp_red_[ i ] = esp_ini_[ i ];
        }
        esp_red_[ n_ ] = esp_ini_[ n_ ];
    }
    if ( current_var > 0 ){
        if( h[ current_var - 1 ] == 0 )
            theta_acum_not_in_A -= (long double) theta[ current_var - 1];//the set of fixed points is A. If the last position was sampled as Unfix , update set
        else
            b_ ++;//otherwise, (current_var-1) \in B
    }
    a = marginal_order - b_; // ojo: h[current var] = 0 => current_var \in A
    //split the ESP by variable current_var:
    //esp -> esp_no_a + esp_yes_a
    //since esp in the next iteration we be esp_no_a of the current iter (esp=esp_no_a) then
    //esp_no_a is directly stored in esp
    esp_red_yes_a_[1]= t_sampling_[ current_var ];
    for (int k = 1 ; k < num_vars ; k ++){
        esp_red_[k]= esp_red_[k] - esp_red_yes_a_[k];
        esp_red_yes_a_[ k + 1 ] =  esp_red_[k]* t_sampling_[ current_var ];
        res +=  g_n_[ n_ -a -k ][ b_ ] * esp_red_[k];
    }
    esp_red_[num_vars]= esp_red_[num_vars] - esp_red_yes_a_[num_vars];
    res +=  g_n_[n_ - a ][ b_ ] ;//* esp_red_[0]; // iter k=0
    if (num_vars != 0)
        res +=  g_n_[n_ -a - num_vars][ b_ ] * esp_red_[num_vars]; // iter k= num_vars
    
    result = (long double) exp( - theta_acum_not_in_A + theta[ current_var ] )* res;
    if ( result < 0 ){
        //cout<<"ERROR Negative marginal probability, maybe theta is too large?"<<endl;
        //exit(0);
    }
    return result;
}

long double Hamming::compute_marginal(int *h , double * theta ){
    Generic gen;
    double  * theta_red = new double[ n_ ];
    double  theta_not_in_a = 0;
    int     a = 0 , b = 0 , num_vars = 0;
    long double res =  0, psi = 0 ;

    for (int i = 0 ; i < n_ ; i++ ){
        if( h[ i ] == 0) a++;
        else theta_not_in_a += theta[ i ];
        if( h[ i ] == 1 ) b++;        
        if( h[ i ] != 1 &&  h[ i ] != 0 ) {
            theta_red[ num_vars ++ ] = theta[ i ];
        }
    }
    psi = psi_whm(theta);
    gen.elementary_symmetric_polynomial( theta_red, num_vars , t_, aux_esp_, esp_red_ );
    for (int k = 0 ; k < num_vars + 1 ; k ++)
        res +=  g_n_[ n_ -a -k ][ b ] * esp_red_[k];
    
    delete [] theta_red;
    return (long double) exp( - theta_not_in_a  )* res / psi;
}


void Hamming::multistage_sampling(int m, double *theta, int **samples){
    Generic gen;
    int     *h = new int[ n_ ]; for ( int i = 0 ; i < n_ ; i ++ ) h[ i ] = 0;
    long double  marg = 0, marg_0 = 0, marg_1 = 0, rand_double = 0 ;
    long double  marg_ini = 0;
    
    for ( int i = 0 ; i < n_ ; i ++ ) t_sampling_[ i ] = (long double) exp( theta[ i ]) - 1 ;

    marg_ini = psi_whm(theta);
    
    //initialize esp_ini_;
    gen.elementary_symmetric_polynomial( theta, n_ ,t_, aux_esp_, esp_ini_ );
    for (int s = 0 ; s < m ; s ++) {
        marg = marg_ini;
        for (int items_set = 0 ; items_set < n_; items_set++) {
            h[ items_set ] = 0; //for the marginal computation
            marg_0 = rand_double = 0;
            marg_0 = compute_marginal_iterative(h , theta, items_set + 1);//8
            marg_1 = marg - marg_0;
            //rand_double = marg * (double)rand() / (RAND_MAX);//2
            rand_double = marg * unif_rand();
            if (rand_double < marg_0) {
                marg = marg_0;
                h[ items_set ] = 0;
            }else{
                marg = marg_1;
                h[ items_set ] = 1;
            }
        }
        samples[ s ] = new int[ n_ ];
        dist_decomp_vector2perm(h , samples[ s ]);
    }
    delete [] h;
 
}

void Hamming::gibbs_sampling(int m, double *theta, int model, int **samples){
    int burning_period_samples = n_*n_;
    int*sigma = new int[ n_ ];
    Generic  gen;
    gen.generate_random_permutation( n_ , 1, sigma);
    double ratio = 0 , rand_double = 0 ;
    int h_i = 0, h_j = 0, h_i_new = 0, h_j_new = 0;
    
    for ( int sample = 0 ; sample < m + burning_period_samples ; sample ++){
        int i, j;
        do{
            i = (int) (unif_rand() * n_);//rand() % (n_);
            j = (int) (unif_rand() * n_);//rand() % (n_);
        }while ( i == j );
        h_i     = (i == sigma[i] - 1) ? 0 : 1 ;
        h_j     = (j == sigma[j] - 1) ? 0 : 1 ;
        h_i_new = (i == sigma[j] - 1) ? 0 : 1 ;
        h_j_new = (j == sigma[i] - 1) ? 0 : 1 ;
        
        ratio = exp(- h_j_new * theta[j]) * exp(- h_i_new * theta[i]) / (exp(- h_j * theta[j]) * exp(- h_i * theta[i]));
        //rand_double = (double)rand()/RAND_MAX;
        rand_double = unif_rand();
        if (rand_double < ratio ) {
            int aux = sigma[i];
            sigma[i] = sigma[j];
            sigma[j] = aux;
        }
        if(sample>=burning_period_samples){
            samples[sample-burning_period_samples]=new int[ n_ ];
            for ( int i = 0  ; i < n_ ; i ++)   samples[ sample - burning_period_samples ][ i ] = sigma[ i ];
        }
    }
    delete [] sigma;
}

int Hamming::distance_to_sample(int **samples, int m, int *sigma){
    int dist= 0;
    for(int s = 0 ; s < m ; s ++){
        for ( int i = 0 ; i < n_ ; i ++) if(samples[ s ][ i ] != sigma[ i ]) dist ++;
    }
    return dist;
}

void Hamming::estimate_consensus_exact_mm(int m, int**samples, int*sigma){
    int**freq = new int*[n_];
    Lap lap;
    int*rows=new int[n_], *cols=new int[n_], *u=new int[n_], *v=new int[n_];
    for (int i = 0 ; i < n_ ; i++ ){freq[i]=new int[n_]; for (int j = 0 ; j < n_ ; j++)freq[i][j] = 0 ;}
    
    for (int i = 0 ; i < m ; i ++)
        for (int j = 0 ; j < n_ ; j ++)
            freq[j][samples[i][j]-1]--;
    
    //int cost = -1 * lap.lap(n_, freq, rows, cols, u, v);
    lap.lap(n_, freq, rows, cols, u, v);
    for (int i = 0; i < n_; i++)
        sigma[i] = rows[i] + 1;
    delete [] rows;
    delete [] cols;
    delete [] u;
    delete [] v;
    for (int i = 0 ; i < n_ ; i++) delete [] freq[ i ];
    delete [] freq;
}


double Hamming::expectation(double theta){
    double x_n = 0 , x_n_1 = 0, aux = 0 ;
    for (int k = 0 ; k <= n_ ; k ++){
        aux = pow (exp(theta )-1, k) / facts_[ k ];
        x_n += aux;
        if (k < n_ )
            x_n_1 += aux ;//pow (exp(theta )-1, k) / facts_[ k ];
    }
    return (double )(n_ * x_n - x_n_1 * exp( theta )) / x_n;
}

void Hamming::expectation(double *theta, double*h_expected){
    double numer = 0 , denom = 0 ;
    long double ** esp_no_a_  = new long double*[ n_ + 1 ];
    long double ** esp_yes_a_ = new long double*[ n_ + 1 ];
    for (int k = 0 ; k <= n_ ; k ++){
        esp_no_a_ [ k ]= new long double [ n_ ];
        esp_yes_a_[ k ]= new long double [ n_ ];
        for (int i = 0 ; i < n_ ; i ++){
            esp_no_a_ [ k ][ i ] = 0;
            esp_yes_a_[ k ][ i ] = 0;
        }
    }
    Generic gen;
    gen.elementary_symmetric_polynomial(theta, n_, t_, aux_esp_,esp_red_ );
    gen.split_elementary_symmetric_polynomial (esp_red_, theta , n_, esp_no_a_, esp_yes_a_);

    for ( int i = 0 ; i < n_ ; i++ ){
        numer = denom = 0 ;
        for (int k = 0 ; k <= n_ ; k ++){
            if (k > 0 )numer += facts_[n_ - k ]  * esp_no_a_[ k - 1 ][ i ];
            denom += facts_[n_ - k ] * esp_red_[ k ];
        }
        h_expected[ i ] = 1 - (numer * exp( theta[ i ])) / denom ;
    }
    for (int k = 0 ; k <= n_ ; k ++){
        delete [] esp_no_a_ [ k ];
        delete [] esp_yes_a_[ k ];
    }
    delete [] esp_no_a_;
    delete [] esp_yes_a_;
}

void Hamming::sample_to_h_vector(int **samples, int m, int * sigma, double *h_avg){
    //if sigma != NULL => h( samples sigma^{-1} )
    for (int i = 0 ; i < n_ ; i++) h_avg[ i ] = 0;
    for (int s = 0 ; s < m ; s ++)
        for (int i = 0 ; i < n_ ; i++)
            if ( sigma == NULL ){
                if (samples[s][ i ] != i + 1)
                    h_avg[ i ] ++;
            }else {//right compose with the inverse of the sample
                //h_j = 0  <=> sigma^{-1}(j) = sigma_0^{-1}(j) <=> sigma(i) = sigma_0(i) = j
                if ( sigma[ i ] != samples[ s ][ i ])
                    h_avg[ samples[ s ][ i ] - 1] ++;
            }
    for (int i = 0 ; i < n_ ; i++) h_avg[ i ] = (double) h_avg[ i ] / m ;
}

void Hamming::estimate_consensus_approx_gmm(int m, int **samples, int *sigma_0, double * best_likelihood){
    Lap     lap;
    Newton_raphson nr(n_);
    Generic gen;
    int     * sigma_0_inv = new int[ n_ ],* sigma_0_inv_neig_best = new int[ n_ ];
    int     ** freq = new int*[ n_ ];
    int     * rows=new int[n_], *cols=new int[n_], *u=new int[n_], *v=new int[n_];
    double  * h_avg = new double[ n_ ],*theta= new double[ n_ ];
    double  a1 = 0;
    for (int i = 0 ; i < n_ ; i++){ freq[ i ]= new int[n_ ]; for(int j = 0 ; j < n_ ;  j++) freq[ i ][ j ] = 0 ;}
    
    for (int i = 0 ; i < m ; i ++)
        for (int j = 0 ; j < n_ ; j ++)
            freq[ j ][ samples[ i ][ j ] - FIRST_ITEM ]--;//for the lap, minimize
    
    //freq[i,j]=sample[i,j]
    
    //int cost = -1 * lap.lap(n_, freq, rows, cols, u, v);
    lap.lap(n_, freq, rows, cols, u, v);
    for (int i = 0 ; i < n_ ; i++){
        sigma_0[ i ] = rows[i] + FIRST_ITEM;
        sigma_0_inv[ (rows[i] + FIRST_ITEM ) - FIRST_ITEM ] = i + FIRST_ITEM;
    }
    
    for (int i = 0 ; i < n_ ; i++)
        for(int j = 0 ; j < n_ ;  j++)
            freq[ i ][ j ] = freq[ i ][ j ] * -1;
    for (int i = 0 ; i < n_ ; i++){
        h_avg[ i ] = (double ) 1 - (double)freq[ sigma_0_inv[ i ] - FIRST_ITEM ][ i ] / m ;
    }
    nr.mle_theta_weighted_mallows_hamming( m , h_avg, theta);
    for ( int i = 0 ; i < n_ ; i++) a1 += h_avg[ i ] * theta[ i ] ;
    *best_likelihood = (double)-m * (a1 + log ( psi_whm(theta)));
    
    for (int i = 0 ; i < n_ ; i++) delete [] freq[ i ];
    delete [] theta;
    delete [] rows;
    delete [] h_avg;
    delete [] cols;
    delete [] u;
    delete [] v;
    delete [] freq;
    delete [] sigma_0_inv_neig_best;
    delete [] sigma_0_inv;
}
/*
void Hamming::variable_neighborhood_search(int m, int **freq, int *sigma_0_inv, double *f_eval){
    double f_eval_ini;
    Generic gen;
    do{
        f_eval_ini = (*f_eval);
       // cout<<" distance ini"<<(*f_eval)<<"\t "; gen.print_int_vector(sigma_0_inv, n_);
        local_search_swap_gmm  (m, freq, sigma_0_inv, f_eval);
       // cout<<" distance new swap "<<(*f_eval)<<"\t "; gen.print_int_vector(sigma_0_inv, n_);
        local_search_insert_gmm(m, freq, sigma_0_inv, f_eval);
       // cout<<" distance new ins "<<(*f_eval)<<"\t "; gen.print_int_vector(sigma_0_inv, n_);
    }while( (f_eval_ini) < *f_eval );//(improve);
//    cout<<"end vns"<<endl;
}

void Hamming::local_search_insert_gmm(int m, int **freq, int *sigma_0_inv, double *f_eval){
    Newton_raphson nr(n_);
    Generic gen;
    int     * sigma_0_inv_neig_best = new int[ n_ ];
    int     * next_sol_inv      = new int[ n_ ];
    double  * theta             = new double[ n_ ];
    double  * h_avg              = new double [ n_ ];
    double  a1 = 0 ;
    double  likelihood_neig  =0, likelihood_neig_best=0;
    bool halt_local_search;
    do {
        likelihood_neig_best = 0;
        halt_local_search = true;

        for(int from = 0 ; from < n_ ; from ++){
            for ( int to = 0 ; to < n_ ; to++){
                if(from != to){
                    gen.insert_at (sigma_0_inv, n_, from, to, next_sol_inv);
                    for (int i = 0 ; i < n_ ; i++)
                    {h_avg[ i ] = (double ) 1 - (double)freq[ next_sol_inv[ i ] - FIRST_ITEM ][ i ] / m ;
                    }
                    nr.mle_theta_weighted_mallows_hamming( m , h_avg, theta);
                    a1 = 0;
                    for ( int i = 0 ; i < n_ ; i++) a1 += h_avg[ i ] * theta[ i ] ;
                    likelihood_neig = (double)-m * (a1 + log ( psi_whm(theta)));
                    if (likelihood_neig > likelihood_neig_best || likelihood_neig_best == 0){
                        likelihood_neig_best = likelihood_neig;
                        for ( int i = 0 ; i < n_ ; i++) sigma_0_inv_neig_best [ i ] = next_sol_inv[ i ];
                    }
                    
                }
            }
        }
        
        if( likelihood_neig_best > *f_eval){
            *f_eval = likelihood_neig_best;
            for (int i = 0 ; i < n_ ; i ++) sigma_0_inv [ i ] = sigma_0_inv_neig_best[ i ];
            halt_local_search = false;
        }
    } while ( ! halt_local_search);
    delete [] h_avg;
    delete [] sigma_0_inv_neig_best;
    delete [] theta;
    delete [] next_sol_inv;
}

void Hamming::local_search_swap_gmm(int m, int **freq, int *sigma_0_inv, double *f_eval){
    int     * sigma_0_inv_neig_best = new int[ n_ ];
    int     *rows=new int[n_], *cols=new int[n_], *u=new int[n_], *v=new int[n_];
    double  *theta = new double [ n_ ];
    double   * h_avg = new double[ n_ ];
    double   likelihood_neig_best, likelihood_neig, a1 = 0;

    Newton_raphson nr(n_);
    Generic gen;

    for (int i = 0 ; i < n_ ; i++){
        //h_vec[ i ] = m - freq[ sigma_0_inv[ i ]][ i ];
        h_avg[ i ] = (double ) 1 - (double)freq[ sigma_0_inv[ i ] - FIRST_ITEM ][ i ] / m ;
    }
    
    bool halt_local_search;
    do{
        halt_local_search = true;
        likelihood_neig_best = 0;
        for (int i = 0 ; i < n_ - 1; i++)
            for (int j = i+1 ; j < n_ ; j++){
                    int aux = sigma_0_inv[ i ];
                    sigma_0_inv [ i ] = sigma_0_inv[ j ];
                    sigma_0_inv[ j ] = aux;
                    h_avg[ i ] = (double) 1 - (double) freq[ sigma_0_inv[ i ] - FIRST_ITEM ][ i ] / m;
                    h_avg[ j ] = (double) 1 - (double) freq[ sigma_0_inv[ j ] - FIRST_ITEM ][ j ] / m;

                    nr.mle_theta_weighted_mallows_hamming( m , h_avg, theta);
                    a1 = 0;
                    for ( int r = 0 ; r < n_ ; r++) a1 += h_avg[ r ] * theta[ r ] ;
                    likelihood_neig = (double) -m * (a1 + log ( psi_whm(theta)));
                    if (likelihood_neig > likelihood_neig_best || likelihood_neig_best == 0){
                        likelihood_neig_best = likelihood_neig;
                        for ( int i = 0 ; i < n_ ; i++) sigma_0_inv_neig_best [ i ] = sigma_0_inv[ i ];
                    }
                    
                    aux = sigma_0_inv[ i ];
                    sigma_0_inv [ i ] = sigma_0_inv[ j ];
                    sigma_0_inv[ j ] = aux;
                    h_avg[ i ] = (double) 1 - (double) freq[ sigma_0_inv[ i ] - FIRST_ITEM ][ i ] / m;
                    h_avg[ j ] = (double) 1 - (double) freq[ sigma_0_inv[ j ] - FIRST_ITEM ][ j ] / m;
                    
                }
        if ( *f_eval < likelihood_neig_best ){
            //cout<<"better sol found in swap "<<likelihood_neig_best<<endl;
            halt_local_search = false;
            *f_eval = likelihood_neig_best;
            for ( int i = 0 ; i < n_ ; i++) sigma_0_inv[ i ] = sigma_0_inv_neig_best[ i ];
        }
    }while ( ! halt_local_search) ;
    //for (int i = 0 ; i < n_ ; i++) sigma_0[ sigma_0_inv[ i ] - FIRST_ITEM ] = i + FIRST_ITEM;

    delete [] theta;
    delete [] rows;
    delete [] h_avg;
    delete [] cols;
    delete [] u;
    delete [] v;
    delete [] sigma_0_inv_neig_best;
}
*/

long double Hamming::get_likelihood(int m, int **samples, int model, int * sigma_0){
    Newton_raphson nr(n_);
    long double likelihood ;
    double  * theta = new double[ n_ ];
    double  psi, a1 = 0;
    if(model == MALLOWS_MODEL){
        int dist = distance_to_sample(samples, m, sigma_0);
        theta[ 0 ] = nr.Newton_raphson_method( (double) dist/m, 0.0, HAMMING_DISTANCE, MALLOWS_MODEL, -1, NULL);
        psi = psi_hm(theta[ 0 ]);
        likelihood =  - theta[0] * dist - m * log( psi );
    }else{
        Generic gen;
        double  * h_avg = new double[ n_];
        sample_to_h_vector(samples, m, sigma_0, h_avg);
        nr.mle_theta_weighted_mallows_hamming( m , h_avg, theta);
        for ( int i = 0 ; i < n_ ; i++) a1 +=(double)h_avg[i]  * theta[ i ] ;
        delete [] h_avg;
        likelihood = (long double)-m * (a1 + log ( psi_whm(theta)));
    }
    delete [] theta;
    return likelihood;
}

long double Hamming::get_likelihood_from_h(int m, int model, double *theta, double * h_avg){
    long double likelihood ;
    double psi, a1 = 0;
    if(model == MALLOWS_MODEL){
        int dist = 0;
        for (int i = 0 ; i < n_ ;i ++) dist += h_avg[ i ];
        dist *= m;
        psi = psi_hm(theta[ 0 ]);
        likelihood =  - theta[0] * dist - m * log( psi );
    }else{
        Generic gen;
        int   *h_sum = new int[ n_];
        for ( int i = 0 ; i < n_ ; i++) h_sum[ i ] = h_avg[ i ]* m;
        for ( int i = 0 ; i < n_ ; i++) a1 +=(double)h_sum[i] / m * theta[ i ] ;
        delete [] h_sum;
        likelihood = (long double)-m * (a1 + log ( psi_whm(theta)));
    }
    return likelihood;
}


void Hamming::estimate_theta(int m, int *sigma_0, int **samples, int model, double *theta){
    Newton_raphson nr(n_);
    if (model == MALLOWS_MODEL){
        //Newton_raphson nr(n_);
        double  dist = 0;
        dist = distance_to_sample( samples, m, sigma_0);
        *theta = nr.Newton_raphson_method( (double) dist/m, 0.0, HAMMING_DISTANCE, MALLOWS_MODEL, -1, NULL);
    }
    else{
        Generic gen;
        double   *h_avg         = new double [ n_ ];
        sample_to_h_vector(samples, m, sigma_0, h_avg);
        nr.mle_theta_weighted_mallows_hamming( m , h_avg, theta);
        //for ( int j = 0 ; j < n_ ; j ++) cout<< theta[j]<<" "; cout<<" theta (estimate_theta_gmm) "<<endl;
        delete [] h_avg;
    }
}















