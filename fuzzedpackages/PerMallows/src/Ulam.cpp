//
//  Ulam.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 08/07/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//
#include <R.h>
#include "Ulam.h"

#include <vector>
#include <set>
double Ulam::probability(int *s, int *s_0, double * theta){
    int dist = distance(s, s_0);
    double *proba = new double[ n_ ];
    calculate_probas_at_each_distance( theta[ 0 ], proba);
    double prob =  exp(-dist * theta[0] )/proba[ n_ - 1 ];
    delete [] proba;
    return prob;
}

int Ulam::search_LIS(int* M, int* A, int i, int L ) {
    int j = 0;
    int k = L-1;
    while( j <= k ) {
        int m = ( j + k ) / 2;
        if( A[M[m]] <= A[i] ) j = m + 1;
        else k = m - 1;
    }
    
    return k;
}

int Ulam::longest_increasing_subsequence(int*sigma){
    // http://robentan.blogspot.com.es/2011/11/more-efficient-algorithm-for-longest.html
    M[0] = 0;
    P[0] = -1;
    int L = 1;
    
    for(int i=1; i<n_; ++i) {
        int j = search_LIS( M, sigma, i, L );
        if( j == -1 ) P[i] = -1;
        else P[i] = M[j];
        
        if( j == L-1 || sigma[i] < sigma[M[j+1]] ) {
            M[j+1] = i;
            if( j+2 > L ) L = j+2;
        }
    }
    return L;
 
}

int Ulam::distance(int *s1, int *s2){
    //int*comp=new int[n_], *inv = new int[ n_];
    for(int i = 0 ; i < n_ ; i ++) inv_[ s2[ i ] - 1 ] = i + 1;
    for (int i=0; i<n_ ; i++)   comp_[i]=s1[inv_[i]-1];
    int dist = n_ - longest_increasing_subsequence(comp_);
    //delete [] comp;
    //delete [] inv;
    return dist;
}

int Ulam::gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k){
    int j; //index
    //test for special cases
    if(n == 0)    {
        (*k) = 0;
        return(GEN_EMPTY);
    }
    //initialize: vector[0] = n, vector[1, ..., n - 1] are 1
    vector[0] = n;
    for(j = 1; j < n; j++)        vector[j] = 1;
    (*k) = 1;
    return(GEN_NEXT);
}

int Ulam::gen_part_next(unsigned char *vector, unsigned char *k, int bound){
    //bound == 0 => no bound
    static int j = 0; //remember index of the rightmost part which is greater than 1
    int        r;     //temporary remainder
    int        temp;  //auxiliary element
    
    //easy case
    if(j >= 0 && vector[j] == 2){ // j >= 0 &&  mod to avoid illegal read
        vector[j] = 1;
        j--;
        (*k)++;
        //terminate if the num of columns is smaller than 'bound'
        /*if ((int)vector[0] < bound){
            j=0;
            return (GEN_TERM);
        }*/
        return(GEN_NEXT);
    }
    //terminate if all parts are 1
    if(vector[0] == 1){
        j = 0;
        return(GEN_TERM);
    }
    //decrease
    vector[j]--;
    temp = vector[j];
    r = *k - j;
    //set right-hand elements
    while(r > temp){
        j++;
        vector[j] = temp;
        r -= temp;
    }
    *k = j + 2;
    //set rightmost element
    if(r > 1){
        j++;
        vector[j] = r;
    }
    //terminate if the num of columns is smaller than 'bound'
    /*if ((int)vector[0] < bound){
        j=0;
        return (GEN_TERM);
    }*/
    return(GEN_NEXT);
}

void Ulam::random_sample_at_dist(int dist, int m, int **samples){
    fill_shapes_of_n();
    for (int s = 0 ; s < m ; s++) {
        samples[ s ] = new int[ n_ ];
        generate_permu_with_given_LIS(n_ - dist , samples[ s ]);
    }
}

int Ulam::set_median(int m, int **samples, int *sigma_0){
    //TODO
    /*int min_dist = n_ * m;
    int set_median = -1;
    for(int i = 0; i < m ; i++){
        int dist = distance_to_sample(samples, m, samples[ i ]);
        if (dist < min_dist) {
            min_dist = dist;
            set_median = i;
        }
    }
    for(int i = 0; i < n_ ; i++) sigma_0[i] = samples[set_median][i];
    return min_dist;*/
	int min_dist = 0, min_pos = 0 ;
	int * sum_dist = new int[ m ];
	for(int i = 0; i < m ; i++) sum_dist [ i ] = 0;
    
	for(int i = 0; i < m - 1 ; i++)
		for(int j = i + 1 ; j < m ; j++){
			int d = distance(samples[ i ], samples[ j ]);
			sum_dist[ i ] += d;
			sum_dist[ j ] += d;
		}
	min_dist = sum_dist[ 0 ];
	min_pos = 0 ;
	for(int i = 1; i < m ; i++) if ( min_dist > sum_dist[ i ] ) {min_dist = sum_dist[ i ]; min_pos = i;}
	for(int i = 0; i < n_ ; i++) sigma_0[ i ] = samples[ min_pos][ i ];
	delete [] sum_dist;
	return min_dist;
}

void Ulam::estimate_theta(int m, int *sigma_0, int **samples, int model, double *theta){
    Newton_raphson newton(n_);
    fill_shapes_of_n();
    int dist_avg = distance_to_sample(samples, m, sigma_0);
    *theta = newton.Newton_raphson_method( (double)dist_avg/m, -1.001,ULAM_DISTANCE, MALLOWS_MODEL, -1, num_permus_per_dist_);
}

long double Ulam::get_likelihood(int m, int** samples, int model, int * sigma_0) {
    Newton_raphson newton(n_);
    long double likelihood, psi = 0;
    int     dist_avg = 0 ;
    double  theta;

    fill_shapes_of_n();
    dist_avg = distance_to_sample(samples, m, sigma_0);
    theta    = newton.Newton_raphson_method( (double)dist_avg/m, -1.001,ULAM_DISTANCE, MALLOWS_MODEL, -1, num_permus_per_dist_);
    for (int i = 0 ; i < n_ ; i ++ ) psi += num_permus_per_dist_[ i ] * exp (-theta * i );
    likelihood = - dist_avg * theta - m* log ( psi );
    return likelihood;
}

int Ulam::distance_to_sample(int **samples, int m, int *sigma){
    int dist= 0;
    int *comp = new int[ n_ ], *sigma_inv = new int[ n_ ];
    for(int j = 0 ; j < n_ ; j ++) sigma_inv[sigma[ j ] - 1 ] = j + 1;
    for(int s = 0 ; s < m ; s ++){
        for(int i = 0 ; i < n_ ; i ++) comp[ i ] = samples[ s ][ sigma_inv [ i ] - 1 ];
        dist += ( n_ - longest_increasing_subsequence(comp));//distance(comp);
    }
    delete []sigma_inv;
    delete []comp;
    return dist ;
}

void Ulam::fill_shapes_of_n(){
    /* adapted from ZS1: Zoghbi, Antoine and Stojmenovic, Ivan: Fast Algorithms for Generating Integer Partitions. International Journal of Computer Mathematics, 70, 1998, 319-332.
     generaes partitions in  anti-lexicographic order */
    
    if (shapes_of_n_->size()==0 ){
        unsigned char k;                  //length of figures
        unsigned char *vector     = NULL; //where the current figure is stored
        int           gen_result;         //return value of generation functions
        //int           x;                //iterator
        num_partitions_of_n_ = 0;
        int           prev_distance = -1 , dist ;
        long double   cont = 0 ;
        int           part_len;
        
        //alloc memory for vector
        /*vector = (unsigned char *)malloc(sizeof(unsigned char) * n_ );
        if(vector == NULL)    {
            fprintf(stderr, "error: insufficient memory\n");
            exit(EXIT_FAILURE);
        }*/
        vector = new unsigned char [ n_ ];
        Ferrers_diagram*f;
        gen_result = gen_part_init(vector, n_, &k);
        while(gen_result == GEN_NEXT ) {//&& (int)vector[0] >= bound
            //for(x = 0; x < k; x++)            printf("___%u ", vector[x]); cout<<endl;
            
            part_len = (int)k;
            int*part = new int[part_len];//DO NOT delete, member of Ferrers_diagram
            for(int i = 0 ; i < part_len; i++) part[i]=(int)vector[i];
            f = new Ferrers_diagram(n_, part , part_len);
            shapes_of_n_->push_back(f);
            num_partitions_of_n_ ++;
            f->calculate_hook_length(facts_[ n_ ]);
            dist = f->get_resulting_distance();
            num_permus_per_dist_[ dist ] += f->get_num_permus();
            if ( dist != prev_distance ){
                first_index_at_dist_ [ dist ] = cont;
                num_permus_at_shape_acumul_.push_back( f->get_num_permus()) ;
                //cout<<"Generating shape at distance "<<dist<<endl;
            }else
                num_permus_at_shape_acumul_.push_back( num_permus_at_shape_acumul_.at( cont - 1 ) + f->get_num_permus() ) ;
            prev_distance = dist;
            cont ++;
            gen_result = gen_part_next(vector, &k, 0);
        }
        delete [] vector;
        //free(vector);
    }
}



void Ulam::generate_permu_with_given_LIS(int l, int *sigma){//TODO cambiar a 2 syt como out-param
    //STEP 1
    int     d = n_ - l;
    int     to_insert;
    int     col, row, aux, new_col, new_row;
    int     *col_index = new int[n_], *row_index = new int[n_];
    long double fs_index ;//ferrer shape index
    int fs_length;//ferrer shape len
    //double permus = (double)rand() / (double)(RAND_MAX) * num_permus_per_dist_[ d ];
    double permus = unif_rand() * num_permus_per_dist_[ d ];
    fs_index = first_index_at_dist_[ d ];
    while(num_permus_at_shape_acumul_[ fs_index ] <= permus){
        fs_index++;
        //cout<<num_permus_at_shape_acumul_[ target_shape ]<<endl;
    }
    fs_length = shapes_of_n_->at(fs_index)->get_ferrers_shape_length();
    
    //STEP 2
    int* shape1 = new int [fs_length];//member of f1
    int* shape2 = new int [fs_length];
    for (int i = 0 ; i < shapes_of_n_->at(fs_index)->get_ferrers_shape_length() ; i ++){
        shape1[ i ] = shapes_of_n_->at(fs_index)->get_ferrers_shape()[ i ];
        shape2[ i ] = shapes_of_n_->at(fs_index)->get_ferrers_shape()[ i ];
    }
    Ferrers_diagram* f1 = new Ferrers_diagram(n_   , shape1 , fs_length);
    Ferrers_diagram* f2 = new Ferrers_diagram(n_   , shape2 , fs_length);
    
    f1->random_SYT();
    f2->random_SYT();
    
    
    //STEP 3
    int ** tableau1 = f1->get_syt();//
    int ** tableau2 = f2->get_syt();
    for (int i = 0 ; i < f2->get_ferrers_shape_length() ; i++){
        for (int j =  0 ; j < f2->get_ferrers_shape()[i] ; j++) {
            row_index[ tableau2[ i ][ j ] - 1 ] = i;
            col_index[ tableau2[ i ][ j ] - 1 ] = j;
        }
    }
//    Generic gen;
//    gen.print_int_matrix(tableau1, f1->get_num_rows(), f1->get_num_cols());
//    gen.print_int_matrix(tableau2, f1->get_num_rows(), f1->get_num_cols());
    //Schensted algorithm
    for (int index = n_ - 1 ; index >= 0 ; index --){
        col = col_index[ index ];
        row = row_index[ index ];
        to_insert = tableau1[ row ][ col ];
        while (row != 0) {
            new_col=0, new_row = row - 1;
            while ( f1->get_ferrers_shape()[new_row] > new_col+1
                    && tableau1[new_row][new_col + 1 ] < to_insert)
                new_col++;
            aux = tableau1[new_row][new_col];
            tableau1[new_row][new_col] = to_insert;
            to_insert = aux;
            row = new_row;
            col = new_col;
        }
        sigma[ index ] = to_insert;
        tableau1[row_index[ index ]][ col_index[ index ]] = n_ + 1;
//        gen.print_int_vector(sigma, n_);
//        gen.print_int_matrix(tableau1, f1->get_num_rows(), f1->get_num_cols());
//        gen.print_int_matrix(tableau2, f1->get_num_rows(), f1->get_num_cols());
    }
    
    
    delete [] col_index;
    delete [] row_index;
    delete f1;
    delete f2;
}

void Ulam::distances_sampling(int m, double theta, int **samples){
    double  distance_acum = 0, rand_distance = 0;
    double  *proba_acumul        = new double[ n_ ];
    int     target_distance;
    
    fill_shapes_of_n();
    //calculate_probas_at_each_distance( theta, proba, bound_l);
    proba_acumul[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < n_ ; i++)//acumulate the number of permus at each distance
        proba_acumul[i] = num_permus_per_dist_[i] * exp ( -theta * i ) + proba_acumul[ i - 1 ];
    for (int i = 0 ; i < m ; i ++){
        //rand_distance = (double) rand() / (double)(RAND_MAX) * proba_acumul[ n_ - 1 ];
        rand_distance = unif_rand()* proba_acumul[ n_ - 1 ];
        target_distance = 0;
        while(proba_acumul[ target_distance ] <= rand_distance) target_distance++;
        samples[ i ] = new int[ n_ ];
        //cout<<"ulam distance greene_niejenhuis_wilf "<<target_distance<<endl;
        generate_permu_with_given_LIS( n_ - target_distance, samples[i]);
        distance_acum += (double)target_distance;
    }
//cout<<"end"<<endl;
    //cout<<"average distance (dist_sampling) "<<distance_acum/m<<endl;
    delete [] proba_acumul;
}

void Ulam::gibbs_sampling(int m, double *theta, int model, int **samples){
    //fill_shapes_of_n();
    
    int     burning_period_samples = n_*log(n_);
    int     *sigma = new int[ n_ ];
    int     *sigma_prime = new int[ n_ ];
    Generic  gen;
    //void Generic::insert_at(int *sigma, int n, int move, int to, int*res){
    gen.generate_random_permutation( n_ , 1, sigma);
    
    for(int sample = 0 ; sample < m + burning_period_samples ; sample ++){
        int a;
        int b;
        do {
            a = (int) (unif_rand() * n_);//rand() % ( n_ );
            b = (int) (unif_rand() * n_);//rand() % ( n_ );
        }while (a == b);
        gen.insert_at(sigma, n_ , a , b , sigma_prime);
        bool make_swap = false;
        if(  distance(sigma) > distance(sigma_prime) )  make_swap = true;
        else{
            double rand_double = unif_rand();//(double)rand()/RAND_MAX;
                if(rand_double < exp(-theta[0])) make_swap = true;
        }
        if(make_swap){
            for(int i = 0  ; i < n_ ; i ++) sigma[ i ] = sigma_prime[ i ];
        }
        if(sample>=burning_period_samples){
            samples[sample-burning_period_samples]=new int[ n_ ];
            for(int i = 0  ; i < n_ ; i ++)   samples[ sample - burning_period_samples ][ i ] = sigma[ i ];
        }
    }
    delete [] sigma_prime;
}

double Ulam::expectation(double theta){
    fill_shapes_of_n();
    long double numer = 0, denom = 0;
    for (int d = 0 ; d < n_ - 1; d++){
        long double aux = num_permus_per_dist_[ d ] * exp(-theta *d ) ;
        numer += aux * d;
        denom += aux;
    }
    return (double)numer / denom;
}



long double Ulam::num_permus_at_distance_approx(int d){
    //gordon, a measure of agreement between arrays
    return facts_[d] * pow ((facts_[n_]/(facts_[d]* facts_[n_-d] )),2);
}

long double Ulam::num_permus_at_distance(int d){
    fill_shapes_of_n();
    return num_permus_per_dist_[ d ];
}



double Ulam::psi(double theta){
    double  *proba = new double[ n_ ];
    calculate_probas_at_each_distance( theta, proba);
    return proba[ n_ - 1 ];
}

void Ulam::calculate_probas_at_each_distance(double theta, double *proba){
    fill_shapes_of_n();
    //calculate_probas_at_each_distance( theta, proba, bound_l);
    proba[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < n_ ; i++)//acumulate the number of permus at each distance
        proba[i] = num_permus_per_dist_[i] * exp ( -theta * i ) + proba[ i - 1 ];
    ///end
    
}

