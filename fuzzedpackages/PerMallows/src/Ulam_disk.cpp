//
//  Ulam_disk.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 31/01/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//


#include "Ulam_disk.h"
#include <stdlib.h>
#include <stdio.h>
#include <R.h>

void    Ulam_disk::read_permus_per_dist(){
    if ( num_permus_per_dist_[ 0 ] == 0 ){//not initialized
        //num_permus_per_dist_num_permus_per_dist_ = new long double[ n_ ];
        char integer_string[5];
        sprintf(integer_string, "%d", n_);
        char str_permus_per_dist[600];
        strcpy(str_permus_per_dist, str_base_path);
        strcat(str_permus_per_dist, "permus_per_dist_");
        strcat(str_permus_per_dist, integer_string);
        ifstream  file;
        file.open(str_permus_per_dist);
        if(!((file ))){
            //cout << "Cannot read input file. Have you generated the files? : "<<str_permus_per_dist<<endl;
            //checked from the R code
            return;
        }
        for (int i = 0 ; i < n_ ; i ++)
            file >> num_permus_per_dist_[ i ];
        file.close();
    }
}


void Ulam_disk::save_counts_to_file_bin(){
    /* adapted from ZS1: Zoghbi, Antoine and Stojmenovic, Ivan: Fast Algorithms for Generating Integer Partitions. International Journal of Computer Mathematics, 70, 1998, 319-332.
     generaes partitions in  anti-lexicographic order */
    ofstream file_permus_per_dist;
//    ofstream file_permus_per_shape;
    FILE *file_permus_per_shape = NULL;
    //file = fopen(str_permus_per_shape, "r");
    
    char integer_string[5];
    sprintf(integer_string, "%d", n_);
    char str_permus_per_shape[500] ;//= "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_shape_"; // un file e estos por cada n y distance
    strcpy(str_permus_per_shape, str_base_path);
    strcat(str_permus_per_shape, "permus_per_shape_bin_");
    
    char str_permus_per_shape_n_d[600];
    char str_permus_per_dist[600] ;//= "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_dist_";
    strcpy(str_permus_per_dist, str_base_path);
    strcat(str_permus_per_dist, "permus_per_dist_");
    
    strcat(str_permus_per_shape, integer_string);
    strcat(str_permus_per_dist, integer_string);
    
    unsigned char k;                  //length of figures
    unsigned char *vector     = NULL; //where the current figure is stored
    int           gen_result;         //return value of generation functions
    int           prev_distance = -1 , dist ;
    int           part_len;
    long double     permus_per_shape = 0;
    vector = (unsigned char *)malloc(sizeof(unsigned char) * n_ );
    if(vector == NULL)    {
        //ERROR fprintf(stderr, "error: insufficient memory\n");exit(EXIT_FAILURE);
    }
    for (int i = 0 ; i < n_ ;i++) num_permus_per_dist_[ i ] = 0 ;
    
    Ferrers_diagram*f;
    gen_result = gen_part_init(vector, n_, &k);
    int * info = new int [ n_ + 3 ];
    while(gen_result == GEN_NEXT ) {
        part_len = (int)k;
        int*part = new int[part_len];//DO NOT delete, member of Ferrers_diagram
        for(int i = 0 ; i < part_len; i++) part[i]=(int)vector[i];
        dist = part[ 0 ];
        
        f = new Ferrers_diagram(n_, part , part_len);
        f->calculate_hook_length(facts_[ n_ ]);
        dist = f->get_resulting_distance();
        num_permus_per_dist_[ dist ] += f->get_num_permus();
        if ( dist != prev_distance ){
            permus_per_shape = 0;
            //cout<<"Generating shape at distance "<<dist<<endl;
//            if (file_permus_per_shape.is_open() ) file_permus_per_shape.close();
            if (file_permus_per_shape != NULL ) fclose(file_permus_per_shape);
            strcpy(str_permus_per_shape_n_d, str_permus_per_shape);
            strcat(str_permus_per_shape_n_d, "_");
            sprintf(integer_string, "%d", dist );
            strcat(str_permus_per_shape_n_d, integer_string);
//            file_permus_per_shape.open(str_permus_per_shape_n_d);
            file_permus_per_shape = fopen(str_permus_per_shape_n_d, "w");
        }
        permus_per_shape += f->get_num_permus();
        /*file_permus_per_shape << permus_per_shape <<" ";
        for ( int i =  0; i < part_len ; i++ )file_permus_per_shape<<part[ i ]<<" ";
        file_permus_per_shape<<endl;*/
        
        int np_exp = (int) log10(permus_per_shape);
        long double aux = pow(10, np_exp);
        int np_significant_left = (int) floor( (double) permus_per_shape / aux);
        double np_significant_dou =  (double) permus_per_shape / aux - np_significant_left;
        int np_significant_right = np_significant_dou * 1000000;
        info [0] = np_significant_left;
        info [1] = np_significant_right;
        info [2] = np_exp;
        for (int i = 0 ; i < part_len ; i++)
            info[ i + 3 ] = part[ i ];
        for (int i = part_len + 3 ; i < n_ + 3 ; i++) info[ i ] = 0;

        //fwrite(info, sizeof(int [n_ + 3 ]), 1, file_permus_per_shape);
        fwrite(info, sizeof(int) * ( n_ + 3 ), 1, file_permus_per_shape);
        
        prev_distance = dist;
        gen_result = gen_part_next(vector, &k, 0);
        
        delete f ;
    }
    free(vector);
    delete [] info;
    fclose(file_permus_per_shape);

    //cout<<str_permus_per_dist<<endl;
    file_permus_per_dist.open (str_permus_per_dist);
    for (int i = 0 ; i < n_; i++ ){file_permus_per_dist<< num_permus_per_dist_[ i ]<<endl ;}
    file_permus_per_dist.close();
}


double Ulam_disk::expectation(double theta){
    read_permus_per_dist();//OJO WARNING - si es de clase
    long double numer = 0, denom = 0;
    for (int d = 0 ; d < n_ - 1; d++){
        long double aux = num_permus_per_dist_[ d ] * exp(-theta *d ) ;
        numer += aux * d;
        denom += aux;
    }
    return (double)numer / denom;
}

void Ulam_disk::distances_sampling(int m, double theta, int **samples){
    double  rand_distance = 0;
    int     target_distance, num_perm_at_dist_to_gen;
    double  *proba_acumul       = new double[ n_ ];
    double  *m_random_numbers   = new double[ m ];
    int     *m_target_distances = new int [ n_ ]; for (int i = 0 ; i < n_ ; i ++) m_target_distances[ i ] = 0;
    int     *m_shapes_lengths   = new int[ m ];
    int     **m_shapes           = new int*[ m ]; for (int i=0;i< m ;i++) m_shapes[i]= new int[ n_ ];
    
    read_permus_per_dist();//OJO WARNING - si es de clase
    proba_acumul[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < n_ ; i++)//acumulate the number of permus at each distance
        proba_acumul[i] = num_permus_per_dist_[i] * exp ( -theta * i ) + proba_acumul[ i - 1 ];
    
    for (int s=0;s< m ;s++){
        //rand_distance = ((double)rand() / (double)(RAND_MAX)) * proba_acumul[ n_ -1];
        rand_distance = unif_rand() * proba_acumul[ n_ -1];
        //cout<<"rand_dist "<<rand_distance<<endl;
        target_distance = 0;
        while(proba_acumul[ target_distance ] < rand_distance)//ekhine pone <=
            target_distance++;
            m_target_distances[target_distance]++;
    }
    
    
    int s=0;
    for (target_distance=0;target_distance < n_  ;target_distance++){
        if (m_target_distances[target_distance] != 0){
            num_perm_at_dist_to_gen = m_target_distances[target_distance];
            for (int j=0;j<num_perm_at_dist_to_gen;j++){
                //m_random_numbers[j] = ((double)rand() / (double)(RAND_MAX)) * (double)num_permus_per_dist_[ target_distance ];
                m_random_numbers[j] = unif_rand() * (double)num_permus_per_dist_[ target_distance ];
                m_shapes_lengths[j]=0;
            }
            if (num_perm_at_dist_to_gen > 1)
                sort(m_random_numbers,m_random_numbers + num_perm_at_dist_to_gen);
            
            read_mutiple_shapes_from_file_bin (target_distance, m_random_numbers,num_perm_at_dist_to_gen, m_shapes, m_shapes_lengths);
            
            for (int j=0;j<num_perm_at_dist_to_gen;j++){
                samples[ s ] = new int[ n_ ] ;
                generate_permu_with_given_LIS( n_ - target_distance, samples[ s ], m_shapes[j], m_shapes_lengths[j]);
                s++;
            }
        }
    }
    
    delete [] proba_acumul;
    delete [] m_random_numbers;
    delete [] m_target_distances;
    delete [] m_shapes_lengths;
    for (int i = 0 ; i < m ; i ++) delete [] m_shapes[ i ];
    delete [] m_shapes;
}

void Ulam_disk::read_mutiple_shapes_from_file_bin(int d, double * bounds, int num_bounds, int ** shapes, int*shape_len){
    
    char integer_string[5];
    
    sprintf(integer_string, "%d", n_ );
    char str_permus_per_shape[500] ;
    strcpy (str_permus_per_shape, str_base_path);
    strcat(str_permus_per_shape, "permus_per_shape_bin_");
    strcat(str_permus_per_shape, integer_string);
    strcat(str_permus_per_shape, "_");
    sprintf(integer_string, "%d", d);
    strcat(str_permus_per_shape, integer_string);
    
    //ifstream  file;
    FILE *file;
    file = fopen(str_permus_per_shape, "r");
    if (file == NULL){
        //cout << "Cannot read input file - permus_per_dist  at read_mutiple_shapes_from_file: "<<str_permus_per_shape<<endl;
        //exit(1);
    }
    double  bound;
    int     bound_exp;
    double  bound_significant;
    bool    line_found = false;
    bool    aux1 = false;
    int     lines_read  = 0;
    //int     lines_in_block = 128; global
    //global line_block_ to store the read data
    int     * line_mid, *line_pre, * line_last;
    int     ini = 0, end , mid;
    
    for (int i=0;i<num_bounds;i++){
        bound=bounds[i];
        //bound = bound_significant e+ bound_exp
        bound_exp = (int) log10(bound);
        bound_significant = (double) bound / (double) pow(10.0, bound_exp);
        
        
        do{
            line_last = line_block_ + ( (lines_read - 1) * (n_ + 3 ) );
            //            cout<<line_block_[0]<<"."<<line_block_[1]<<"e+"<<line_block_[2]<<"\t"<<line_last[0]<<"."<<line_last[1]<<"e+"<<line_last[2]<<endl;
            if (lines_read > 0 && (line_last[ 2 ] > bound_exp || (line_last[ 2 ] == bound_exp && (double) line_last[0]+(0.000001*line_last[1]) > bound_significant))){
                //in block
                ini = 0, end = lines_read-1;
                line_found = false;
                do {
                    
                    mid = (end - ini )/ 2 + ini;
                    line_mid = line_block_ + ( (mid ) * (n_ + 3 ) );
                    line_pre = line_block_ + ( (mid - 1) * (n_ + 3 ) );
                    aux1 = (line_mid[ 2 ] > bound_exp || (line_mid[ 2 ] == bound_exp && (double) line_mid[0]+(0.000001*line_mid[1]) > bound_significant) );
                    if (aux1  &&
                        ( mid == 0 ||
                        (line_pre[ 2 ] < bound_exp || (line_pre[ 2 ] == bound_exp && (double) line_pre[0]+(0.000001*line_pre[1]) < bound_significant ) ))){
                        line_found = true;
                    }else if ( aux1 )
                        end = mid - 1;
                    else
                        ini = mid + 1;
                } while ( !line_found );
                //cout<<"Found "<<line_pre[0]+(0.000001*line_pre[1])<<" "<<bound_significant<<" "<<line_mid[0]+(0.000001*line_mid[1])<<";"<<line_pre[2]<<";"<<line_mid[2]<<";"<<bound_exp<<endl;
            }else
                //lines_read = (int) fread(line_block_  ,sizeof(int [ n_ + 3]), lines_in_block_ , file);
                lines_read = (int) fread(line_block_, sizeof(int) * (n_ + 3), lines_in_block_, file);
        }while ( ! line_found );
        
        while (line_mid[shape_len[i] + 3] != 0  && shape_len[i] < n_ ){
            shapes[i][ shape_len[i] ] = (int) line_mid[shape_len[i]+3];
            (shape_len[i])++;
        }
    }
    //delete [] line_block;
    fclose(file);
}


void Ulam_disk::generate_permu_with_given_LIS (int l, int *sigma, int* shape, int shape_len){
    int     to_insert;
    int     col, row, aux, new_col, new_row;
    int     *col_index = new int[n_], *row_index = new int[n_];
    
    int* shape1 = new int [ shape_len ];//member of f1
    int* shape2 = new int [ shape_len ];
    memcpy(shape1,shape,sizeof(int)*shape_len);
    memcpy(shape2,shape,sizeof(int)*shape_len);
    
    Ferrers_diagram * f1 = new Ferrers_diagram( n_, shape1 , shape_len);
    Ferrers_diagram * f2 = new Ferrers_diagram( n_ , shape2 , shape_len);
    
    f1->random_SYT();
    f2->random_SYT();
    int ** tableau1 = f1->get_syt();
    int ** tableau2 = f2->get_syt();
    
    
    for (int i = 0 ; i < f2->get_ferrers_shape_length() ; i++){
        for (int j =  0 ; j < f2->get_ferrers_shape()[i] ; j++) {
            row_index[ tableau2[ i ][ j ] - 1 ] = i;
            col_index[ tableau2[ i ][ j ] - 1 ] = j;
        }
    }
    
    for (int index = n_ - 1 ; index >= 0 ; index --){
        col = col_index[ index ];
        row = row_index[ index ];
        to_insert = tableau1[ row ][ col ];
        while (row != 0) {
            new_col=0, new_row = row - 1;
            while (f1->get_ferrers_shape()[new_row] > new_col+1
                   && tableau1[new_row][new_col + 1 ] < to_insert)
                new_col++;
            aux = tableau1[new_row][new_col];
            tableau1[new_row][new_col] = to_insert;
            to_insert = aux;
            row = new_row;
            col = new_col;
        }
        sigma[index ] = to_insert;
        tableau1[row_index[ index ]][ col_index[ index ]] = n_ + 1;
        //gen.print_int_vector(sigma, n_);
    }
    delete [] col_index;
    delete [] row_index;
    //delete [] shape;
    delete f1;
    delete f2;
}


long double Ulam_disk::get_likelihood(int m, int **samples, int model, int *sigma_0){
    Newton_raphson newton(n_);
    long double likelihood, psi = 0;
    int     dist_avg = 0 ;
    double  theta;
    
    read_permus_per_dist();//OJO WARNING - si es de clase
    dist_avg = distance_to_sample(samples, m, sigma_0);
    theta    = newton.Newton_raphson_method( (double)dist_avg/m, -1.001,ULAM_DISTANCE, MALLOWS_MODEL, -1, num_permus_per_dist_);
    for (int i = 0 ; i < n_ ; i ++ ) psi += num_permus_per_dist_[ i ] * exp (-theta * i );
    likelihood = - dist_avg * theta - m* log ( psi );
    return likelihood;

}

void Ulam_disk::estimate_theta(int m, int *sigma_0, int **samples, int model, double *theta){
    Newton_raphson newton(n_);
    read_permus_per_dist();//OJO WARNING - si es de clase
    int dist_avg = distance_to_sample(samples, m, sigma_0);
    *theta = newton.Newton_raphson_method((double)dist_avg/m, -1.001,ULAM_DISTANCE, MALLOWS_MODEL, -1, num_permus_per_dist_);
}















