//
//  Ferrers_diagram.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 08/07/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Ferrers_diagram__
#define __perms_mallows__Ferrers_diagram__

#include <iostream>
#include <vector>
#include <cmath>
#include "Generic.h"

class Ferrers_diagram;
using namespace std;
class Ferrers_diagram{
    
private:
    int     partition_of_n_;
    int     *ferrers_shape_;
    int     ferrers_shape_length_;

    long double hook_length_;
    long double number_of_permutations_;

    int     **syt_;//greene_niejenhuis_wilf
    int     *ferrers_shape_dynamic_; // for the random generation of SYT
    int     ferrers_shape_length_dynamic_;

public:
    Ferrers_diagram(int n, int*partition, int k){
        partition_of_n_         = n;
        ferrers_shape_length_   = k;
        ferrers_shape_          = partition;
        hook_length_            = -1;
        number_of_permutations_ = -1;
        syt_                    = NULL;
        ferrers_shape_dynamic_  = NULL; //
    }
    ~Ferrers_diagram(){
        delete [] ferrers_shape_;
        if ( syt_ != NULL){
            for ( int i = 0 ; i < ferrers_shape_length_ ; i ++)
                delete [] syt_[ i ];
            delete [] syt_;
            delete [] ferrers_shape_dynamic_;
        }
    }
    int ** get_syt(){return syt_;}
    void walk(int rand_cell, int *cell_i, int* cell_j);
    
    void init_tables_for_random_SYT_generation(){
        ferrers_shape_length_dynamic_ = ferrers_shape_length_ ;
        syt_                   = new int * [ferrers_shape_length_ ];
        ferrers_shape_dynamic_ = new int [ferrers_shape_length_ ];
        for ( int i = 0 ; i < ferrers_shape_length_  ; i++){
            ferrers_shape_dynamic_[ i ] = ferrers_shape_[ i ];
            syt_[ i ] = new int[ferrers_shape_[i]  ];
            for (int j = 0 ; j < ferrers_shape_[i]  ; j++ )syt_[i][j]=0;
        }
    }
    
    void delete_tables_for_random_SYT_generation(){
            for ( int i = 0 ; i < ferrers_shape_length_ ; i ++)
                delete [] syt_[ i ];
            delete [] syt_;
            delete [] ferrers_shape_dynamic_;
    }

    int get_num_cols(){return ferrers_shape_[0];}//equivalent to the LIS of the gen permu
    int get_num_rows(){return ferrers_shape_length_;}//
    int get_resulting_distance(){return partition_of_n_ - ferrers_shape_[0];}//equivalent to the distance of the gen permu
    
    int get_ferrers_shape_length(){return ferrers_shape_length_;}
    int*get_ferrers_shape(){return ferrers_shape_;}
    
    long double get_num_permus(){return number_of_permutations_;}
    
    long double get_hook_len(){return hook_length_;}
    
    long double calculate_hook_length(long double n_factorial);
    
    void random_SYT();
    
    int get_num_cells_down(int i , int j);
    
    /*void SYT_toStr(){
        for ( int i = 0 ; i < ferrers_shape_length_  ; i++){
            for ( int j = 0 ; j < ferrers_shape_[i]  ; j++)
                cout<< syt_[i ][ j ]<<" ";
            cout<<endl;
        }
    }
    
    void toStr(){
        Generic gen;
        cout<<"Partition of "<<partition_of_n_<<": ";
        gen.print_int_vector(ferrers_shape_, ferrers_shape_length_);
        cout<<"\t hook length  "<<hook_length_<<"\n";
    }*/

};

#endif /* defined(__perms_mallows__Ferrers_diagram__) */
