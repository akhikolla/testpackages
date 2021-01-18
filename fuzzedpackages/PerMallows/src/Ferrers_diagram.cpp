//
//  Ferrers_diagram.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 08/07/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include "Ferrers_diagram.h"
#include <R.h>
void Ferrers_diagram::random_SYT(){
    //￼A Probabilistic Proof of a Formula for the Number of Young Tableaux of a Given Shape
    //by Greene, Nijenhuis and Wilf
    int     n_count = partition_of_n_;
    int     rand_cell;
    int     corner_i, corner_j;
    init_tables_for_random_SYT_generation();

    do{
        //rand_cell = rand() % n_count;
        rand_cell = (int) (unif_rand() * n_count);
        walk(rand_cell, &corner_i, &corner_j);
        syt_[corner_i][ corner_j ] = n_count;
        n_count --;
        //if just one box in last row, delete row
        if (corner_i + 1 == ferrers_shape_length_dynamic_ && corner_j == 0) ferrers_shape_length_dynamic_--;
        //delete box
        ferrers_shape_dynamic_[corner_i ] --;
        //SYT_toStr();
    }while (n_count > 0);
    
}

void Ferrers_diagram::walk(int rand_cell, int *cell_i, int* cell_j){
    //get i j for the random cell
    int c_i = 0 , c_j = 0  ; // cell position
    int hook_right, hook_bottom; //to count the number of cells to the right and bottom
    int cont = 0;
    int ran ;
    bool enc = false, is_dynamic_corner;
    for (int i = 0 ; (!enc && i < ferrers_shape_length_dynamic_ ) ; i++)
        if (rand_cell  < ferrers_shape_dynamic_[ i ]){
            c_i = i ;
            c_j = rand_cell;
            enc = true;
        }else rand_cell -=  ferrers_shape_dynamic_[ i ];

    is_dynamic_corner  = ferrers_shape_dynamic_[ c_i ] ==  ( c_j + 1)
    && (c_i == ferrers_shape_length_dynamic_ - 1 || ferrers_shape_dynamic_[ c_i + 1 ] <=  c_j );
    while ( !is_dynamic_corner ){
        hook_right = ferrers_shape_dynamic_[ c_i ] - c_j - 1;
        cont = c_i + 1;
        hook_bottom = 0 ;
        while (cont < ferrers_shape_length_dynamic_ && c_j < ferrers_shape_dynamic_[ cont++ ] )
            hook_bottom ++;
        //ran = rand () % (hook_right + hook_bottom ) ;
        ran = (int) (unif_rand() * (hook_right + hook_bottom ) ) ;
        if ( ran < hook_bottom ) c_i += ran + 1;
        else c_j += ran - hook_bottom + 1;
        is_dynamic_corner  = ferrers_shape_dynamic_[ c_i ] ==  ( c_j + 1)
            && (c_i == ferrers_shape_length_dynamic_ - 1 || ferrers_shape_dynamic_[ c_i + 1 ] <=  c_j );
    }
    *cell_i = c_i;
    *cell_j = c_j;
    
}

long double Ferrers_diagram::calculate_hook_length(long double n_factorial){
//    Generic gen;
    if( hook_length_ == -1 ){
        hook_length_ = 1 ;
        long double cell_hook;
        for (int i = 0 ; i < ferrers_shape_length_; i ++){
            for (int j = 0 ; j < ferrers_shape_[ i ]; j++) {
                cell_hook = ferrers_shape_[ i ] - j + get_num_cells_down(i, j) - 1;
                hook_length_ = hook_length_ * cell_hook;
            }
        }
        hook_length_ = n_factorial/ hook_length_;
        number_of_permutations_ = hook_length_ * hook_length_;
    }
    return hook_length_;
}

int Ferrers_diagram::get_num_cells_down(int i, int j){
    //returns only the ones below AND itself
    int num_cols_down = 0 ;
    while ( i < ferrers_shape_length_ && ferrers_shape_[ i ] > j ){
        num_cols_down++;
        i++;
    }
    return num_cols_down;
}





















