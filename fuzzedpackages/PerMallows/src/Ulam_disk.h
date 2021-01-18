//
//  Ulam_disk.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 31/01/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Ulam_disk__
#define __perms_mallows__Ulam_disk__

#define GEN_NEXT  0 //ok, print and continue
#define GEN_TERM  1 //ok, terminate
#define GEN_EMPTY 2 //ok, print EMPTY SET and continue
#define GEN_ERROR 3 //an error occured, print an error message and terminate

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include "Ulam.h"
#include "Newton_raphson.h"
#include "Ferrers_diagram.h"
#include "Exponential_model.h"


class Ulam_disk;
using namespace std;
class Ulam_disk: public Ulam{
protected:

    void    read_permus_per_dist();

    void    read_mutiple_shapes_from_file_bin(int d, double * bounds, int num_bounds, int ** shapes, int*shape_len);
    
    void    distances_sampling(int m, double theta, int**samples);
    void    generate_permu_with_given_LIS (int l, int *sigma, int* shape, int shape_len);
    void    estimate_theta(int m, int *sigma_0, int **samples, int model, double *theta);
    long double get_likelihood(int m, int** samples, int model, int * sigma_0);
    
    double  right_ ;
    int     left_, aux_,cont_;
    ///
    int     * line_block_;
    int     lines_in_block_;



    
public:
    //int n_;
    char str_base_path [500] ;//= "/Users/eki/Desktop/aux_files/";
    Ulam_disk(int n) : Ulam(n){
        strcpy(str_base_path, "./");///Users/eki/Library/Developer/Xcode/DerivedData/perms_mallows-fxsohiuulqqerjgnrdpitiwnjkin/Build/Products/Debug/permus_per_shape_5_0
        lines_in_block_ = 128;
        line_block_ = new int[lines_in_block_ * (n_ + 3 )]; for (int i = 0 ; i < lines_in_block_ * (n_ + 3 ) ; i ++ ) line_block_[ i ] = 0 ;
    }

    ~Ulam_disk(){
        delete [] line_block_;
    };
    
    void    save_counts_to_file_bin();
    double  expectation(double theta);

    /*************          ABSTRACT        ***************/
    //// local
    long double num_permus_at_distance(int d) {
        read_permus_per_dist();
        return num_permus_per_dist_[ d ];
    }
  
    /*************          ABSTRACT   NOT SUPPORTED     ***************/




};



#endif /* defined(__perms_mallows__Ulam_disk__) */


















