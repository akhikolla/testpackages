//
// Created by Jiandong Wang on 2/14/20.
//
// Copyright (c) NMSU Song lab


#ifndef JOINT_GRID_JOINT_GRID_H
#define JOINT_GRID_JOINT_GRID_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <climits>

#include "Clusters.h"

using namespace std;

vector< vector<double > > Find_Grid(Cluster &clusters, int min_bin_limit = 2);

vector<double> Find_1D_Grid(Cluster &clusters, int dim, int min_bin_limit);

vector<vector<double> > prep_index(vector<double> &c1, vector<double> &c2, double middle_1, double middle_2);


double binary_search_index(const vector<vector<double> > &c_index, const int left, const int right, const int size_c1,
                        const int size_c2, bool &overlap, vector<double> &err_sum);

#endif //JOINT_GRID_JOINT_GRID_H
