//
//  OptFramedClust.h
//  OptClust_Log_Linear
//
//  Created by Tathagata Debnath on 6/15/20.
//  Copyright © 2020 Tathagata Debnath and Joe Song.
//  All rights reserved.
//  Revised by Joe Song.

// #include <stdio.h>
// #include <cassert>
// #include <cstring>
// #include <iostream>
// #include <string>
// #include <list>

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

using namespace std;

struct frame_info
{
  double ssq = std::numeric_limits<double>::infinity(); // 1.79769e+308;
  int Frame_ID = -1;
};

struct clustering
{
  int Frame_ID = -1;
  double ssq = std::numeric_limits<double>::infinity(); // 1.79769e+308;
  std::vector<int> Borders;
  std::vector<double> centers;
  std::vector<double> size;
  double totss = std::numeric_limits<double>::infinity(); //1.79769e+308;
  std::vector<double> withinss;


};

clustering MFC(
    std::vector<double> & Data_Points,
    int width, int K,
    int First, int Last,
    int Prev, int Next);

frame_info BDP(
    int width, int K,
    int First, int Last,
    int Prev, int Next,
    std::vector< std::vector< double > > &  S,
    std::vector< std::vector< int > > & J,
    const std::vector<double> & sum_x,
    const std::vector<double> & sum_x_sq,
    std::vector< std::vector<int> > & Cluster_Border);


void linear_clustering(
    std::vector< std::vector< double > > & S,
    std::vector< std::vector< int > > & J,
    int Prev, int Next, int Middle_Frame,
    const std::vector<double> & sum_x,
    const std::vector<double> & sum_x_sq,
    std::vector< std::vector<int> > & Cluster_Border
);


void fill_row_k(
    int imin, int imax, int k, int Middle_Frame,
    int jmin, int jmax,
    std::vector< std::vector<double> > & S,
    std::vector< std::vector<int> > & J,
    const std::vector<double> & sum_x,
    const std::vector<double> & sum_x_sq);


/*
 void fill_row_q_2020_07_22(int imin, int imax, int q, int imin_0, int imax_0, int Middle_Frame,
 int jmin, int jmax,
 std::vector< std::vector<double> > & S,
 std::vector< std::vector<int> > & J,
 const std::vector<double> & sum_x,
 const std::vector<double> & sum_x_sq);
 */


double ssq(
    const int j, const int i, const int Middle_Frame,
    const std::vector<double> & sum_x, // running sum of xi
    const std::vector<double> & sum_x_sq // running sum of xi^2
);




void backtrack(
    std::vector< std::vector<int> > & J,
    std::vector<int> & B, int K, int N
);






