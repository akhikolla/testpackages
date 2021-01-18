/* This file is part of SMITIDvisu package.
* Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>

* SMITIDvisu is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SMITIDvisu is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with SMITIDvisu. If not, see <https://www.gnu.org/licenses/>.
*/


#include <vector>
#include <utility>
#include <algorithm>

using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

#ifndef __KRUSKAL_HPP__
#define __KRUSKAL_HPP__

//' @title compute the minimum spanning tree
//' @name mstCompute
//' @description compute the minimum spanning tree of a matrix representing edges between nodes (of a graph)
//' @param mat weighted matrix representing nodes connection (edges weight)
//' @return a matrix with 1 if nodes are linked, 0 otherwise.
// [[Rcpp::export]]
NumericMatrix mstCompute(const NumericMatrix &mat);

double kruskal(vector< pair <long double, pair<int, int> > > edges, NumericMatrix &out);

bool edgesCompare(const pair <long double, pair <int, int> > &i, const pair <long double, pair <int, int> > &j);

int rootTree( int * tree, int node );

void addTree( int * tree, int node1, int node2);


#endif
