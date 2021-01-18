// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include <Rcpp.h>
#include "Joint_Grid.h"
#include "Clusters.h"

// [[Rcpp::export]]
Rcpp::List findgrid(Rcpp::List cluster_info, int k, int nobs, int ndims, int bin_limit) {

  // obtain cluster centers
  // Rcpp::NumericMatrix centers = Rcpp::as<Rcpp::NumericMatrix>(cluster_info["centers"]);

  // obtain cluster labels
  Rcpp::IntegerVector labels = Rcpp::as<Rcpp::IntegerVector>(cluster_info["clusters"]);

  // obtain data
  Rcpp::NumericMatrix data = Rcpp::as<Rcpp::NumericMatrix>(cluster_info["data"]);

  // type-cast centers to std::vector<std::vector<double > >
  // vector<vector<double> > _centers(k, vector<double>(ndims, 0));
  // for(int i=0; i<k; i++){
  //   for(int j=0; j<ndims; j++){
  //     _centers[i][j] = centers(i,j);
  //   }
  // }

  // type-cast labels to std::vector<int>
  vector<int> _labels = Rcpp::as<vector<int> >(labels);

  // type-cast data to std::vector<std::vector<double > >
  vector<vector<double> > _data(nobs, vector<double>(ndims,0));
  for(int i=0; i<nobs; i++){
    for(int j=0; j<ndims; j++){
      _data[i][j] = data(i,j);
    }
  }

  // make a 'cluster' object
  Cluster clust_obj(k, _labels, _data);

  // get grid lines
  vector<vector<double> > _grid_lines = Find_Grid(clust_obj, bin_limit);

  // type cast it to List
  Rcpp::List lines(ndims);
  for(int i=0; i<ndims; i++){
    lines[i] = Rcpp::wrap(_grid_lines[i]);
  }

  return lines;
}
