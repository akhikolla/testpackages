/*
 * cluster.h
 *
 *  Created on: Sep 18, 2015
 *      Author: root
 */

#ifndef SRC_CLUSTER_H_
#define SRC_CLUSTER_H_

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

List clusterl0(vector<mat> data, uvec svs, int sz, int seed, int maxIter, double thres, int debug);
List clusterl1(vector<mat> data, fvec lus, fvec lvs, float lz, int maxOuter, double thresOuter, int maxInner, double thresInner, int debug);
List cMvSsvd(vector<mat> data, fvec lvs, float lz, int maxOuter, double thresOuter, int maxInner, double thresInner, int debug);

#endif /* SRC_CLUSTER_H_ */
