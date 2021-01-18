/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// +++ grid functions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
arma::Mat<double> grid_init(arma::Col<double> X, arma::Col<double> Y){

	// kernel density class matrix
	arma::Mat<double> grid(2, 6);

	// size, min, max, range, and factor=size/range of X
	grid(0, 0) = X.size();
	grid(0, 1) = std::floor(arma::min(X));
	grid(0, 2) = std::ceil(arma::max(X));
	grid(0, 3) = grid(0, 2) - grid(0, 1);
	grid(0, 4) = grid(0, 0) / grid(0, 3);
	// size, min, max, range, and factor=size/range of Y
	grid(1, 0) = Y.size();
	grid(1, 1) = std::floor(arma::min(Y));
	grid(1, 2) = std::ceil(arma::max(Y));
	grid(1, 3) = grid(1, 2) - grid(1, 1);
	grid(1, 4) = grid(1, 0) / grid(1, 3);

	// cell size
	grid(0, 5) = std::sqrt(std::pow((X[1]-X[0]), 2) + std::pow((Y[1]-Y[0]), 2));
	grid(1, 5) = grid(0, 5);

	return grid;
}

// [[Rcpp::export]]
arma::Col<int> grid_p2cell(double x, double y, arma::Mat<double> grid){
	// Convert embedded point xy-coordinates to grid cell coordinates
	arma::Col<int> p2c(2);
	p2c[0] = (int) std::floor((x - grid(0, 1)) * grid(0, 4));
	p2c[1] = (int) std::floor((y - grid(1, 1)) * grid(1, 4));
	return p2c;
}

// [[Rcpp::export]]
arma::Mat<double> grid_D2cell(arma::Mat<double> D, arma::Mat<double> grid){
	// Convert a vector of embedded data-points to grid cell coordinates
	arma::Mat<double> D2c(D.n_rows,2);
	for (size_t i = 0; i < D.n_rows; ++i){
		arma::Col<int> p2c = grid_p2cell(D(i, 0), D(i, 1), grid);
		D2c(i, 0) = p2c[0];
		D2c(i, 1) = p2c[1];
	}
	return D2c;
}

// [[Rcpp::export]]
arma::Col<int> grid_n2cell(int n, arma::Mat<double> grid){

	// convert node (grid ordinal number) to grid cell coordinates
	// given x=X.size() and y=Y.size(),
	// conversion goes from left to right and from bottom to top:

	//		(y-1)x	(y-1)x+1	(y-1)x+2	...		yx-1
	//
	//		2x		2x+1		2x+2		...		3x-1
	//		x		x+1			x+2			...		2x-1
	//		0		1			2			...		x-1

	arma::Col<int> n2c(2);
	n2c[1] = (int) std::floor(n / grid(0, 0));
	n2c[0] = (int) std::floor(n - (n2c[1] * grid(0, 0)));

	return n2c;
}

// [[Rcpp::export]]
arma::Mat<int> grid_N2cell(arma::Mat<double> grid){
	// convert all grid nodes to grid cell coordinates
	int grid_size = grid(0, 0) * grid(1, 0);
	arma::Mat<int> N2c(grid_size, 2);
	for (int n = 0; n < grid_size; ++n){
		N2c(n, 1) = (int) std::floor(n / grid(0, 0));
		N2c(n, 0) = (int) std::floor(n - (N2c(n, 1) * grid(0, 0)));
	}
	return N2c;
}

// [[Rcpp::export]]
arma::Mat<int> grid_M2cell(arma::Col<int> M, arma::Mat<double> grid){
	// convert vector of nodes to grid cell coordinates
	arma::Mat<int> M2c(M.size(), 2);
	for (size_t i = 0; i < M.size(); ++i){
		arma::Col<int> n2c = grid_n2cell(M[i], grid);
		M2c(i, 0) = n2c[0];
		M2c(i, 1) = n2c[1];
	}
	return M2c;
}

// [[Rcpp::export]]
arma::Col<int> grid_bound(int n, arma::Mat<double> grid){

	// get boundaring nodes of node n

	arma::Col<int> oBound(9);
	oBound.zeros();
	int oBound_size = 0;

	// 8-connectivity (8 Chamfer neighbourhood on a square grid)
	arma::Col<int> iBound(9);
	iBound[0] = n -grid(0, 0) -1;
	iBound[1] = n -grid(0, 0);
	iBound[2] = n -grid(0, 0) +1;
	iBound[3] = n -1;
	iBound[4] = n;
	iBound[5] = n +1;
	iBound[6] = n +grid(0, 0) -1;
	iBound[7] = n +grid(0, 0);
	iBound[8] = n +grid(0, 0) +1;

	// remove not valid neighbours
	arma::Col<int> n2c = grid_n2cell(n, grid);
	for (size_t i = 0; i < iBound.size(); ++i){
		arma::Col<int> b2c = grid_n2cell(iBound[i], grid);
		// neighbours must be contiguous
		if (std::abs(n2c[0] - b2c[0]) > 1) continue;
		if (std::abs(n2c[1] - b2c[1]) > 1) continue;
		// neighbour cell coordinates must lay in the grid limits
		if (b2c[0]>=0 && b2c[1]>=0 && b2c[0]<grid(0,0) && b2c[1]<grid(1,0)){
			oBound[oBound_size] = iBound[i];
			++oBound_size;
		}
	}

	if (oBound_size > 0) oBound.resize(oBound_size);
	return oBound;
}

// [[Rcpp::export]]
arma::Col<int> grid_cross(int n, arma::Mat<double> grid){

	// get cross nodes of node n

	arma::Col<int> oBound(5);
	oBound.zeros();
	int oBound_size = 0;

	// 4-connectivity (4 Chamfer neighbourhood on a square grid)
	arma::Col<int> iBound(5);
	iBound[0] = n -grid(0, 0);
	iBound[1] = n -1;
	iBound[2] = n;
	iBound[3] = n +1;
	iBound[4] = n +grid(0, 0);

	arma::Col<int> n2c = grid_n2cell(n, grid);
	for (size_t i = 0; i < iBound.size(); ++i){
		arma::Col<int> b2c = grid_n2cell(iBound[i], grid);
		// neighbours must be contiguous
		if (std::abs(n2c[0] - b2c[0]) > 1) continue;
		if (std::abs(n2c[1] - b2c[1]) > 1) continue;
		// neighbour cell coordinates must lay in the grid limits
		if (b2c[0]>=0 && b2c[1]>=0 && b2c[0]<grid(0,0) && b2c[1]<grid(1,0)){
			oBound[oBound_size] = iBound[i];
			++oBound_size;
		}
	}

	if (oBound_size > 0) oBound.resize(oBound_size);
	return oBound;
}

// [[Rcpp::export]]
arma::Col<int> grid_peaks(arma::Mat<double> Z, arma::Mat<double> grid){

	// arma::vec X = arma::sort(vectorise(Z));
	arma::Col<double> V = arma::vectorise(Z);
	arma::Col<int> M(V.size());

	int M_size = 0;
	for (size_t p = 0; p < V.size(); ++p){
		arma::Col<int> p_bound = grid_bound(p, grid);
		size_t i = 0;
		while (i < p_bound.size() && V[p] >= V[p_bound[i]]) ++i;
		if (i == p_bound.size()){
			M[M_size] = p;
			++M_size;
		}
	}

	if (M_size > 0) M.resize(M_size);

	return M;
}


// +++ DHWS functions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
Rcpp::List wtt_cpp(arma::Col<double> X, arma::Col<double> Y, arma::Mat<double> Z){

	arma::Mat<double> grid = grid_init(X, Y);
	arma::Col<double> S = arma::vectorise(Z);
	arma::Col<uword> J = arma::sort_index(S);

	// initialize peaks vector
	arma::Col<int> M(S.size());

	// pre-assign first cluster centroid to the highest peak
	M[0] = J[(J.size() -1)];
	int M_size = 0;

	// initialize cell cluster assignments
	arma::Col<int> K(S.size());
	K = -K.ones();

	for (int j = (J.size()-1); j >= 0 ; --j){

		int n = J[j];
		arma::Col<int> n_bound = grid_bound(n, grid);

		// consider the node as the highest of its neighbourhood
		// (pre-assign node to itself, as a new cluster)
		int highest = n;
		for (size_t i = 0; i < n_bound.size(); ++i){
			if (K[n_bound[i]] >= 0 && S[n_bound[i]] >= S[highest]){
			 	highest = n_bound[i];
				K[n] = K[highest];
			}
		}

		// if K[n]==-2 it is the first (representative) node of a regional maxima
		if (K[n] == -2){
			K[n] = M_size;
			M[M_size] = n;
			++M_size;
		}
		// if K[n]==-1 this node has no other labelled node in its neighborhood,
		// thus either:
		// 1. it is a peak (the centroid of a new cluster);
		if (K[n] == -1){
			// pre-assign it to itself as a new cluster
			K[n] = M_size;
			// check it is the highest of its neighbourhood;
			for (size_t i = 0; i < n_bound.size(); ++i){
				if (n_bound[i] != n && S[n_bound[i]] >= S[n]) K[n] = -1;
			}
			// add it to the list of peaks (clusters)
			if (K[n] != -1){
				M[M_size] = n;
				++M_size;
			}
		}
		// 2. it belongs to a plateau and either:
		// 2.1 it is a set of regional maxima then there is no labelled node in the extended neighbourhood of n;
		// 2.2 it is a set of regional minima then there must be some labelled node in the extended neighborhood of n;
		// in any case, move it to the last position with respect to other nodes with same S[n];
		if (K[n] == -1){
			// eventually mark it as representative
			K[n] = -2;
			int u = j-1;
			while (u > -1 && S[n] == S[J[u]]){
				// check no other node has already been marked as representative
				if (K[J[u]] != -1) K[n] = -1;
				J[u+1] = J[u];
				J[u] = n;
				--u;
			}
			// set pointer 1 position back and retry
			++j;
		}

	}

	if (M_size>0) M.resize(M_size);
	arma::Mat<int> M2c = grid_M2cell(M, grid);

	return Rcpp::List::create(
		Named("grid")=grid,
		Named("P")=M,
		Named("M")=M2c,
		Named("C")=K );

}
