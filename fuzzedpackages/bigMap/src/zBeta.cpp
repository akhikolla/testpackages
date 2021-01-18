/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <bigmemory/MatrixAccessor.hpp>

using namespace arma;
using namespace Rcpp;
using namespace std;

// point-wise squared distances from/to k
// [[Rcpp::export]]
arma::Col<double> distk(int k, SEXP X, bool is_distance)
{
	// input data
	Rcpp::XPtr<BigMatrix> ptrX(X);
	MatrixAccessor<double> mtxX(*ptrX);
	int n = ptrX -> nrow();
	int m = ptrX -> ncol();

	// output: vector of distances
	arma::Col<double> Dk(n);
	Dk.zeros();

	// input is distance.matrix
	if (is_distance) {
		for (int i = 0; i < n; ++i) Dk[i] = mtxX[i][k] * mtxX[i][k];
	}
	// input is data.matrix
	else {
		double Sk = .0;
		for(int d = 0; d < m; d++) Sk += (mtxX[d][k] * mtxX[d][k]);
		for (int i = 0; i < n; ++i) {
			double Si = .0;
			for(int d = 0; d < m; d++) {
				Si += (mtxX[d][i] * mtxX[d][i]);
				Dk[i] -= (mtxX[d][k] * mtxX[d][i]);
			}
			Dk[i] *= 2.0;
			Dk[i] += (Sk + Si + FLT_MIN);
		}
	}

	return Dk;
}

// compute chunk of Betas
// [[Rcpp::export]]
arma::Col<double> zBeta(int thread_rank, int threads, SEXP X, bool is_distance, double ppx, double tol, int mxI)
{
	// printf(" rank %3d/%3d \n", thread_rank, threads);

	// input data
	Rcpp::XPtr<BigMatrix> ptrX(X);
	MatrixAccessor<double> mtxX(*ptrX);
	int n = ptrX -> nrow();

	// get chunk breaks
	std::vector<int> breaks (threads+1, 0);
	for (int b = 0; b < threads; b++) breaks[b] = (int) b *(n +1.0) /threads;
	breaks[threads] = n;

	// printf(" rank %3d/%3d , chunk: %6d ... %6d \n", thread_rank, threads, breaks[thread_rank], breaks[thread_rank+1]-1);

	// output data (chunk Betas)
	int chunk_size = breaks[thread_rank] - breaks[thread_rank-1];
	arma::Col<double> Beta(chunk_size);
	Beta.ones();

	// internal parameters
	double minCk = n *DBL_EPSILON;
	double logppx = std::log(ppx);

	int k = 0;
	for (int i = breaks[thread_rank-1]; i < breaks[thread_rank]; i++){

		// +++ get point-wise squared distances to/from i
		arma::Col<double> Dk = distk(i, X, is_distance);

		double maxBeta = DBL_MAX, minBeta = 0.0, Hdff = 0.0;
		// +++ compute Beta[k]
		for (int itr = 0; itr < mxI; ++itr) {

			// increase/decrease Beta[k]
			if (Hdff > 0.0) {

				minBeta = Beta[k];
				// on first iterations avoid exploring whole range up to DBL_MAX
				if (maxBeta == DBL_MAX || maxBeta == 0.0) Beta[k] *= 2.0;
				else Beta[k] = (Beta[k] + maxBeta) / 2.0;
			}
			else if (Hdff < 0.0){
				maxBeta = Beta[k];
				// on first iterations avoid exploring whole range down to 0.0
				if (minBeta == 0.0 || minBeta == DBL_MAX) Beta[k] /= 2.0;
				else Beta[k] = (Beta[k] + minBeta) / 2.0;
			}

			// as Dk[k] = 0, Pk[k] is going to be 1
			// but we don't want to take this into account in Ck
			double Ck = -1.0, Hk =  0.0;
			for (int j = 0; j < n; ++j){
				double Pj = std::exp(- Beta[k] * Dk[j]);
				Hk += Dk[j] * Pj;
				Ck += Pj;
			}

			// check that not all density is at point k
			if (Ck > minCk){
				Hdff = std::log(Ck) + Beta[k] * Hk/Ck - logppx;
			}
			// if it is, decrease betak
			else Hdff = -1.0;

			// check stop criterium
			if (std::abs(Hdff) < tol) break;
		}

		// +++ next k
		k++;
	}
	return Beta;
}
