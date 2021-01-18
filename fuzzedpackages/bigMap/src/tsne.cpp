/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <Rcpp.h>

#include "tsne.h"

// ony for X2P dgemm_ version !!
// extern "C" {
//     #include <R_ext/BLAS.h>
// }

using namespace std;

// Perform t-SNE (only for 2D embedding !!)
void TSNE::run2D(int n, double* P, double* Y, double &Cost, double alpha, int max_iter) {

	// Only tested for 2D embeddings !!!!
	int dimY = 2;

	// +++ Parameters:
	// 1. momentum; alpha = alpha
	// 2. learning rate; eta = etaX * range(Y)/2 * log(n-1) * 4
	//    etaX : learning rate boost factor (dropped)
	//    log(n-1) : sample size scaling of eta
	//    *4.0 comes from the gradient solution

	std::vector<double> Y_range (dimY*2, 0);
	for (int d = 0; d < dimY; d++){
		// set minimum
		Y_range[d*dimY + 0] = -1.0;
		// check maximum
		Y_range[d*dimY + 1] = 1.0;
	}

	std::vector<double> eta (dimY, 0);
	double etaX = std::log(n-1) * 2.0;

	// . input similarities (P distribution)
	// P[i, j] == 0,  means j not in neighbourhood of i

	// . output similarities (Q distribution)
	// . Q normalization factor
	double Q = .0;

	// +++ Allocate memory

	// . gradient forces
	double* atrF = (double*) malloc(n * dimY * sizeof(double));
	if (atrF == NULL) Rcpp::stop("Memory allocation failed! \n");
	double* repF = (double*) malloc(n * dimY * sizeof(double));
	if (repF == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . update of mapped positions
	double* uY = (double*) calloc(n * dimY, sizeof(double));
	if (uY== NULL) Rcpp::stop("Memory allocation failed! \n");

	// +++ Embedding optimization
	for (int iter = 0; iter < max_iter; iter++) {
		// Compute gradient
		Gradient(Y, n, dimY, P, Q, atrF, repF);
		// set current value for learning rate
		for (int d=0; d<dimY; d++){
			eta[d] = (Y_range[d *dimY +1] - Y_range[d *dimY +0]) * etaX;
		}
		// update (with momentum and learning rate)
		for (int i = 0; i < n; i++){
			for (int d = 0; d < dimY; d++){
				int k = i*dimY + d;
				uY[k] = alpha * uY[k] - eta[d] * (atrF[k] - repF[k] / Q);
				Y[k] += uY[k];
				// check minimum
				if (Y[k] < Y_range[d*dimY + 0]) Y_range[d*dimY + 0] = Y[k];
				// check maximum
				else if (Y[k] > Y_range[d*dimY + 1]) Y_range[d*dimY + 1] = Y[k];
			}
		}
	}

	// +++ cost function
	getCost(Y, n, dimY, P, Q, Cost);

	// +++ Clean up memory
	free(atrF); atrF  = NULL;
	free(repF); repF = NULL;
	free(uY); uY = NULL;
}

// Compute gradient forces of the t-SNE cost function
void TSNE::Gradient(double* Y, int n, int dimY, double* P, double &Q, double* atrF, double* repF)
{
	// Set gradient forces to zero
	for(int i = 0; i < n; i++) {
		for(int d = 0; d < dimY; d++){
			atrF[i*dimY + d] = .0;
			repF[i*dimY + d] = .0;
		}
	}
	// Set Q normalization factor to zero
	Q = FLT_MIN;
	// Compute new gradient forces
	std::vector<double> L (dimY, 0);
	for(int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++){
			double Lij = FLT_MIN;
			for(int d = 0; d < dimY; d++){
				L[d] = Y[i*dimY + d] - Y[j*dimY + d];
				Lij +=  L[d] * L[d];
			}
			double Qij = 1.0 / (1.0 + Lij);
			Q += Qij;
			int ij = tIdx(n, i, j);
			for(int d = 0; d < dimY; d++) {
				// attractive gradient
				atrF[i*dimY + d] += P[ij] * Qij * L[d];
				atrF[j*dimY + d] -= P[ij] * Qij * L[d];
				// repulsive gradient
				repF[i*dimY + d] += Qij * Qij * L[d];
				repF[j*dimY + d] -= Qij * Qij * L[d];
			}
		}
	}
	Q *= 2.0;
}

// Evaluate t-SNE cost function (exactly)
void TSNE::getCost(double* Y, int n, int dimY, double* P, double &Q, double &Cost)
{
	Q = .0;
	double jxH = .0;	// join cross entropy
	for (int i = 0; i < n; i++) {
		double Si = .0;
		for (int d = 0; d < dimY; d++) Si += Y[i *dimY +d] * Y[i *dimY +d];
		for (int j = i+1; j < n; j++){
			double Sj = .0;
			double Lij = .0;
			for(int d = 0; d < dimY; d++){
				Sj += Y[j *dimY +d] * Y[j *dimY +d];
				Lij -= Y[i*dimY + d] * Y[j*dimY + d];
			}
			Lij *= 2.0;
			Lij += (Si + Sj + FLT_MIN);
			double Qij = 1.0 / (1.0 + Lij);
			int ij = tIdx(n, i, j);
			jxH += P[ij] * std::log(Qij);
			Q += Qij;
		}
	}
	// normalize cross entropy
	jxH *= 2.0;
	jxH += std::log(2.0 *Q);
	// cost
	Cost = jxH /std::log(n *(n -1));
}

void TSNE::X2P(double* X, int n, int m, double* Beta, double* P)
{
	// allocate memory
	// . squared components
	double* S = (double*) malloc(n * sizeof(double));
	if (S == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . row similarity distribution
	double* M = (double*) malloc(n * sizeof(double));
	if (M == NULL) Rcpp::stop("Memory allocation failed! \n");
	// compute squared components
	for (int i = 0; i < n; i++) {
		S[i] = .0;
		for(int v = 0; v < m; v++) S[i] += X[i *m +v] * X[i *m +v];
	}
	// compute similarities & marginals
	for (int i = 0; i < n; i++) {
		double Zi = .0;
		for (int j = 0; j < i; j++){
			double Lji = .0;
			for(int v = 0; v < m; v++) Lji -= X[i *m +v] * X[j *m +v];
			Lji *= 2.0;
			Lji += S[i] + S[j] + FLT_MIN;
			M[j] = std::exp(-Beta[j] * Lji);
			Zi += M[j];
		}
		for (int j = i+1; j < n; j++) {
			double Lij = .0;
			for(int v = 0; v < m; v++) Lij -= X[i *m +v] * X[j *m +v];
			Lij *= 2.0;
			Lij += S[i] + S[j] + FLT_MIN;
			M[j] = std::exp(-Beta[i] * Lij);
			Zi += M[j];
		}
		for (int j = 0; j < i; j++) {
			int ji = tIdx(n, j, i);
			P[ji] += M[j] /Zi /(2 *n);
		}
		for (int j = i+1; j < n; j++) {
			int ij = tIdx(n, i, j);
			P[ij] += M[j] /Zi /(2 *n);
		}
	}
	// deallocate memory
	free(S); S = NULL;
	free(M); M = NULL;
}

// // transform input similarities into probabilities
// // FROM INPUT DATA
// void TSNE::X2P(double* X, int n, int m, double* Beta, double* P)
// {
// 	// allocate memory
// 	// . distances
// 	double* L = (double*) malloc(n * (n-1) / 2 * sizeof(double));
// 	if (L == NULL) Rcpp::stop("Memory allocation failed! \n");
// 	// . row marginals
// 	double* Z = (double*) malloc(n * sizeof(double));
// 	if (Z == NULL) Rcpp::stop("Memory allocation failed! \n");
// 	for (int i = 0; i < n; i++) Z[i] = FLT_MIN;
// 	// compute similarities & marginals
// 	for (int i = 0; i < n; i++) {
// 		double Si = .0;
// 		for(int v = 0; v < m; v++) Si += X[i *m +v] * X[i *m +v];
// 		for (int j = i+1; j < n; j++) {
// 			double Sj = .0;
// 			int ij = tIdx(n, i, j);
// 			L[ij] = .0;
// 			for(int v = 0; v < m; v++){
// 				Sj += X[j *m +v] * X[j *m +v];
// 				L[ij] -= X[i *m +v] * X[j *m +v];
// 			}
// 			L[ij] *= 2.0;
// 			L[ij] += (Si + Sj + FLT_MIN);
// 			Z[i] += std::exp(-Beta[i] * L[ij]);
// 			Z[j] += std::exp(-Beta[j] * L[ij]);
// 		}
// 	}
// 	// normalization & symmetrization
// 	for (int i = 0; i < n; i++) {
// 		for (int j = i+1; j < n; j++) {
// 			int ij = tIdx(n, i, j);
// 			P[ij] = (std::exp(-Beta[i] * L[ij]) / Z[i] + std::exp(-Beta[j] * L[ij]) / Z[j]) / (2*n);
// 		}
// 	}
// 	// deallocate memory
// 	free(L); L = NULL;
// 	free(Z); Z = NULL;
// }

// transform input similarities into probabilities
// FROM FULL-DISTANCE-MATRIX
void TSNE::D2P(double* D, int n, double* Beta, double* P)
{
	// allocate memory
	// . row marginals
	double* Z = (double*) malloc(n * sizeof(double));
	if (Z == NULL) Rcpp::stop("Memory allocation failed! \n");
	for (int i = 0; i < n; i++) Z[i] = FLT_MIN;
	// compute similarities & marginals
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			int ij = tIdx(n, i, j);
			double Lij = D[i*n + j] * D[i*n + j];
			Z[i] += std::exp(-Beta[i] * Lij);
			Z[j] += std::exp(-Beta[j] * Lij);
		}
	}
	// normalization & symmetrization
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			int ij = tIdx(n, i, j);
			double Lij = D[i*n + j] * D[i*n + j];
			P[ij] = (std::exp(-Beta[i] * Lij) / Z[i] + std::exp(-Beta[j] * Lij) / Z[j]) / (2*n);
		}
	}
	// deallocate memory
	free(Z); Z = NULL;
}
