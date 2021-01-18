/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>

#include <math.h>
#include <float.h>
#include "tsne.h"

using namespace Rcpp;

// +++ SCKT -------------------------------------------------------------------

// Exact implementation of t-SNE (SCKT)
// [[Rcpp::export]]
double sckt_zTSNE(int thread_rank, int threads, int layers, SEXP X, SEXP B, SEXP Y, SEXP I, double iters, double alpha, bool isDistance)
{
	// input data
	Rcpp::XPtr<BigMatrix> bmX(X);
	MatrixAccessor<double> maX(*bmX);
	int n = bmX -> nrow();
	int m = bmX -> ncol();
	// input Betas
	Rcpp::XPtr<BigMatrix> bmB(B);
	MatrixAccessor<double> maB(*bmB);
	// current mapping positions
	Rcpp::XPtr<BigMatrix> bmY(Y);
	MatrixAccessor<double> maY(*bmY);
	// sampled row indexes
	Rcpp::XPtr<BigMatrix> bmI(I);
	MatrixAccessor<int> maI(*bmI);

	// get chunk breaks
	std::vector<int> breaks (threads+1, 0);
	for (int b = 0; b < threads; b++) breaks[b] = (int) b *(n+1.0)/threads;
	breaks[threads] = n;

	// get thread-size
	int thread_size = 0;
	for (int l = 0; l < layers; l++)
	{
		int a = (thread_rank + l) % threads;
		thread_size += (breaks[a+1] - breaks[a]);
		//printf(" rank %3d/%3d , chunk %3d : %6d ... %6d \n", thread_rank, threads, l, maI[0][breaks[a]], maI[0][(breaks[a+1]-1)]);
	}

	// get thread-indexes & thread-layers
	std::vector<int> zIdx (thread_size, 0);
	std::vector<int> zLay (thread_size, 0);
	int z = 0;
	for (int l = 0; l < layers; l++)
	{
		int a = (thread_rank + l) % threads;
		for (int i = breaks[a]; i < breaks[a+1]; i++){
			zIdx[z] = maI[0][i];
			zLay[z] = l;
			z++;
		}
	}

	// . chunk of input data
	if (isDistance) m = thread_size;
	double* thread_X = (double*) malloc(thread_size *m *sizeof(double));
	if (thread_X == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . chunk of Betas
	double* thread_B = (double*) malloc(thread_size *1 *sizeof(double));
	if (thread_B == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . starting mapping positions
	double* thread_Y = (double*) malloc(thread_size *2 *sizeof(double));
	if (thread_Y == NULL) Rcpp::stop("Memory allocation failed! \n");

	if (isDistance)
	{
		for (int i = 0; i < thread_size; i++){
			for (int j = 0; j < thread_size; j++){
				thread_X[i *thread_size + j] = maX[zIdx[j]][zIdx[i]];
			}
			for (int d = 0; d < 2; d++){
				thread_Y[i *2 +d] = maY[zLay[i] *2 +d][zIdx[i]];
			}
			thread_B[i] = maB[0][zIdx[i]];
			z++;
		}

	}
	else
	{
		for (int i = 0; i < thread_size; i++){
			for (int j = 0; j < m; j++){
				thread_X[i *m + j] = maX[j][zIdx[i]];
			}
			for (int d = 0; d < 2; d++){
				thread_Y[i *2 +d] = maY[zLay[i] *2 +d][zIdx[i]];
			}
			thread_B[i] = maB[0][zIdx[i]];
			z++;
		}
	}

	// . input similarities (P distribution)
	// P[i, j] == 0,  means j not in neighbourhood of i
	double* thread_P = (double*) calloc(thread_size * (thread_size-1) / 2, sizeof(double));
	if (thread_P == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . cost function value
	double thread_Cost = .0;

	// +++ TSNE instance
	TSNE* tsne = new TSNE();

	// +++ Compute input similarities
	if (isDistance){
		tsne -> D2P(thread_X, thread_size, thread_B, thread_P);
	} else {
		tsne -> X2P(thread_X, thread_size, m, thread_B, thread_P);
	}

	// +++ Run tsne
	tsne -> run2D(thread_size, thread_P, thread_Y, thread_Cost, alpha, iters);

	// update mapping positions
	z = 0;
	for (int l = 0; l < layers; l++)
	{
		int a = (thread_rank + l) % threads;
		for (int b = 0; b < (breaks[a+1] -breaks[a]); b++)
		{
			int j = maI[0][breaks[a] +b];
			for (int d = 0; d < 2; d++) maY[l *2 +d][j] = thread_Y[z *2 +d];
			z++;
		}
	}

	// deallocate memory
	delete(tsne);
	free(thread_X); thread_X = NULL;
	free(thread_B); thread_B = NULL;
	free(thread_Y); thread_Y = NULL;
	free(thread_P); thread_P = NULL;

	return thread_Cost;

}

// +++ MPI --------------------------------------------------------------------

// get chunks by thread: indexes and mapping-positions
// [[Rcpp::export]]
void zChnks(Rcpp::List& Z_list, const arma::Mat<double>& Y, const arma::Col<int>& I, const Rcpp::List& brks_list)
{
	for (int z = 0; z < brks_list.length(); z++)
	{
		arma::Mat<int> brks = brks_list[z];
		arma::Mat<double> zChnk = Z_list[z];
		size_t j = 0;
		for (size_t l = 0; l < brks.n_rows; l++) {
			for (size_t i = brks(l, 0); i < brks(l, 1); i++) {
				zChnk(j, 0) = I[i];
				zChnk(j, 1) = Y(I[i], 2*l + 0);
				zChnk(j, 2) = Y(I[i], 2*l + 1);
				j++;
			}
		}
		Z_list[z] = zChnk;
	}
}

// restructure global mapping
// [[Rcpp::export]]
void updateY(arma::Mat<double>& Y, const arma::Col<int>& I, const Rcpp::List& zMap_list, const Rcpp::List& brks_list)
{
	for (int z = 0; z < zMap_list.length(); z++)
	{
		arma::Mat<int> thrd_brks = brks_list[z];
		arma::Mat<double> zY = zMap_list[z];
		size_t k = 0;
		for (size_t l = 0; l < thrd_brks.n_rows; l++){
			int j1 = l *2, j2 = l *2 +1;
			for (int i = thrd_brks(l, 0); i < thrd_brks(l, 1); i++)
			{
				Y(I[i], j1) = zY(k, 0);
				Y(I[i], j2) = zY(k, 1);
				++k;
			}
		}
	}
}

// compute embedding size
// [[Rcpp::export]]
arma::Col<double> eSize(arma::Mat<double>& Y)
{
	arma::Mat<double> col_range(2, Y.n_cols);
	col_range.zeros();
	for (size_t i = 0; i < Y.n_rows; i++){
		for (size_t j = 0; j < Y.n_cols; j++)
		{
			if (Y(i, j) < col_range(0, j)) col_range(0, j) = Y(i, j);
			else if (Y(i, j) > col_range(1, j)) col_range(1, j) = Y(i, j);
		}
	}
	arma::Col<double> lyr_size(Y.n_cols /2);
	lyr_size.zeros();
	for (size_t l = 0; l < lyr_size.size(); l++)
	{
		lyr_size[l] += std::pow((col_range(1, 2 *l) - col_range(0, 2 *l)), 2);
		lyr_size[l] += std::pow((col_range(1, 2 *l +1) - col_range(0, 2 *l +1)), 2);
		lyr_size[l] = std::sqrt(lyr_size[l]);
	}
	return lyr_size;
}

// Exact implementation of t-SNE (MPI)
// [[Rcpp::export]]
double mpi_zTSNE(SEXP X, SEXP B, arma::Mat<double>& Y, const arma::Col<int>& I, double iters, double alpha, bool isDistance)
{
	// input data
	Rcpp::XPtr<BigMatrix> bmX(X);
	MatrixAccessor<double> maX(*bmX);
	// input Betas
	Rcpp::XPtr<BigMatrix> bmB(B);
	MatrixAccessor<double> maB(*bmB);

	// input dimensions
	int thread_size = Y.n_rows;
	int m = bmX -> ncol();

	// chunk of data
	double* thread_X = (double*) malloc(thread_size *m *sizeof(double));
	if (thread_X == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . chunk of Betas
	double* thread_B = (double*) malloc(thread_size *1 *sizeof(double));
	if (thread_B == NULL) Rcpp::stop("Memory allocation failed! \n");
	// // . starting mapping positions
	double* thread_Y = (double*) malloc(thread_size *2 *sizeof(double));
	if (thread_Y == NULL) Rcpp::stop("Memory allocation failed! \n");

	if (isDistance)
	{
		for (int i = 0; i < thread_size; i++) {
			for (int d = 0; d < 2; d++) {
				thread_Y[i *2 +d] = Y(i, d);
			}
			for (int j = 0; j < thread_size; j++) {
				thread_X[i *thread_size + j] = maX[I[j]][I[i]];
			}
			thread_B[i] = maB[0][I[i]];
		}
	}
	else
	{
		for (int i = 0; i < thread_size; i++) {
			for (int d = 0; d < 2; d++) {
				thread_Y[i *2 +d] = Y(i, d);
			}
			for (int j = 0; j < m; j++) {
				thread_X[i *m + j] = maX[j][I[i]];
			}
			thread_B[i] = maB[0][I[i]];
		}
	}

	// . input similarities (P distribution)
	// P[i, j] == 0,  means j not in neighbourhood of i
	double* thread_P = (double*) calloc(thread_size * (thread_size-1) / 2, sizeof(double));
	if (thread_P == NULL) Rcpp::stop("Memory allocation failed! \n");
	// . cost function value
	double thread_Cost = .0;

	// +++ TSNE instance
	TSNE* tsne = new TSNE();

	// +++ Compute input similarities
	if (isDistance) {
		tsne -> D2P(thread_X, thread_size, thread_B, thread_P);
	} else {
		tsne -> X2P(thread_X, thread_size, m, thread_B, thread_P);
	}

	// +++ Run tsne
	tsne -> run2D(thread_size, thread_P, thread_Y, thread_Cost, alpha, iters);

	// update mapping positions
	for (int i = 0; i < thread_size; i++) {
		for (int d = 0; d < 2; d++) {
			Y(i, d) = thread_Y[i *2 +d];
		}
	}

	// deallocate memory
	delete(tsne);
	free(thread_X); thread_X = NULL;
	free(thread_B); thread_B = NULL;
	free(thread_Y); thread_Y = NULL;
	free(thread_P); thread_P = NULL;

	return thread_Cost;

}
