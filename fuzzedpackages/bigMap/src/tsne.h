/*
 The bigMap Package for R.

 Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

 bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
*/

#ifndef TSNE_H
#define TSNE_H

static inline int tIdx(int n, int i, int j) {
	return ( (n*i) - (i+1)*(i+2)/2 + j );
}

class TSNE
{
public:
	// transform input similarities into P from INPUT-DATA
	void X2P(double* X, int n, int m, double* Beta, double* P);
	// transform input similarities into P from DISTANCES Triangular-Matrix
	void D2P(double* X, int n, double* Beta, double* P);
	// run TSNE
	void run2D(int n, double* P, double* Y, double &Cost, double alpha, int max_iter);
private:
	void Gradient(double* Y, int n, int dimY, double* P, double &Q, double* atrF, double* repF);
	void getCost(double* Y, int n, int dimY, double* P, double &Q, double &Cost);
};

#endif
