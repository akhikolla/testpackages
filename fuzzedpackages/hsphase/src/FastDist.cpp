// Copyright (C) 2014 Mohammad H. Ferdosi
//
// HSPhase is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// HSPhase program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


/*
 * FastDist.cpp
 *
 *  Created on: 09/09/2013
 *      Author: mhf
 */

#include "FastDist.h"
using namespace std;

int dist(int const * matrix, const int* nrow, const int* ncol, int* result)
{
	//For matrix

	int const**pRowsMat = new int const*[*nrow];

	for (int i = 0; i < *nrow; i++)
	{
		pRowsMat[i] = matrix + (*ncol) * i;
	}

	// For result
	int **pRowsRes = new int*[*nrow];

	for (int i = 0; i < *nrow; i++)
	{
		pRowsRes[i] = result + (*nrow) * i;
	}

	int frq = 0;
#pragma omp parallel for private(frq)  schedule(dynamic)  num_threads(2)
	for (int i = 0; i < *nrow; i++)
	{
		for (int j = i; j < *nrow; j++)
		{
			for (int k = 0; k < *ncol; k++)
			{
				/*if ((pRowsMat[i][k] == 2 && pRowsMat[j][k] == 0) || (pRowsMat[i][k] == 0 && pRowsMat[j][k] == 2))
				 {
				 frq = frq + 1;
				 }*/
				frq = frq + abs(pRowsMat[i][k] - pRowsMat[j][k]);
			}
			pRowsRes[i][j] = frq;
			frq = 0;
		}
	}


	delete[] pRowsMat;
	delete[] pRowsRes;
	return 0;
}

extern "C"
{

void fastDist(int const * matrix, int *nrow, int *ncol, int* result)
{

	dist(matrix, nrow, ncol, result);

}
}
