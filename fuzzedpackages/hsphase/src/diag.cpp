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
 * diag.cpp
 *
 *  Created on: 17/01/2013
 *  Author: Mohammad H. Ferdosi
 */

#include "diag.h"
#include "block.h"
void diagnostic(int const * matrix, int *nrow, int *ncol, int * result,
		vector<int> &hetSite, int* zeroFrq, int* twoFrq)
{
	for (int i = 0; i < (*nrow) * (*ncol); i++)
	{
		result[i] = matrix[i];
	}

	//For matrix
	int const**pRowsMat = new int const*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pRowsMat[i] = matrix + (*nrow) * i;
	}

	// For result
	int **pRowsRes = new int*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pRowsRes[i] = result + (*nrow) * i;
	}

	int *tempCol1 = new int[*nrow];
	int *tempCol2 = new int[*nrow];
	int *tempCol2R = new int[*nrow];
	int *keeptempCol2 = new int[*nrow];
	vector<int> switchA, switchB;
	memoryCLS tbackwardMemory(nrow);

	for (int i = 0; i < *nrow; i++)
	{
		tempCol1[i] = pRowsMat[hetSite[0]][i];
	}

	strandOrigin(tempCol1, nrow);
	for (int j = 0; j < *nrow; j++)
		pRowsRes[hetSite[0]][j] = tempCol1[j];

	// Initialise backward memory
	tbackwardMemory.memoryMaker(tempCol1);

	// Main loop
	int n = 1;
	vector<int>::iterator iElement;

	for (iElement = hetSite.begin() + 1; iElement != hetSite.end(); ++iElement)
	{

		double zero =
				(*(twoFrq + *iElement) < *(zeroFrq + *iElement)) ?
						*(zeroFrq + *iElement) : *(twoFrq + *iElement);
		double two =
				(*(twoFrq + *iElement) > *(zeroFrq + *iElement)) ?
						*(zeroFrq + *iElement) : *(twoFrq + *iElement);

		if (((double) (two / zero)) < .4
				|| ((double) (zero / (*nrow)) > .8))
		{
			continue;
		}

		for (int i = 0; i < *nrow; i++)
		{
			*(keeptempCol2 + i) = *(*(pRowsMat + *iElement) + i);
			*(tempCol2 + i) = *(*(pRowsMat + *iElement) + i);
			*(tempCol2R + i) = *(*(pRowsMat + *iElement) + i);
		}

		strandOrigin(tempCol2, nrow);
		strandOrigin(tempCol2R, nrow);
		reverseConvert(tempCol2R, nrow);
		switchDetector(tbackwardMemory.pMemory, tempCol2, switchA, nrow);
		switchDetector(tbackwardMemory.pMemory, tempCol2R, switchB, nrow);

		if (switchA.size() < switchB.size())
		{
			for (int j = 0; j < (*nrow); j++)
			{
				if (tempCol2[j] == 3 || tempCol2[j] == 4)
				{
					*(*(pRowsRes + *iElement) + j) = *(tempCol2 + j);
				}
			}
			tbackwardMemory.memoryMaker(tempCol2);
		}
		else if (switchA.size() > switchB.size())
		{
			for (int j = 0; j < (*nrow); j++)
			{
				if (tempCol2[j] == 3 || tempCol2[j] == 4)
				{
					*(*(pRowsRes + *iElement) + j) = *(tempCol2R + j);
				}
			}
			tbackwardMemory.memoryMaker(tempCol2R);
		}

		n = n + 1;
		switchA.clear();
		switchB.clear();
	}

	delete[] pRowsMat;
	delete[] pRowsRes;
	delete[] tempCol1;
	delete[] tempCol2;
	delete[] tempCol2R;
	delete[] keeptempCol2;

}
