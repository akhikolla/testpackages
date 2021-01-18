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
 * File:   hsp.cpp
 * Author: Mohammad H. Ferdosi
 *
 * Created on 13 October 2012, 6:14 PM
 */

#include "memory.h"
#include "hsp.h"

using namespace std;

/**
 *
 * @param col
 * @param nrow
 * @return convert 4 to 3 and 3 to 4 (3 and 4 are strand A and B of sire)
 */
int reverseConvert(int *col, int *nrow)
{
	for (int i = 0; i < *nrow; i++)
	{
		if (col[i] == 4)
			col[i] = 5;
		if (col[i] == 3)
			col[i] = 4;
		if (col[i] == 5)
			col[i] = 3;

	}
	return (0);
}

int reverseConvert(int *col, const int *nrow)
{
	for (int i = 0; i < *nrow; i++)
	{
		if (col[i] == 4)
			col[i] = 5;
		if (col[i] == 3)
			col[i] = 4;
		if (col[i] == 5)
			col[i] = 3;

	}
	return (0);
}

/**
 *
 * @param matrix
 * @param nrow
 * @param ncol
 * @param result
 * @param MaxBlock
 * @return convert 3 to 1 and 4 to 2 and the rest to 0 (1 and 2 are strand A and B of sire)
 */
int c2rBlocks(int const * matrix, int *nrow, int *ncol, int* result,
		int* MaxBlock)
{

	c2rBlocks2(matrix, nrow, ncol, result);

	for (int i = 0; i < (*nrow) * (*ncol); i++)
	{
		if (result[i] == 3)
			result[i] = 1;
		else if (result[i] == 4)
			result[i] = 2;
		else
			result[i] = 0;
	}

	return (0);

}

int c2rBlocks2(int const* matrix, int *nrow, int const*ncol, int* result)
{
	int fAnchor = 0, loc1 = 0, loc2 = 0;

	for (int i = 0; i < (*nrow) * (*ncol); i++)
	{
		*(result + i) = *(matrix + i);
	}

	for (int i = 0; i < (*nrow) * (*ncol); i = i + *nrow)
	{
		loc1 = loc2 = i;
		fAnchor = 0;

		for (int j = i; j < (*nrow + i); j++)
		{

			if (*(matrix + j) == 3 || *(matrix + j) == 4)
			{
				if (fAnchor == 0)
				{
					fAnchor = *(matrix + j); //fAnchor = 3 or 4
					loc1 = j;
					for (int k = i; k < loc1; k++)
					{
						*(result + k) = *(matrix + j);
					}

				}
				else
				{
					if (fAnchor == *(matrix + j))
					{
						loc2 = j;
						for (int k = loc1; k < loc2; k++)
						{
							*(result + k) = *(matrix + j);
						}
						loc1 = j;
					}
					else
					{
						fAnchor = *(matrix + j);
						loc1 = j;
					}
				}
			}
		}
		for (int z = (*nrow + i) - 1; z != i; z--)
		{

			if (*(matrix + z) == 3 || *(matrix + z) == 4)
			{

				int temp = *(matrix + z);
				for (int k = z; k < (*nrow + i); k++)
				{
					*(result + k) = temp;
				}
				break;
			}
		}
	}

	for (int i = 0; i < (*nrow) * (*ncol); i++)
	{
		if (*(result + i) != 3 && *(result + i) != 4)
			*(result + i) = 0;
	}

	return (0);
}

int memMaker(int *lastMemory, int* newCol, int *nrow)
{
	for (int i = 0; i < *nrow; i++)
	{
		if (*(newCol + i) == 3 || *(newCol + i) == 4)
		{
			*(lastMemory + i) = *(newCol + i);
		}
	}
	return (0);

}


int strandOrigin(int *col, int *nrow)
{
	for (int i = 0; i < *nrow; i++)
	{
		if (*(col + i) == 0)
			*(col + i) = 3;
		else if (*(col + i) == 2)
			*(col + i) = 4;
	}
	return (0);

}

int strandOrigin(int *col, const int *nrow)
{
	for (int i = 0; i < *nrow; i++)
	{
		if (*(col + i) == 0)
			*(col + i) = 3;
		else if (*(col + i) == 2)
			*(col + i) = 4;
	}
	return (0);

}

int switchDetector(int *Memory, int *tempCol2, int *nrow)
{
	int sw = 0;
	for (int i = 0; i < *nrow; i++)
	{
		if ((*(Memory + i) == 3 || *(Memory + i) == 4)
				&& (*(tempCol2 + i) == 3 || *(tempCol2 + i) == 4))
		{
			if (*(Memory + i) != *(tempCol2 + i))
				sw = sw + 1;
		}
	}

	return (sw);
}

int switchDetector(int *Memory, int *tempCol2, const int *nrow)
{
	int sw = 0;
	for (int i = 0; i < *nrow; i++)
	{
		if ((*(Memory + i) == 3 || *(Memory + i) == 4)
				&& (*(tempCol2 + i) == 3 || *(tempCol2 + i) == 4))
		{
			if (*(Memory + i) != *(tempCol2 + i))
				sw = sw + 1;
		}
	}

	return (sw);
}


int switchDetector(int *Memory, int *tempCol2, vector<int> &switches, int *nrow)
{
	switches.clear();
	switches.reserve(*nrow); // Initialise the variable

	for (int i = 0; i < *nrow; i++)
	{
		if ((*(Memory + i) == 3 || *(Memory + i) == 4)
				&& (*(tempCol2 + i) == 3 || *(tempCol2 + i) == 4))
		{
			if (*(Memory + i) != *(tempCol2 + i))
			{
				switches.push_back(i);
			}
		}
	}

	return (0);
}

int switchDetector(int *Memory, int *tempCol2, vector<int> &switches,
		const int *nrow)
{
	switches.clear();
	switches.reserve(*nrow);
	for (int i = 0; i < *nrow; i++)
	{
		if ((*(Memory + i) == 3 || *(Memory + i) == 4)
				&& (*(tempCol2 + i) == 3 || *(tempCol2 + i) == 4))
		{
			if (*(Memory + i) != *(tempCol2 + i))
			{
				switches.push_back(i);
			}
		}
	}

	return (0);
}
int c2rStrandF(const int* matBlock, const int* matGenotype, const int* nrow,
		const int* ncol, double* result)
{

	//

	for (int i = 0; i < 2 * (*ncol); i++)
	{
		result[i] = 9;
	}
	//For Block
	int const**pRowsBolck = new int const*[*ncol];
	for (int i = 0; i < *ncol; i++)
	{
		pRowsBolck[i] = matBlock + (*nrow) * i;
	}

	//For Genotype
	int const**pRowsGen = new int const*[*ncol];
	for (int i = 0; i < *ncol; i++)
	{
		pRowsGen[i] = matGenotype + (*nrow) * i;
	}

	//For result
	double **pRowsRes = new double*[2];
	for (int i = 0; i < 2; i++)
	{
		pRowsRes[i] = result + (*ncol) * i;
	}

	double tempSumStrand1 = 9, tempSumStrand2 = 9;
	vector<int> firstAllele;

	for (int k = 0; k < (*ncol); k++)
	{
		for (int j = 0; j < *nrow; j++)
		{
			if (*(*(pRowsBolck + k) + j) != 0)
			{
				//Save the ind with 0 and 1
				if (*(*(pRowsGen + k) + j) == 0 || *(*(pRowsGen + k) + j) == 2)
				{
					firstAllele.push_back(j);
				}
			}
		}

		vector<int>::iterator iElement;
		double block11 = 0, block12 = 0, block21 = 0, block22 = 0;
		for (iElement = firstAllele.begin(); iElement != firstAllele.end();
				++iElement)
		{
			if (*(*(pRowsBolck + k) + *iElement) == 1
					&& *(*(pRowsGen + k) + *iElement) != 1)
			{
				if (*(*(pRowsGen + k) + *iElement) == 0)
					block11 = block11 + 1;
				if (*(*(pRowsGen + k) + *iElement) == 2)
					block12 = block12 + 1;
			}
			if (*(*(pRowsBolck + k) + *iElement) == 2
					&& *(*(pRowsGen + k) + *iElement) != 1)
			{
				if (*(*(pRowsGen + k) + *iElement) == 0)
					block21 = block21 + 1;
				if (*(*(pRowsGen + k) + *iElement) == 2)
					block22 = block22 + 1;
			}
		}

		tempSumStrand1 = floor(
				((double) ((block12) / (block11 + block12))) + .5);
		tempSumStrand2 = floor(
				((double) ((block22) / (block21 + block22))) + .5);

		if ((block11 + block12) == 0 || (block21 + block22) == 0)
		{
			tempSumStrand1 = 9;
			tempSumStrand2 = 9;
		}

		pRowsRes[0][k] = (double) tempSumStrand1;
		pRowsRes[1][k] = (double) tempSumStrand2;

		tempSumStrand1 = 9;
		tempSumStrand2 = 9;
		block11 = 0, block12 = 0, block21 = 0, block22 = 0;
		firstAllele.clear();
	}

	delete[] pRowsBolck;
	delete[] pRowsGen;
	delete[] pRowsRes;

	return (0);
}

int frequencyVector(vector<int>::iterator strat, int seachNumber,
		vector<int>::iterator end)
{
	int frq = 0;
	vector<int>::iterator iElement;
	for (iElement = strat; iElement != end; ++iElement)
	{

		if (*iElement == seachNumber)
			frq = frq + 1;
	}
	return (frq);
}

int phaseFunction(int const * genotypeMat, int const *nrow, int const *ncol,
		int const* blockMat, int const* sirePhasedMat, int* result)
{
	for (int i = 0; i < (*ncol) * (*nrow); i++)
	{
		result[i] = 9;
	}

	//for genotype
	int const**pGenotypeMat = new int const*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pGenotypeMat[i] = genotypeMat + (*nrow) * i;
	}

	int const**pBlockMat = new int const*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pBlockMat[i] = blockMat + (*nrow) * i;
	}

	int const**pSirePhasedMat = new int const*[2];

	for (int i = 0; i < 2; i++)
	{
		pSirePhasedMat[i] = sirePhasedMat + (*ncol) * i;
	}

	int **pResMat = new int *[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pResMat[i] = result + (*nrow) * i;
	}

	for (int i = 0; i < *ncol; i++)
	{
		for (int j = 0; j < *nrow; j++)
		{
			if (*(*(pBlockMat + i) + j) == 1)
			{
				*(*(pResMat + i) + j) = *(*(pSirePhasedMat) + i);
			}
			if (*(*(pBlockMat + i) + j) == 2)
			{
				*(*(pResMat + i) + j) = *(*(pSirePhasedMat + 1) + i);
			}

		}
	}

	for (int i = 0; i < *ncol; i++)
	{
		for (int j = 0; j < *nrow; j++)
		{
			if (*(*(pGenotypeMat + i) + j) == 0)
				*(*(pResMat + i) + j) = 0;
			if (*(*(pGenotypeMat + i) + j) == 2)
				*(*(pResMat + i) + j) = 1;
		}

	}

	delete[] pGenotypeMat;
	delete[] pBlockMat;
	delete[] pSirePhasedMat;
	delete[] pResMat;

	return (0);
}

int  phaseFunctionNoGenotype(int const *nrow, int const *ncol,
		int const* blockMat, int const* sirePhasedMat, int* result)
{
	for (int i = 0; i < (*ncol) * (*nrow); i++)
	{
		result[i] = 9;
	}


	int const**pBlockMat = new int const*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pBlockMat[i] = blockMat + (*nrow) * i;
	}

	int const**pSirePhasedMat = new int const*[2];

	for (int i = 0; i < 2; i++)
	{
		pSirePhasedMat[i] = sirePhasedMat + (*ncol) * i;
	}

	int **pResMat = new int *[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pResMat[i] = result + (*nrow) * i;
	}

	for (int i = 0; i < *ncol; i++)
	{
		for (int j = 0; j < *nrow; j++)
		{
			if (*(*(pBlockMat + i) + j) == 1)
			{
				*(*(pResMat + i) + j) = *(*(pSirePhasedMat) + i);
			}
			if (*(*(pBlockMat + i) + j) == 2)
			{
				*(*(pResMat + i) + j) = *(*(pSirePhasedMat + 1) + i);
			}

		}
	}


	delete[] pBlockMat;
	delete[] pSirePhasedMat;
	delete[] pResMat;

	return (0);
}
int c2rRecombinations(const uint* matrix, const uint* nrow, const uint* ncol, const uint* method, double* result)
{
	uint fAnchor = 0, loc1 = 0, loc2 = 0;
	double ratio = 0;

	for (uint i = 0; i < (*nrow) * (*ncol); i++)
	{
		result[i] = 0;
	}

	for (uint i = 0; i < (*nrow) * (*ncol); i = i + *nrow)
	{
		loc1 = loc2 = i;
		fAnchor = 0;
		ratio = 0;

		for (uint j = i; j < *nrow + i; j++)
		{

			if (matrix[j] == 1 || matrix[j] == 2)
			{
				if (fAnchor == 0)
				{
					fAnchor = matrix[j];
					loc1 = j;
				}
				else
				{

					if (fAnchor != matrix[j])
					{
						loc2 = j;
						fAnchor = matrix[j];
					}
					if (*method == 1)
						ratio = 1;
					else if (*method == 2)
						ratio = 1 / (double) (loc2 - loc1);
					for (uint k = loc1; k < loc2; k++)
					{
						result[k] = ratio;
					}
					loc1 = j;
				}
			}
		}

	}

	return (0);
}
int recombinationFun(int const *matrix, int const * nrow, int const * ncol,
		int * result)
{
	for (int i = 0; i < (*nrow); i++)
	{
		result[i] = 9;
	}


	//for genotype
	int const**blockMat = new int const*[*nrow];

	for (int i = 0; i < *nrow; i++)
	{
		blockMat[i] = matrix + (*ncol) * i;
	}

	int n = 0;

	for (int i = 0; i < *nrow; i++)
	{
		for (int j = 0; j < *ncol - 1; j++)
		{
			if (blockMat[i][j] != blockMat[i][j+1] && blockMat[i][j] != 0)
			{
				n = n + 1;
			}
		}
		result[i] = n;

		n = 0;
	}

	return (0);
}

