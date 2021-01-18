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
 * File:   memory.cpp
 * Author: Mohammad H. Ferdosi
 *
 * Created on 27 October 2012, 3:25 PM
 */

#include "memory.h"

memoryCLS::memoryCLS(int* nrow)
{
	pMemory = new int[*nrow];
	pPreviousMemory = new int[*nrow];
	pValidity = new int[*nrow];
	pRows = nrow;

	fill_n(pMemory, *nrow, 9);
	fill_n(pPreviousMemory, *nrow, 0);
	fill_n(pValidity, *nrow, 0);
}

memoryCLS::~memoryCLS()
{
	delete[] pMemory;
	delete[] pPreviousMemory;
	delete[] pValidity;
}

int memoryCLS::freeMemory()
{
	delete[] pMemory;
	delete[] pPreviousMemory;
	delete[] pValidity;
	return (0);
}

int memoryCLS::memoryMaker(int* col)
{
	for (int i = 0; i < *pRows; i++)
	{
		if (col[i] == 3 || col[i] == 4)
		{
			pMemory[i] = col[i];
		}
	}
	return (0);
}

int memoryCLS::memoryMaker(int* col, int *validation)
{
	for (int i = 0; i < *pRows; i++)
	{
		if ((col[i] == 3 || col[i] == 4) && validation[i] == 1)
		{
			pMemory[i] = col[i];
		}
	}
	return (0);
}

int memoryCLS::evalution(int *swichLocations)
{
	for (int i = 0; i < *pRows; i++)
	{
		if (swichLocations[i] != 1)
		{
			pValidity[i] = pValidity[i] + 1;
		}
		else
		{
			pValidity[i] = 0;
		}
	}
	return (0);
}
