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
 * File:   block.h
 * Author: Mohammad H. Ferdosi
 *
 * Created on 21 October 2012, 2:17 PM
 */

#ifndef BLOCK_H
#define	BLOCK_H
#include <vector>
#include "hsp.h"

using namespace std;

class block
{
public:

	block(int const **pRowsMat, int ** result, int * tbackwardmemory,
			vector<int> switchLocus, vector<int>::iterator iElement, int *nrow,
			vector<int> &hetsite);
	virtual ~block();

	//methods
	int recombinationDetector(int ur, int ul, int dr, int dl);
	int validitySwitch(int windowsWidth);
	int makeMemory(int *nsap);

	int* getItsBackwardMemoryMain_() const
	{
		return (itsBackwardMemoryMain_);
	}

	int* getCrossover() const
	{
		return crossover;
	}

	void setCrossover(int* crossover)
	{
		this->crossover = crossover;
	}

private:
	const int *itsNrow_;
	vector<int> itsHetsite_;
//	const int *itsWindowsWidth_;
	int *itsValidateMem_;
	int *itsValidateConst_;
	int *itsBackwardMemoryMain_;
	int *itstbackwardmemoryBlock_;
	int const**itspRowsMat_;
	int** itsResult_;
	int* crossover;

	vector<int>::iterator itsiElement_;
	vector<int>::iterator itsHetEnd_;

};
int c2rphaseOPT(int const * matrix, int const *zeroFrq, int const * oneFrq,
		int const* twoFrq, int* nrow, int* ncol, int* result,
		vector<int> &hetSite, int *forwardvectorsize, bool *excludeFP,int *nsap);
int c2rphaseOPTRecombination(int const * matrix, int const *zeroFrq, int const * oneFrq, int const* twoFrq, int* nrow,
		int* ncol, int* result, vector<int> &hetSite, int *forwardvectorsize /* = 0 */, bool *excludeFP, int *nsap,
		int* recombinationMatrix);

#endif	/* BLOCK_H */

