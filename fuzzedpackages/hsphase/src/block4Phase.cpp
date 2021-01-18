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
 * File:   block4Phase.cpp
 * Author: Mohammad H. Ferdosi
 *
 * Created on 27 October 2012, 3:34 PM
 */

#include "block4Phase.h"

block4Phase::block4Phase(const uint *matrix, const uint *nrow, const uint *ncol,
uint *result, const uint *siregenotype, const uint * str)
{
	SNP *halfsibStrands = new SNP[*nrow]; // real row
	SNP sire;
	str_ = *str;

	for (uint i = 0; i < *ncol; i++)
	{
		sire.strand1.push_back(siregenotype[i]);
		sire.strand2.push_back(siregenotype[i + (*ncol)]);

	}

	for (uint j = 0; j < ((*nrow) / 2); j++)
	{
		for (uint i = j * ((*ncol) * 2); i < (j * ((*ncol) * 2) + (*ncol)); i++)
		{
			halfsibStrands[j].strand1.push_back(matrix[i]);
			halfsibStrands[j].strand2.push_back(matrix[i + (*ncol)]);
		}
	}

	for (uint i = 0; i < (*ncol) * (*nrow) / 2; i++)
	{
		result[i] = matrix[i];
	}

	int *block = new int[*ncol];

	for (uint j = 0; j < (*nrow) / 2; j++)
	{

		blockMaker(sire, halfsibStrands[j], block, ncol);

		for (uint i = j * ((*ncol)); i < (j * ((*ncol)) + (*ncol)); i++)
		{
			result[i] = block[i - (j * ((*ncol)))];
		}

		sire.strand1.clear();
		sire.strand2.clear();

		for (uint i = 0; i < *ncol; i++)
		{
			sire.strand1.push_back(siregenotype[i]);
			sire.strand2.push_back(siregenotype[i + (*ncol)]);

		}

	}

}

block4Phase::block4Phase(const block4Phase& orig)
{
}

block4Phase::~block4Phase()
{
}

int block4Phase::sireStrdDetector(const SNP &sire, const SNP &halfsib)
{
    int std1 = 0, std2 = 0;
	for (uint i = 0; i < sire.strand1.size(); i++)
	{
		if (sire.strand1[i] == sire.strand2[i] && halfsib.strand1[i] != halfsib.strand2[i])
		{
			if (sire.strand1[i] != halfsib.strand1[i])
				std1 = std1 + 1;
			if (sire.strand1[i] != halfsib.strand2[i])
				std2 = std2 + 1;
		}
	}

	if (std1 < std2)
		return 1;
	else
		return 2;
	return (0);
}

int block4Phase::blockMaker(SNP& sire, const SNP& halfsib, int * block /* = 0 */, const uint* ncol /* = 0 */)
{
	int strand = 0;
	if (str_ == 0)
		strand = this->sireStrdDetector(sire, halfsib);
	if (str_ == 1)
		strand = 1;
	if (str_ == 2)
		strand = 2;
	if (strand == 0)
	{
		for (uint i = 0; i < *ncol; i++)
		{
			block[i] = 0;
		}
	}
	else if (strand == 1)
	{
		int n = 3, t = 5;
		for (uint i = 0; i < halfsib.strand1.size(); i++)
		{
			if (sire.strand1[i] == sire.strand2[i] && sire.strand1[i] != 9)
			{
				if (t == 0)
					block[i] = n;
				else
					block[i] = 0;

			}
			else if (halfsib.strand1[i] == sire.strand1[i]
					&& halfsib.strand1[i] != 9)
			{
				block[i] = n;
				t = 0;
			}
			else if (halfsib.strand1[i] == sire.strand2[i]
					&& halfsib.strand1[i] != 9)
			{

				if (n == 3)
					n = 4;
				else
					n = 3;
				block[i] = n;
				sire.recombination(i);
				t = 0;
			}
			else
			{
				block[i] = 0;
			}
		}

	}
	else if (strand == 2)
	{
		int n = 3, t = 5;
		for (uint i = 0; i < halfsib.strand1.size(); i++)
		{
			if (sire.strand1[i] == sire.strand2[i] && sire.strand1[i] != 9)
			{
				if (t == 0)
					block[i] = n;
				else
					block[i] = 0;
			}
			else if (halfsib.strand2[i] == sire.strand1[i]
					&& halfsib.strand1[i] != 9)
			{
				block[i] = n;
				t = 0;
			}
			else if (halfsib.strand2[i] == sire.strand2[i]
					&& halfsib.strand1[i] != 9)
			{
				if (n == 3)
					n = 4;
				else
					n = 3;
				block[i] = n;
				sire.recombination(i);
				t = 0;
			}
			else
			{
				block[i] = 0;
			}
		}

	}
	else
	{
		return (1);
	}

	return (0);
}

SNP::SNP()
{

}

SNP::~SNP()
{

}

SNP SNP::recombination(unsigned int index)
{
	if (index > strand1.size() || index > strand2.size())
	{
		index = strand1.size();
	}
	if (index != 0)
	{
		vector<uint> tempStrand1;
		vector<uint> tempStrand2;
		for (uint i = 0; i < index; i++)
		{
			tempStrand1.push_back(strand1[i]);
			tempStrand2.push_back(strand2[i]);
		}

		for (unsigned int i = index; i < strand1.size(); i++)
		{
			tempStrand1.push_back(strand2[i]);
			tempStrand2.push_back(strand1[i]);
		}

		strand1 = tempStrand1;
		strand2 = tempStrand2;
	}
	else
	{
		vector<uint> temp;
		temp = strand1;
		strand1 = strand2;
		strand2 = temp;
	}
	return (*this);
}


