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
 * Ohd.cpp
 *
 *  Created on: 02/08/2013
 *      Author: mhf
 */

#include "Ohd.h"
#define mat imat
#define uint unsigned int

using namespace Rcpp;
using namespace std;
/**
 *
 * @param genotpe
 * @param elementSearch
 * @param atLeast
 * @return
 */
arma::vec frequency(arma::mat genotpe, int elementSearch, bool atLeast)
{
	arma::vec hetsite(genotpe.n_cols);
	hetsite.zeros();

	for (uint i = 0; i < genotpe.n_cols; i++)
	{
		for (uint j = 0; j < genotpe.n_rows; j++)
		{
			if (genotpe(j, i) == elementSearch)
			{
				hetsite[i] = hetsite[i] + 1;
				if (atLeast)
					break;
			}
		}
	}

	return (hetsite);
}

/**
 *
 * @param zeroFrq
 * @param twoFrq
 * @param total
 * @param genotype
 * @return
 */
arma::vec hetIndDetector(arma::vec zeroFrq, arma::vec twoFrq, arma::vec total, arma::mat genotype, bool uniq_check)
{
	arma::vec ohFrq(genotype.n_rows);
	ohFrq.zeros();
	for (uint i = 0; i < genotype.n_cols; i++)
	{
		if (total[i] == 2)
		{
			if (twoFrq[i] == 1 || uniq_check)
			{
				for (uint j = 0; j < genotype.n_rows; j++)
				{
					if (genotype(j, i) == 2)
					{
						ohFrq[j] = ohFrq[j] + 1;
					}
				}

			}
			if (zeroFrq[i] == 1 || uniq_check)
			{
				for (uint j = 0; j < genotype.n_rows; j++)
				{
					if (genotype(j, i) == 0)
					{
						ohFrq[j] = ohFrq[j] + 1;
					}
				}
			}

		}
	}
	return (ohFrq);
}
/**
 *
 * @param geno
 * @return
 */
SEXP ohd(SEXP geno, SEXP unique)
{
	NumericMatrix genotype(geno);
	int n = genotype.nrow(), k = genotype.ncol();
	arma::mat Genotype(as<arma::mat>(genotype).begin(), n, k, TRUE);

	bool u_check = as<bool> (unique);

	arma::vec zeroFreq = frequency(Genotype, 0, TRUE);
	arma::vec twoFreq = frequency(Genotype, 2, TRUE);
	arma::vec total = zeroFreq + twoFreq;

	zeroFreq = frequency(Genotype, 0, FALSE);
	twoFreq = frequency(Genotype, 2, FALSE);

	arma::vec oh = hetIndDetector(zeroFreq, twoFreq, total, Genotype, u_check);

	return wrap(oh);
}

