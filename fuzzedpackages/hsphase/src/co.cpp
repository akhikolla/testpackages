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
 * co.cpp
 *
 *  Created on: 03/08/2013
 *      Author: mhf
 */

#include "co.h"

#define mat imat
#define uint unsigned int

using namespace Rcpp;
using namespace std;

/**
 *
 * @param vector1 first individual
 * @param vector2 second individual
 * @return recombination sites
 */
arma::mat twoVecHetDet(const arma::mat &vector1, const arma::mat & vector2)
{
	arma::mat sum2vec(1, vector1.size());
	sum2vec.zeros();
	sum2vec = vector1 + vector2;

	arma::mat result(1, vector1.size());
	result.zeros();

	int previous = 0, previousLocation = 0;
	for (uint i = 0; i < vector1.size(); i++)
	{
		if (sum2vec[i] == 0 || sum2vec[i] == 2 || sum2vec[i] == 4)
		{
			previous = sum2vec[i];
			previousLocation = i;
			break;
		}
	}
	for (uint i = previousLocation; i < vector1.size(); i++)
	{
		if (sum2vec[i] == 0 || sum2vec[i] == 2 || sum2vec[i] == 4)
		{
			if (abs(previous - sum2vec[i]) == 2)
			{
				for (uint j = previousLocation; j < i; j++)
					result[j] = 1;
				previous = sum2vec[i];
				previousLocation = i;
			}
			else if (sum2vec[i] == 0 || sum2vec[i] == 2 || sum2vec[i] == 4)
			{
				previous = sum2vec[i];
				previousLocation = i;
			}

		}

	}

	return (result);
}

/**
 *
 * @param genotype
 * @param hetsite heterozygote loci in the sire
 * @return
 */
SEXP co(SEXP geno, SEXP hetsite)
{
	NumericMatrix genotype(geno);
	int n = genotype.nrow(), k = genotype.ncol();
	arma::mat Genotype(as<arma::mat>(genotype).begin(), n, k, TRUE);

	NumericVector heterozygoteSite(hetsite);
	arma::vec HeterozygoteSite(as<arma::vec>(heterozygoteSite).begin(), heterozygoteSite.length(), TRUE);

	arma::mat crossOver(as<arma::mat>(genotype).begin(), n, k, TRUE);
	crossOver.zeros();

	arma::mat tempCrossOver(as<arma::mat>(genotype).begin(), 1, k, TRUE);
	tempCrossOver.zeros();

	for (int i = 0; i < genotype.nrow() - 1; i++)
	{
		for (int j = i + 1; j < genotype.nrow(); j++)
		{
			tempCrossOver = twoVecHetDet(Genotype.row(i), Genotype.row(j));
			crossOver.row(i) = crossOver.row(i) + tempCrossOver;
			crossOver.row(j) = crossOver.row(j) + tempCrossOver;
		}
	}
	return (wrap(crossOver));
}

