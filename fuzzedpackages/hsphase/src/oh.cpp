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


#include "oh.h"

#define mat imat

using namespace Rcpp;

int make_na(arma::mat &genotype)
{
	for (unsigned int i = 0; i < genotype.n_rows; i++)
	{
		for (unsigned int j = 0; j < genotype.n_cols; j++)
		{
			if (genotype(i, j) == 1)
			{
				genotype(i, j) = 9;
			}
		}

	}
	return 0;
}

arma::mat twoFreq(arma::mat data)
{
	// Rprintf("%d\n", data.n_rows);
	/**
	 * data.n_row: NO. sire
	 * data.n_col: NO. animals
	 */

	arma::mat result = arma::zeros<arma::mat>(1, data.n_rows);
	int frequency = 0;
	for (unsigned int i = 0; i < data.n_rows; i++)
	{
		frequency = 0;
		for (unsigned int j = 0; j < data.n_cols; j++)
		{
			if (data(i, j) == 2)
			{
				frequency = frequency + 1;
			}
		}
		result(0, i) = frequency;
	}
	return result;
}

arma::mat vecMinusMat(arma::mat genotype, arma::mat sire, int index)
{
	arma::mat result = sire;
	arma::mat a = genotype.row(index);

	for (unsigned int i = 0; i < sire.n_rows; i++)
	{
		result.row(i) = abs(a - sire.row(i));
	}
	return result;
}

SEXP oh(SEXP geno, SEXP sire)
{
	NumericMatrix genotype(geno);
	int n = genotype.nrow(), k = genotype.ncol();
	arma::mat Genotype(as<arma::mat>(genotype).begin(), n, k, TRUE);

	NumericMatrix sireGeno(sire);
	n = sireGeno.nrow(), k = sireGeno.ncol();
	arma::mat SireGeno(as<arma::mat>(sireGeno).begin(), n, k, TRUE);

	/**
	 *  saveing the final results
	 */
	arma::mat result = arma::zeros<arma::mat>(genotype.nrow(), sireGeno.nrow());

	arma::mat res2 = arma::zeros<arma::mat>(sireGeno.nrow(), genotype.nrow());
	make_na(Genotype);

    #pragma omp parallel for private(res2) schedule(dynamic)
	for (int i = 0; i < genotype.nrow(); i++)
	{
		//  Rprintf("%d\n", i);
		res2 = vecMinusMat(Genotype, SireGeno, i);
		result.row(i) = twoFreq(res2);
	}
	List res;
	res["OH"] = result;
	return res;
}

