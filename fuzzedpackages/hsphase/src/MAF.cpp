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
 * MAF.cpp
 *
 *  Created on: 16/01/2014
 *      Author: mhf
 */

#include "MAF.h"
SEXP MAFC(SEXP snp)
{
	NumericVector SNPs(snp);

	double z = 0, o = 0, t = 0, result = 0;
	for (int i = 0; i < SNPs.length(); i++)
	{
		if (SNPs[i] == 0)
		{
			z = z + 1;
		}
		if (SNPs[i] == 1)
		{
			o = o + 1;
		}
		if (SNPs[i] == 2)
		{
			t = t + 1;
		}
	}
	result = (z * 2 + o) / ((z + o + t) * 2);
	if (result > 0.5)
		result = 1 - result;
	return wrap(result);
}

