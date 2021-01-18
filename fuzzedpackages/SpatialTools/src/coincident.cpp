#include "coincident_cpp.h"

using namespace Rcpp;

SEXP coincident_cpp(SEXP coords1, SEXP coords2, SEXP eps){

		#include <math.h>

		NumericMatrix x(coords1), y(coords2);
		NumericVector tol(eps);
		NumericMatrix coin(x.nrow(), 2);

		int count = 0;

		for(int i = 1; i <= x.nrow(); i++)
		{
			for(int j = 1; j <= y.nrow(); j++)
			{
				if(fabs(x(i - 1, 0) - y(j - 1, 0)) < tol[0])
				{
					if(fabs(x(i - 1, 1) - y(j - 1, 1)) < tol[0])
					{
						coin(i - 1, 0) = i;
						coin(i - 1, 1) = j;
						count++;
					}
				}
			}
		}

		NumericMatrix return_coin(count, 2);

		for(int i = 0; i < x.nrow(); i++)
		{
			if(coin(i, 0) > 0)
			{
				return_coin(return_coin.nrow() - count, 0) = coin(i, 0);
				return_coin(return_coin.nrow() - count, 1) = coin(i, 1);
				count--;
			}
		}

	 	return return_coin;
}
