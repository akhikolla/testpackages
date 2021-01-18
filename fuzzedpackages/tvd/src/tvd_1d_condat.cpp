// tvd_1d_condat.cpp -- Implementation of the Condat algorithm for 1D TVD
// 
// Copyright 2014 Mark Pinese
//
// This file is distributed under the terms of the Eclipse Public 
// License v1.0, available at:
// https://www.eclipse.org/org/documents/epl-v10.html

#include <Rcpp.h>


// Rcpp still uses ints as indices and size objects internally.
// Hopefully this will be upgraded soon.
typedef int Rcpp_index;


// Rcpp implementation of the 1D TVD algorithm in:
// Condat, L (2012) A Direct Algorithm for 1D Total Variation Denoising
// http://hal.inria.fr/docs/00/67/50/43/PDF/condat_killer_tv.pdf
// Accessed 11 August 2014.
// [[Rcpp::export]]
Rcpp::NumericVector tvd_1d_condat_worker(Rcpp::NumericVector& y, double lambda)
{
	// NOTE: in wrapping R code, verify that y.size() <= 2^32-1.  Although
	// newer versions of R get around this index limit, Rcpp is still 
	// limited to int indices.
	Rcpp_index N = y.size();
	Rcpp_index k, k0, kn, kp, i;
	double vmin, vmax, umin, umax;
	Rcpp::NumericVector x(N);

	k = k0 = kn = kp = 1;
	vmin = y[0] - lambda;
	vmax = y[0] + lambda;
	umin = lambda;
	umax = -lambda;

	while (true)
	{
		// b
		if (k == N)
		{
			x[k-1] = vmin + umin;
			break;
		}

		if (y[k] + umin < vmin - lambda)
		{
			// b1
			for (i = k0-1; i < kn; i++)
				x[i] = vmin;
			kn++;
			k = k0 = kp = kn;
			vmin = y[k-1];
			vmax = y[k-1] + 2*lambda;
			umin = lambda;
			umax = -lambda;
		}
		else if (y[k] + umax > vmax + lambda)
		{
			// b2
			for (i = k0-1; i < kp; i++)
				x[i] = vmax;
			kp++;
			k = k0 = kn = kp;
			vmin = y[k-1] - 2*lambda;
			vmax = y[k-1];
			umin = lambda;
			umax = -lambda;
		}
		else
		{
			// b3
			k++;
			umin = umin + y[k-1] - vmin;
			umax = umax + y[k-1] - vmax;
			if (umin >= lambda)
			{
				// b31
				vmin = vmin + (umin - lambda)/(static_cast<double>(k - k0 + 1));
				umin = lambda;
				kn = k;
			}
			if (umax <= -lambda)
			{
				// b32
				vmax = vmax + (umax + lambda)/(static_cast<double>(k - k0 + 1));
				umax = -lambda;
				kp = k;
			}
		}

		// c
		if (k == N)
		{
			if (umin < 0)
			{
				// c1
				for (i = k0-1; i < kn; i++)
					x[i] = vmin;
				kn++;
				k = k0 = kn;
				vmin = y[k-1];
				umin = lambda;
				umax = y[k-1] + lambda - vmax;
			}
			else if (umax > 0)
			{
				// c2
				for (i = k0-1; i < kn; i++)
					x[i] = vmax;
				kp++;
				k = k0 = kp;
				vmax = y[k-1];
				umax = -lambda;
				umin = y[k-1] - lambda - vmin;
			}
			else
			{
				// c3
				for (i = k0-1; i < N; i++)
					x[i] = vmin + umin/(static_cast<double>(k - k0 + 1));
				break;
			}
		}
	}

	return x;
}
