#include <math.h>
#include "rcor_bm.h"
#include "tnorms.h"

using namespace Rcpp;

RcppExport SEXP rcor_matrix_classical(SEXP vx, SEXP r)
{
    int i, j, n;
	
    NumericVector x(vx);
	
    n = x.length();
	
    NumericMatrix Rx(n, n);
		
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(i) < x(j))
	    {
		Rx(i, j) = 1;
		Rx(j, i) = 0;
	    }
	    else if (x(i) > x(j))
	    {
		Rx(j, i) = 1;
		Rx(i, j) = 0;
	    }
	}
    }
	
    return Rx;
}

RcppExport SEXP rcor_matrix_linear(SEXP vx, SEXP r)
{
    int i, j, n;
    double rx;
    
    NumericVector x(vx);
    
    n = x.length();

    NumericMatrix Rx(n, n);
	
    rx = NumericVector(r)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(i) < x(j))
	    {
		Rx(i, j) = MIN(1, (x(j) - x(i)) / rx);
		Rx(j, i) = 0;
	    }
	    else
	    {
		Rx(j, i) = MIN(1, (x(i) - x(j)) / rx);
		Rx(i, j) = 0;
	    }
	}
    }
	
    return Rx;
}

RcppExport SEXP rcor_matrix_exp(SEXP vx, SEXP r)
{
    int i, j, n;
    double rx;
	
    NumericVector x(vx);
	
    n = x.length();

    NumericMatrix Rx(n, n);
	
    rx = NumericVector(r)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(i) < x(j))
	    {
		Rx(i, j) = 1 - exp((x(i) - x(j)) / rx);
		Rx(j, i) = 0;
	    }
	    else 
	    {
		Rx(j, i) = 1 - exp((x(j) - x(i)) / rx);
		Rx(i, j) = 0;
	    }
	}
    }
	
    return Rx;
}

RcppExport SEXP rcor_matrix_gauss(SEXP vx, SEXP r)
{
    int i, j, n;
    double rx, x_res, x_diff;

    NumericVector x(vx);
	
    n = x.length();

    NumericMatrix Rx(n, n);
	
    rx = NumericVector(r)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    x_diff = x(i) - x(j);		
	    x_res = 1 - exp(-0.5 * (x_diff / rx) * (x_diff / rx));
		
	    if (x(i) < x(j))
	    {
		Rx(i, j) = x_res;
		Rx(j, i) = 0;
	    }
	    else
	    {
		Rx(j, i) = x_res;
		Rx(i, j) = 0;
	    }
	}
    }
	
    return Rx;
}

RcppExport SEXP rcor_matrix_epstol(SEXP vx, SEXP r)
{
    int i, j, n;
    double rx;

    NumericVector x(vx);
	
    n = x.length();

    NumericMatrix Rx(n, n);
	
    rx = NumericVector(r)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(j) > (x(i) + rx))
	    {
		Rx(i, j) = 1;
		Rx(j, i) = 0;
	    }
	    else if (x(i) > (x(j) + rx))
	    {
		Rx(j, i) = 1;
		Rx(i, j) = 0;
	    }
	    else
	    {
		Rx(i, j) = Rx(j, i) = 0;
	    }
	}
    }
	
    return Rx;
}

RcppExport SEXP rcor_matrices_classical(SEXP vx, SEXP vy, SEXP r1, SEXP r2)
{
    int i, j, n;
	
    NumericVector x(vx);
    NumericVector y(vy);
	
    n = x.length();
	
    NumericMatrix Rx(n, n);
    NumericMatrix Ry(n, n);	
		
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(i) < x(j))
	    {
		Rx(i, j) = 1;
		Rx(j, i) = 0;
	    }
	    else if (x(i) > x(j))
	    {
		Rx(j, i) = 1;
		Rx(i, j) = 0;
	    }
				
	    if (y(i) < y(j))
	    {
		Ry(i, j) = 1;
		Ry(j, i) = 0;
	    }
	    else if (y(i) > y(j))
	    {
		Ry(j, i) = 1;
		Ry(i, j) = 0;
	    }
	}
    }
	
    List ret = List::create(Named("Rx") = Rx, Named("Ry") = Ry);
	
    return ret;
}

RcppExport SEXP rcor_matrices_linear(SEXP vx, SEXP vy, SEXP r1, SEXP r2)
{
    int i, j, n;
    double rx, ry;
	
    NumericVector x(vx);
    NumericVector y(vy);
	
    n = x.length();

    NumericMatrix Rx(n, n);
    NumericMatrix Ry(n, n);	
	
    rx = NumericVector(r1)[0];
    ry = NumericVector(r2)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(i) < x(j))
	    {
		Rx(i, j) = MIN(1, (x(j) - x(i)) / rx);
		Rx(j, i) = 0;
	    }
	    else
	    {
		Rx(j, i) = MIN(1, (x(i) - x(j)) / rx);
		Rx(i, j) = 0;
	    }
				
	    if (y(i) < y(j))
	    {
		Ry(i, j) = MIN(1, (y(j) - y(i)) / ry);
		Ry(j, i) = 0;
	    }
	    else 
	    {
		Ry(j, i) = MIN(1, (y(i) - y(j)) / ry);
		Ry(i, j) = 0;
	    }
	}
    }
	
    List ret = List::create(Named("Rx") = Rx, Named("Ry") = Ry);
	
    return ret;
}

RcppExport SEXP rcor_matrices_exp(SEXP vx, SEXP vy, SEXP r1, SEXP r2)
{
    int i, j, n;
    double rx, ry;
	
    NumericVector x(vx);
    NumericVector y(vy);
	
    n = x.length();

    NumericMatrix Rx(n, n);
    NumericMatrix Ry(n, n);	
	
    rx = NumericVector(r1)[0];
    ry = NumericVector(r2)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(i) < x(j))
	    {
		Rx(i, j) = 1 - exp((x(i) - x(j)) / rx);
		Rx(j, i) = 0;
	    }
	    else 
	    {
		Rx(j, i) = 1 - exp((x(j) - x(i)) / rx);
		Rx(i, j) = 0;
	    }
				
	    if (y(i) < y(j))
	    {
		Ry(i, j) = 1 - exp((y(i) - y(j)) / ry);
		Ry(j, i) = 0;
	    }
	    else 
	    {
		Ry(j, i) = 1 - exp((y(j) - y(i)) / ry);
		Ry(i, j) = 0;
	    }
	}
    }
	
    List ret = List::create(Named("Rx") = Rx, Named("Ry") = Ry);
	
    return ret;
}

RcppExport SEXP rcor_matrices_gauss(SEXP vx, SEXP vy, SEXP r1, SEXP r2)
{
    int i, j, n;
    double rx, ry;
    double x_res, y_res;
    double x_diff, y_diff;

    NumericVector x(vx);
    NumericVector y(vy);
	
    n = x.length();

    NumericMatrix Rx(n, n);
    NumericMatrix Ry(n, n);	
	
    rx = NumericVector(r1)[0];
    ry = NumericVector(r2)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    x_diff = x(i) - x(j);
	    y_diff = y(i) - y(j);		
		
	    x_res = 1 - exp(-0.5 * (x_diff / rx) * (x_diff / rx));
	    y_res = 1 - exp(-0.5 * (y_diff / ry) * (y_diff / ry));
		
	    if (x(i) < x(j))
	    {
		Rx(i, j) = x_res;
		Rx(j, i) = 0;
	    }
	    else
	    {
		Rx(j, i) = x_res;
		Rx(i, j) = 0;
	    }
				
	    if (y(i) < y(j))
	    {
		Ry(i, j) = y_res;
		Ry(j, i) = 0;
	    }
	    else
	    {
		Ry(j, i) = y_res;
		Ry(i, j) = 0;
	    }
	}
    }
	
    List ret = List::create(Named("Rx") = Rx, Named("Ry") = Ry);
	
    return ret;
}

RcppExport SEXP rcor_matrices_epstol(SEXP vx, SEXP vy, SEXP r1, SEXP r2)
{
    int i, j, n;
    double rx, ry;

    NumericVector x(vx);
    NumericVector y(vy);
	
    n = x.length();

    NumericMatrix Rx(n, n);
    NumericMatrix Ry(n, n);	
	
    rx = NumericVector(r1)[0];
    ry = NumericVector(r2)[0];
			
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    if (x(j) > (x(i) + rx))
	    {
		Rx(i, j) = 1;
		Rx(j, i) = 0;
	    }
	    else if (x(i) > (x(j) + rx))
	    {
		Rx(j, i) = 1;
		Rx(i, j) = 0;
	    }
	    else
	    {
		Rx(i, j) = Rx(j, i) = 0;
	    }
				
	    if (y(j) > (y(i) + ry))
	    {
		Ry(i, j) = 1;
		Ry(j, i) = 0;
	    }
	    else if (y(i) > (y(j) + ry))
	    {
		Ry(j, i) = 1;
		Ry(i, j) = 0;
	    }
	    else
	    {
		Ry(i, j) = Ry(j, i) = 0;
	    }
	}
    }
	
    List ret = List::create(Named("Rx") = Rx, Named("Ry") = Ry);
	
    return ret;
}
