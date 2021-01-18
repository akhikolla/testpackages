//


#include "lmbr.h"




void  Clmbr::set_tol( double tol )
// set accuracy and scale parameters
{
	if ( ISNAN(tol) || tol<=0 || tol>=1 )  stop( _("invalid 'tol' value") );
	

	subints = 5;	// average number of subintervals per data interval for grid searches

	tol_rho = 0.0001;		//  tolerance for finding 'rho' values by 'bisect' or 'rho_inv' 

	tol_sl_abs = tol;					// maximum absolute error in significance level estimates
	tol_sl_rel = min( 10*tol, 0.01 );	// maximum relative error in significance level estimates


	int i;

// maximum error in x- boundaries
	tol_xb = (xs[ns-1] - xs[0])*tol_sl_rel/64.;	
	i=1;  while( tol_xb < ldexp(1.,-i) )  i++;  tol_xb = ldexp(1.,-i);


// maximum error in y- boundaries
	double  maxY= -Inf,  minY= Inf;
	for(i=0;i<n;i++)  {  if( (*py)[i] > maxY )  maxY= (*py)[i];  if( (*py)[i] < minY )  minY= (*py)[i];  }
	const double dY= maxY - minY;
	tol_yb = dY*tol_sl_rel/64.;	
	i=1;  while( tol_yb < ldexp(1.,-i) )  i++;  tol_yb = ldexp(1.,-i);


// precision for integration limits
	inc_x = tol_xb;
	for( i= max( k1, 0 );  i < ns-2;  i++ )  { 
		double  inc= ( xs[i+1] - xs[i] )/subints;  
		if( inc < inc_x )  inc_x= inc; 
	}
	i=1;  while( inc_x < ldexp(1.,-i) )  i++;  inc_x = ldexp(1.,-i);


// increment for x- grid searches, and for printout of confidence regions by 'cr'
	double inc = ( xs[ns-1] - xs[0] )/(ns-1)/subints;
	double inc_seed[3] = { 5., 2., 1. };
	double inc10 = 1.;
	while ( inc > inc10) inc10 *= 10;
	i = 0;  while ( inc < inc_seed[i]*inc10 - zero_eq ) { i++; if(i==3) { i=0; inc10 /= 10; }  }
	inc = inc_seed[i]*inc10;
	xinc = inc;


// starting increment for y- grid searches
	inc_y = dY/128;
	i=1;  while( inc_y < ldexp(1.,-i) )  i++;  inc_y = ldexp(1.,-i);


// output digits and minimum increment
	const int  digits= 6;
	Rcout << setprecision( digits );
	rel_print_eps =  pow( 10., -(digits-1) );


// check for trivial case
	trivial= false;
	if ( variance_unknown  &&  omega/m < zero_eq )  trivial= true;


	return;
}


