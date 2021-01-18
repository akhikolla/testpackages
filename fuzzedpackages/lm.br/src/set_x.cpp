//

#include "lmbr.h"



void Clmbr::set_x( void )
// check the 'x' input,  set 'ns' and 'px'
{
	int i;

	Vector<double>  x(n);
	int xcol;
	if(Model==M3) xcol=0; else xcol=1;
	for (i=0;i<n;i++) x[i] = *(x_in+xcol*n+i);
	const double  min_xdiff= (x[n-1]-x[0])*0.001;
	double  xib,  xi= x[0] - 1 - min_xdiff;

	for(i=0;i<n;i++) {
		xib= xi;
		xi= x[i];
		if ( !R_FINITE( xi ) )  stop( _("invalid 'x' value") );
		if( xib > xi )  stop( _("'x' values must be non-decreasing") );

// some test examples with near identical x values caused errors
// but warning seems more trouble than it's worth, so below warning is off for now
/*
		const double xdiff= xi - xib;
		if ( 0 < xdiff  &&  xdiff < min_xdiff  ) {
			Rcout << _("consider a repeat predictor value instead of the values") << endl;
			Rcout << xib << ",  " << xi << endl;
			Rf_warning( _("predictor values might be too close for reliable computations") );
		}
*/
	}


// count number of seperate 'x' values
	ns= 0;
	for(i=1;i<n;i++)  if( x[i] != x[i-1] )  ns++;
	ns++;


	bool  lack= false;
	if(Model==M1)  if( ns < 4 )  lack= true;
	if(Model==M2)  if( ns < 3 )  lack= true;
	if(Model==M3)  if( ns < 2 )  lack= true;
	if(variance_unknown)  if( m < 3 )  lack= true;
	if(lack)  stop( _("number of seperate 'x' values below minimum for changepoint inference") );

	*px = x;
	if( model_in < 0 )  for (i=0;i<n;i++)  (*px)[i] = -x[n-1-i];

	return;
}

