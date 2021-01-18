//

#include "lmbr.h"



double Clmbr::get_C( int k )  const
{
	int i = 0;
	double d = 1.;

	if ( k % 2 ) {
		i = (k-1)/2;
		while (i>0) { d *= i / (i - 0.5); i--;}
		d /= pi;
	} else {
		i = k/2 - 1;
		while (i>0) { d *= (i + 0.5) / i; i--;}
		d /= 2;
	}

	return d;
}



double Clmbr::fk( int k, double arg )  const
// calculate the f() function defined in K,S&Z eq.(13)
{
	double gv = 0.;
	if ( fabs(arg) < 1. )  gv = C[k-(m-2)]*pow( 1.- arg*arg, k/2. - 1.);
	return gv;
}



double Clmbr::F( int k, double arg )  const
// calculate the F() function defined in K,S&Z eq.(14)
// use its identity with the t-distribution
{

	double Fx;
	if (arg < -1. + zero_eq)
		Fx = 0.;
	else
		if (arg > 1. - zero_eq)
			Fx = 1.;
		else
			Fx = Rf_pt( arg*sqrt(k/(1-arg*arg)), k, 1, 0) ;

	return Fx;
}



double Clmbr::sF( int k, double arg )  const
// calculate the sF() "script F" function defined as integral of Fk(s) from s= -Inf to s= x,
// where Fk() is the F(k,x) function defined in K,S&Z (1991) eq.(14)
{
	if( k < 0 || ISNAN(arg) )  stop( _("'sF': invalid input") );

	double sFx;
	if (arg <= -1. + zero_eq)
		sFx= 0.;
	else
		if (arg >= 1. - zero_eq)
			sFx= 1.;
		else {
			const double r= 1-arg*arg; 
			double sum= 0;
			if( (k%2)==0 ) {
				double ai= r/2;
				for(int i=1;i<=k/2;i++) {
					sum += ai;
					ai *= r*(2*i-1)/(2*i+2);
				}
				sFx= (1 + arg - sum)/2;
			} else {
				double bi= r/3;
				for(int i=1;i<=(k-1)/2;i++) {
					sum += bi;
					bi *= r*(2*i)/(2*i+3);
				}
				sFx= arg/2 + ( arg*asin(arg) + sqrt(r)*(1 - sum) )/pi;
			}
		}

	return sFx;
}


