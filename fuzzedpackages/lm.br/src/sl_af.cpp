//


#include "lmbr.h"




double Clmbr::sl_af( int mode )  const
// estimate significance level for broken line regression by Approximate F method
{
	int q;  
	if (th0ex) q = 2; else q = 1;  
	if (mode==2) q++;

	const double  errors_sq0 = qysq - prime_z*prime_z;

	double sL;
	if (variance_unknown) {
		const double  Fstat = (m-2)*1./q*fabs(errors_sq0/omega - 1.);
		sL =  1. -  Rf_pf( Fstat, q, m-2, 1, 0 ) ;
	}  else  {
		const double  CHIstat = fabs(errors_sq0 - omega);
		sL =  1. -  Rf_pchisq( CHIstat, q ,1,0) ;
	}

	return  sL;
}





double Clmbr::sl_af2( void )  const
{
	int q;
	if (th0ex) q = 3; else q = 2;

	double sL;
	if (variance_unknown) {
		const double  Fstat = (m-2)*1./q*fabs(lambdasq/omega - 1.);
		sL =  1. -  Rf_pf( Fstat, q, m-2, 1, 0 ) ;
	}  else  {
		const double  CHIstat = fabs(lambdasq - omega);
		sL =  1. -  Rf_pchisq( CHIstat, q ,1,0) ;
	}

	return  sL;
}

