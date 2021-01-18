//


#include "lmbr.h"




void  Clmbr::set_alpha0( double a_0,  METHOD met )
// precalculate numbers and vectors based on alpha0, used in routines
{
	if ( !R_FINITE(a_0) )  stop( _("invalid 'alpha0' value") );

	if ( a_0 == alpha0  &&  th0 == th0a0 )  return;

	alpha0 = a_0;
	th0a0 = th0;

	Vector<double>  star_y(n);
	star_y = *psy - alpha0*(*psig1);


	if (Model==M1) {

		if (th0ex) {

			Vector<double>  pf0(n);
			pf0 = *psigx - th0*(*psig1);
			const double  yf0 = star_y*pf0;

			lambdasq = star_y*star_y - yf0*yf0/(pf0*pf0);
			if(lambdasq < 0.)  lambdasq= 0.; 
			lambda = sqrt( lambdasq );

		}  else  {

			Vector<double>  gsm0(n),  gfr0(n),  gbar0(n);
			gsm0 = gsm(th0,k0);  
			gfr0 = gfr(th0,k0);  
			gbar0 = gbar(th0,k0);  

			const double ysm0 = star_y*gsm0;
			const double yfr0 = star_y*gfr0;
			lambdasq = star_y*star_y - ysm0*ysm0 - yfr0*yfr0;
			if(lambdasq < 0.)  lambdasq= 0.; 
			lambda = sqrt( lambdasq );

			if ( ! (met==AF || met==AF2) )  {
				const double ybar0 = star_y*gbar0;
				const double c3 = (*pv1h*gfr0)*(*pxh*gsm0) - (*pv1h*gsm0)*(*pxh*gfr0);
				c1 = -lambda*c3;
				c2 = prime_z + ybar0*c3;
			}
		}
	}



	if (Model==M2) {

		if (th0ex) {

			lambdasq = star_y*star_y;
			lambda = sqrt( lambdasq );

		}  else  {

			Vector<double>  gfr0(n),  gbar0(n);
			gfr0 = gfr(th0,k0);  
			gbar0 = gbar_prime(th0,k0);  

			const double yfr0 = star_y*gfr0;
			lambdasq = star_y*star_y - yfr0*yfr0;
			if(lambdasq < 0.)  lambdasq= 0.; 
			lambda = sqrt( lambdasq );

			if ( ! (met==AF || met==AF2) ) {
				const double ybar0 = star_y*gbar0;
				const double c3 = *pv1h*gfr0;
				c1 = -lambda*c3;
				c2 = prime_z + ybar0*c3;
			}
		}
	}



	if(omega==0) c= 1;  else  c= sqrt(  max( 1 - omega/lambdasq , 0. )  );


	return;
}

