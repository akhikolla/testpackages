//


#include "lmbr.h"




void  Clmbr::set_SL( double cSL )
// set the critical significance level to determine confidence intervals or regions
// traditionally called the alpha-value in statistics literature
// precalculate numbers used in the 'a_af' subroutine 
{
	if ( ISNAN(cSL) || cSL<=0. || cSL>=1. )  stop( _("invalid 'SL' value") );

	if (cSL != prev_SL) {

		SL = cSL;
		prev_SL = SL;

		cFex =Rf_qf(1-SL,3,m-2,1,0); 
		cCHIex =Rf_qchisq(1-SL,3 ,1,0);
		cF =Rf_qf(1-SL,2,m-2,1,0); 
		cCHI =Rf_qchisq(1-SL,2 ,1,0);
	}

	x_vu_ex =(1+cFex*3./(m-2))*omega; 
	x_vk_ex =cCHIex + omega;
	x_vu =(1+cF*2./(m-2))*omega; 
	x_vk =cCHI + omega;


	return;
}

