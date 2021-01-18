//


#include "lmbr.h"





double Clmbr::geo_ex(void) const
// calculate significance level for theta0 outside data range (xs[0], xs[ns-1])
// using Knowles and Siegmund's geometric formulae
{
	if (variance_unknown)  return geo_vu_ex();  else  return geo_vk_ex();
}



double Clmbr::geo_vu_ex(void) const
// case variance unknown
{
	return  2*F(m-1,-w) + pow(1.-w*w, m*0.5-1)*Lgamma/pi;
}



double Clmbr::geo_vk_ex(void) const
// case variance known
{
	return  2*( Rf_pnorm5(-w, 0,1,1,0) + Lgamma/sqrt(2*pi) * Rf_dnorm4(w, 0,1,0) );
}



