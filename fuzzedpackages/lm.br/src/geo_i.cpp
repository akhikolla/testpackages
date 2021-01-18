
//   integrands in  Knowles, Siegmund, Zhang's  geometric-expectation formulae for conditional
//   likelihood-ratio significance levels 


#include "lmbr.h"





double Clmbr::amu_by_Omega( double th,  int k)  const
{
	if ( k==k0 || fabs(th-th0) < zero_eq )  return Inf;

	const double  ro = rho(th,k),  zwr = fabs(z-w*ro);
	const double  rosq = rhosq(th,k),  r = 1 - rosq;

	if( !R_FINITE(th) && !ISNAN(th) )  {
		if( B[k]-rosq < 0 )  return 0.;
		return  zwr * sqrt( (B[k]-rosq) / (1-B[k]) / r );

	}  else  {

		const double  drosq = drhosq(th,k),  OMsq = dgsq(th,k) - drosq/r;

		if( OMsq <= 0 )  return Inf;

		return  sqrt(drosq/OMsq)*zwr/r;
	}
}



double Clmbr::Emupr( double th,  int k)  const
// calculate (E-mu+)*pr, in Knowles,Siegmund,Zhang's geometric-expectation formula for variance unknown
// as a function of theta
{
	if ( k==k0 || fabs(th-th0) < zero_eq || (!R_FINITE(th) && !ISNAN(th)) )  return 0.;

	const double  rosq = rhosq(th,k),  ro = rho(th,k),  r=1-rosq,  zz = 1.-z*z,  wzr = w-z*ro;
	const double  drosq = drhosq(th,k),  zwr = fabs(z-w*ro);
	
	const double  pisq = zz - wzr*wzr/r;
	if ( pisq <= 0. )  return  0.;

	const double  OMsq = dgsq(th,k) - drosq/r;
	if( OMsq <= 0. )  return 0.;

	const double  tau =  sqrt( OMsq * pisq ),  ambt = sqrt(drosq) * zwr / r / tau ;

	if( ambt >= 1. )  return 0.;

	const double  rr= sqrt(r*zz), g= wzr/rr, pr= fk(m-2,g)/rr;

	return  pr * tau * sF(m-3,-ambt) ;
}



double  Clmbr::Emupr_vk( double th,  int k)  const
// calculate (E-mu+)*pr, the integrand in geometric-expectation formula for variance known
// as a function of theta
{
	if ( k==k0  ||  fabs(th-th0) < zero_eq  ||  (!R_FINITE(th) && !ISNAN(th)) )  return 0.;

	const double  rosq = rhosq(th,k),  r=1-rosq,  rr = sqrt(r),  ro= rho(th,k);
	const double  zwr= fabs(z-w*ro);

	const double  drosq = drhosq(th,k),  namu = -zwr*sqrt(drosq)/r ;

	const double  OMsq = dgsq(th,k) - drosq/r;
	if (OMsq <= 0.)  return 0.; 

	const double  OM = sqrt( OMsq ),  mbO = -zwr*sqrt(drosq/OMsq)/r;

	const double  pr = Rf_dnorm4( (w-z*ro)/rr ,0,1,0)/rr;

	const double  pn_mbO =  Rf_pnorm5(mbO ,0,1,1,0) ,   dn_mbO = Rf_dnorm4(mbO ,0,1,0) ;

	return  ( namu*pn_mbO + OM*dn_mbO )*pr;
}







