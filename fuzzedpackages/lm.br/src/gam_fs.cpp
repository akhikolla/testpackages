//

#include "lmbr.h"


 
Vector<double> Clmbr::gam( double th, int k )  const
// calculate gamma
{
	if ( (Model==M2 && (!R_FINITE(th) && !ISNAN(th)) ) || (Model==M1 && th <=xs[0]) || xs[ns-1] <= th) return *nan_m; else
		if (xs[ns-2] <= th  &&  th < xs[ns-1]) return *puqen; else
			if (Model==M1  &&  (xs[0] < th  &&  th <= xs[1]) ) return *puqe1; else
				if (Model==M2  &&  th <= xs[0])  return *puqx;  else
					if ( Model==M3  &&  (!R_FINITE(th) && !ISNAN(th)) )  return  *pv1h;  else

	return  1./sqrt( ff(th,k) ) * q_f(th,k);
}


 
Vector<double> Clmbr::gfr( double th, int k )  const
// gamma frown
{
	if (th >= xs[ns-1]) return *nan_m1; else
		if (xs[ns-2] <= th  &&  th < xs[ns-1]) return *pusen; else

	return  1./sqrt(sf(th,k)*sf(th,k)) * sf(th,k);
}


 
Vector<double> Clmbr::gsm( double th, int k )  const
// this function "gamma-smile" used in Model M1 only
{
	if (th <= xs[0])  return *nan_m1;  else
		if (cov_matrix_diagonal) {
			if (xs[0] < th  &&  th < xs[1])  return *pnuse1;  else
				return  1./sqrt(sfc(th,k)*sfc(th,k)) * sfc(th,k);
		} else {
			if (xs[0] < th  &&  th <= xs[1]) {
				const double e1fr = *pnse1*gfr(th,k);
				return  1./sqrt(se1sq - e1fr*e1fr) * (*pnse1 - e1fr*gfr(th,k));
			} else {
				Vector<double> f1(m1);
				f1 = sfc(th,k) - (sfc(th,k)*sf(th,k))/(sf(th,k)*sf(th,k)) * sf(th,k);
				return  1./sqrt(f1*f1) * f1;
			}
		}
}


 
Vector<double> Clmbr::gbar( double th, int k )  const
// gamma bar
{
	if ( (Model==M1 && th <= xs[0])  ||  xs[ns-1] <= th )  return *nan_m1;  else  {

		Vector<double> fbar(m1);
		fbar = *pv1h - (*pv1h*gsm(th,k))*gsm(th,k) - (*pv1h*gfr(th,k))*gfr(th,k);

		return  1./sqrt(fbar*fbar) * fbar;
	}
}


 
Vector<double> Clmbr::gbar_prime( double th, int k )  const
// gamma bar prime
{
	if ( (Model==M1 && th <= xs[0])  ||  xs[ns-1] <= th )  return *nan_m1;  else  {

		Vector<double> fbar(m1);
		fbar = *pv1h - (*pv1h*gfr(th,k))*gfr(th,k);

		return  1./sqrt(fbar*fbar) * fbar;
	}
}


 
Vector<double> Clmbr::q_f( double th, int k )  const
{
	return  pqx[k] - th*pq1[k];
}


 
Vector<double> Clmbr::sf( double th, int k )  const
{
	return  psx[k] - th*ps1[k];
}



Vector<double> Clmbr::sfc( double th, int k )  const
{
	return  (*psigx - psx[k]) - th*(*psig1 - ps1[k]);
}



double Clmbr::ff( double th, int k )  const
{
	return  qxx[k] + (q11[k]*th - 2*qx1[k])*th;
}


