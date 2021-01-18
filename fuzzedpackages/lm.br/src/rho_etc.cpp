//

#include "lmbr.h"



double Clmbr::rho( double th)  const
{
	int k=0;
	while (k<ns && xs[k]<th) k++;
	return rho(th,k);
}



double Clmbr::rho( double th,  int k )  const
// compute the rho function,  accurately 
{
	if( th>=xs[ns-1] || (Model==M1  &&  th<=xs[0])  )  return NaN;  else  {	//1

// check if th0 or th is on an end-interval
		bool  th1= false, th2= false, th01= false, th02= false;
		if( xs[ns-2] <= th && th < xs[ns-1] ) th2= true;  else  {
			if( (Model==M1 && xs[0]<th && th<=xs[1]) || (Model==M2 && th<=xs[0]) || (Model==M3 && (!R_FINITE(th) && !ISNAN(th))) ) th1= true;
		}
		if( xs[ns-2] <= th0 && th0 < xs[ns-1] ) th02= true;  else  {
			if(  ( Model==M1 && xs[0]<th0 && th0<=xs[1] )  ||  ( Model==M2 && th0<=xs[0] ) || (Model==M3 && (!R_FINITE(th0) && !ISNAN(th0)))  ) th01 = true;
		}

		if( th==th0 || (th1 && th01) || (th2 && th02) )  return 1.;  else  {	//2

			if( (th01 || th02)  &&  !(th1 || th2) )  {	//3

				const double fsq= ff(th,k);
				const double num= qx0[k] - q10[k]*th;
				return  num/sqrt(fsq);

			}  else  {	//3 

				if( th1 || th2 )  {	//4

					if(th1)  return  g0u1;  else  return  g0u2;

				}  else  { 
										
					const double num= f0x[k] - f01[k]*th;
					return  num/sqrt(ff(th0,k0)*ff(th,k));

				}	//4
			} //3
		} //2
	} //1

}



double Clmbr::rhosq( double th,  int k )  const
// compute  rho squared,  accurately 
{
	if( th>=xs[ns-1] || (Model==M1  &&  th<=xs[0])  )  return NaN;  else  {	//1

// check if th0 or th is on an end-interval
		bool  th1= false, th2= false, th01= false, th02= false;
		if( xs[ns-2] <= th && th < xs[ns-1] ) th2= true;  else  {
			if( (Model==M1 && xs[0]<th && th<=xs[1]) || (Model==M2 && th<=xs[0]) || (Model==M3 && (!R_FINITE(th) && !ISNAN(th))) ) th1= true;
		}
		if( xs[ns-2] <= th0 && th0 < xs[ns-1] ) th02= true;  else  {
			if(  ( Model==M1 && xs[0]<th0 && th0<=xs[1] )  ||  ( Model==M2 && th0<=xs[0] ) || (Model==M3 && (!R_FINITE(th0) && !ISNAN(th0)))  ) th01 = true;
		}

		if( th==th0 || (th1 && th01) || (th2 && th02) )  return 1.;  else  {	//2

			if( (th01 || th02)  &&  !(th1 || th2) )  {	//3

				const double num= qx0[k] - q10[k]*th;
				return  num*num/ff(th,k);

			}  else  {	//3 

				if( th1 || th2 )  {	//4

					if(th1)  return  g0u1*g0u1;  else  return  g0u2*g0u2;

				}  else  { 
		
					const double  num= f0x[k] - f01[k]*th;
					return  (num/ff(th0,k0)) * (num/ff(th,k));

				}	//4
			} //3
		} //2
	} //1

}



double Clmbr::drho( double th,  int k )  const
{
	if(th>=xs[ns-1])  return NaN;  else
		if(Model==M1  &&  th<=xs[0])  return NaN;  else  
			if( !R_FINITE(th) && !ISNAN(th) )  return 0.;  else  {
				const double fsq= ff(th,k);
				double dro = (a0[k] - b0[k]*th)/sqrt(fsq)/fsq; 
				if (th < th0)  dro = -dro;
				return dro;
			}
}



double Clmbr::drhosq( double th,  int k )  const
{
	if(th>=xs[ns-1])  return NaN;  else
		if(Model==M1  &&  th<=xs[0])  return NaN;  else  
			if( !R_FINITE(th) && !ISNAN(th) )  return 0.;  else  {
				const double  fsq = ff(th,k);
				const double ab= a0[k] - b0[k]*th; 
				return ab*ab/(fsq*fsq*fsq);
			}
}



double Clmbr::dgsq( double th,  int k )  const
// norm of derivative of gamma, squared
{
	if(th>=xs[ns-1])  return NaN;  else
		if(Model==M1  &&  th<=xs[0])  return NaN;  else  
			if( ck[k]==0. )  return  0.;  else  {
				const double  fsq = ff(th,k);
				return ck[k]/fsq/fsq;
			}
}



double Clmbr::rho_inv( double s,  int k,  int hi_lo )  const
//   Returns 'th' such that rho(th) = s  and 'th' is in (x[k-1], x[k]);
// 'rho_inv' gives the same result for s and for -s, so it checks that  sign( rho(th) ) = sign(s); 
// result 'th' is a quadratic root, so if both roots are in the data interval then 
// it returns the lower value if 'hi_lo' < 0 and the greater value if 'hi_lo' > 0 .
{
	if( k1 <= k  &&  k < ns )  {	//1

		if( k== ns-1 )  {  if( fabs( rho(xs[ns-2],k) - s ) < zero_eq )  return xs[ns-2];  }  else  {	//2

			if( k== k1 )  { 
				if(k>=0) if( fabs( rho(xs[k1],k) - s ) < zero_eq )  return  xs[k1];  
				if(k== -1) if( fabs( rho(-Inf,k) - s ) < zero_eq )  return  -Inf;

			}  else  {	//3

				if( s==1 ) {
					if(k>0) {
						if( xs[k-1] <= th0  &&  th0 <= xs[k] )  return th0;  
					} else {
						if( th0 <= xs[k] )  return th0;  
					} 
					if( k== ns-2  &&  k0==ns-1 )  return xs[ns-2];
					if( k==k1+1  &&  k0==k1 )  return xs[k1];  

				}  else  {	//4

					if( fabs( rho(xs[k],k) - s ) < zero_eq )  return  xs[k];
					const double  thZp = a0[k]/b0[k],  rZp = rho(thZp,k);
					if(k>0)  {
						if( fabs( rho(xs[k-1],k) - s ) < zero_eq )  return  xs[k-1];
						if( xs[k-1] < thZp  &&  thZp < xs[k] )  if( fabs( rZp - s ) < zero_eq )  return  thZp;
					}  else  {
						if( fabs( rho(-Inf,k) - s ) < zero_eq )  return  -Inf;
						if( thZp < xs[k] )  if( fabs( rZp - s ) < zero_eq )  return  thZp;
					}


					if( fabs(s) < zero_eq )  {	
						const double  th1 = f0x[k]/f01[k];   
						if( k > 0 )  {
							if (xs[k-1] <= th1  &&  th1 <= xs[k] )  return th1;
						}  else  {
							if ( th1 <= xs[k] )  return th1;
						}
					}


					double  a,  b,  c;

					if( k0==ns-1  ||  k0==k1 || (k1<0 && k0==0) || (k1>=0 && th0==xs[k1]) )  {
						const double  s2 = s*s;
						a = s2*q11[k] - q10[k]*q10[k];
						b = q10[k]*qx0[k] - s2*qx1[k]; 
						c = s2*qxx[k] - qx0[k]*qx0[k];
					}  else  {
						const double  s2f = s*s*ff(th0,k0);
						a = s2f*q11[k] - f01[k]*f01[k];
						b = f01[k]*f0x[k] - s2f*qx1[k];
						c = s2f*qxx[k] - f0x[k]*f0x[k];
					}

					const double  rad = b*b-a*c;

					if(  rad <= 0.  )   {    
						const double  th1 = -b/a;   
						if( k > 0 )  {
							if (xs[k-1] <= th1  &&  th1 <= xs[k] )  if( fabs( rho(th1,k) - s ) < tol_rho )  return th1;
						}  else  {
							if ( th1 <= xs[k] )  if( fabs( rho(th1,k) - s ) < tol_rho )  return th1;
						}

					}  else  {
						const double  rd = sqrt(rad);
						const double  th1 = (-b-rd)/a,  th2 = (-b+rd)/a;

						bool  i1 = false,  i2 = false;
						if( k > 0 ) {
							if (xs[k-1] <= th1  &&  th1 <= xs[k] )  if( fabs( rho(th1,k) - s ) < tol_rho )  i1 = true;
							if (xs[k-1] <= th2  &&  th2 <= xs[k] )  if( fabs( rho(th2,k) - s ) < tol_rho )  i2 = true;
						}  else  {
							if ( th1 <= xs[k] )  if( fabs( rho(th1,k) - s ) < tol_rho )  i1 = true;
							if ( th2 <= xs[k] )  if( fabs( rho(th2,k) - s ) < tol_rho )  i2 = true;
						}
						if (i1 && i2) { if (hi_lo < 0) return min(th1,th2);  else  return max(th1,th2); }
						if (i1) return th1;
						if (i2) return th2;
					}
				}	//4
			}	//3
		}	//2
	}	//1

// last ditch
	if( 0 <= k  && k < ns )  {
		const double  rb = rho( xs[k], k );
		double  ra;
		if( k > 0 )  ra= rho( xs[k-1], k );  else  ra= rho( -Inf, k );
		if( (ra-s)*(rb-s) < 0. )  {
			const double  thb = xs[k];
			double  tha;
			if( k > 0 )  tha = xs[k-1];  else  {
				tha = min( thb - 1., -1. );
				ra = rho(tha,k);
				while( (ra-s)*(rb-s) < 0. ) { tha *= 2; ra= rho(tha,k); }
			}
			return  bisect( tha, thb, &Clmbr::rho, k, s, tol_rho );
		}
	}


	stop ( _("'rho_inv' no inverse for given 'rho' in data interval") );
	return NaN;
}


