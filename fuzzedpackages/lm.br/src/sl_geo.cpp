//

#include "lmbr.h"

void  igeo2( double * x,  const int n,  void *const ex );



double Clmbr::sl_geo(double * err)
// calculate significance level for changepoint = theta0 by CLR 
// using Knowles, Siegmund, Zhang's geometric formula
// assume model-object's  th0, y, z, w  values already set
{
	if (err!=0)  *err=0.;

	double  sL = 0.;

	if (th0ex)  {

		sL = geo_ex();

	}  else  {

		if ( fabs(z) >= w )  return 1.; 

		double  error=0., er=0, *const per = &er;

		if( th0 < xs[ns-2] )  {
			sL += geo(xs[ns-2], per);  error += *per;
			z = -z;
			sL += geo(xs[ns-2], per);  error += *per;
		}

		double thmin;
		if( k1 < 0 )  thmin= -Inf;  else  thmin= xs[k1];

		if ( th0 > thmin )  {
			sL += geo(thmin, per);  error += *per;
			z = -z;
			sL += geo(thmin, per);  error += *per;
		}

		if (err!=0)  *err = error;
	}

	return  min(sL,1.);
}




double Clmbr::sl_geo2(double * err)
// calculate significance level for  changepoint = (theta0,alpha0)  by CLR 
// using Knowles, Siegmund, Zhang's geometric-expectation formula
// assume  th0  and  y  values already set
{
	double sL;
	if (variance_unknown)  { if (th0ex)  sL =2*F(m,-c);  else  sL =2*F(m-1,-c); } 
		else  sL = 2 * Rf_pnorm5(-lambda*c ,0,1,1,0);

//  'Rdqag' parameters
//  here it integrates another numerical integral, so adds its average error-estimate 

	int  neval =0,  ier =0,  limit =100,  lenw = 4*limit,  last =0;
	int*  iwork= Calloc( limit, int );

	double  lower = -c,  upper = c,  epsabs = tol_sl_abs/2,  epsrel = tol_sl_rel/2,  
				result =0,  abserr =0;
	double *  work= Calloc( lenw, double );

	if (!variance_unknown)  epsabs /= lambda;

	int  ne = 0;
	double  er =0;
	void *  ex[3] = {  this,  &er, &ne  };

	Rdqags( igeo2, ex, &lower, &upper, &epsabs, &epsrel, &result, &abserr, &neval, 
				&ier, &limit, &lenw, &last, iwork, work );

	Free( iwork );  Free( work );

	double  integral= result,  error= abserr + er/ne;

	if (!variance_unknown)   integral *=lambda ,  error *=lambda; 

	sL += integral;

	if( err != 0 )  *err = error;

	return  min( sL, 1. );
}




void  igeo2(double * x, const int n, void *const ex)
// integrand for Rdqags 
{
	Clmbr  **const  ppObj  = static_cast< Clmbr **const > ( ex );
	Clmbr  *const  pObj = ppObj[0];

	double  **const  pper  = static_cast< double **const > ( ex );
	double  *const  per = pper[1];

	int  **const  ppne  = static_cast< int **const > ( ex );
	int  *const  pne = ppne[2];

	for (int i =0; i < n; i++ )  x[i] = (*pObj).prden( x[i], per );
	*pne += n;

	return;
}




double Clmbr::prden( double xi, double * err)
// integrand for  sl( theta, alpha )  generic formula
{
	double den;
	if(variance_unknown)  { if (th0ex) den= fk(m,xi);  else  den= fk(m-1,xi); } 
		else  den = Rf_dnorm4(lambda*xi ,0,1,0) ;


	double  z_tilde;
	if( th0ex )  z_tilde = 0;  else  z_tilde = xi*c1 + c2;

	const double  deltasq = lambdasq*(1-xi*xi) + z_tilde*z_tilde;

	if(variance_unknown)  z= z_tilde/sqrt(deltasq);  else  z= z_tilde;

	double wsq;
	if(variance_unknown)  wsq = 1 - omega/deltasq;  else  wsq = deltasq - omega;
	w = sqrt(max(wsq,0.));

	double  er= 0,  *const per= &er;  

	const double  pr = sl_geo(per);

	if( err != 0 )  *err += (*per)*den;

	return  pr*den;
}




