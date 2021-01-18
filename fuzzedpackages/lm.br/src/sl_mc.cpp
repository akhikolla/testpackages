//

#include "lmbr.h"




double Clmbr::sl_mc(void)  const
// calculate significance level by CLR, Monte Carlo evaluation method
{
	const double  acc= tol_sl_abs;
	double  th0print = th0;   if(model_in  < 0 ) th0print = -th0;
	Rcpp::Function Rflush("flush.console");

	Rcout << endl << _("MC evaluation of conditional likelihood-ratio SL") << endl;
    Rcout << _("for ") << "theta0= " << th0print << ",  " << _("target accuracy =  ") << acc << ":" << endl << endl;
	Rcout << setw(10) << "iteration" << setw(12) << "est. SL" << setw(15) << "est. acc." << endl;
	Rflush();


	const double  wsq= w*w; 
	double  s0, testw;
	if (variance_unknown) { s0 = z/sqrt( 1.- z*z );  testw = wsq/(1-z*z); } 
					else  { s0 = z; testw = wsq; }

	bool shortcut= true;
	const int  N = 10000000;
	int  it,  count= 0;
	const double  tstart = time( NULL );
	double  pstart= tstart;

	GetRNGstate();

	for (it=1; it<N+1; it++)
	{
		Vector<double>  s(m,0.);
		if (th0ex)  s[0] =norm_rand();  else  s[0] =0.;
		for (int i=1;i<m;i++) s[i] =norm_rand();
		if (variance_unknown) s = 1./sqrt(s*s) * s;
		if (!th0ex) s[0] =s0;

		const  Vector<double>&  rs = s;
		if ( m_ge_w(testw, rs) ) count++;

		double  pfinish=  time( NULL ),  ptime= pfinish - pstart;
		if ( it==N/10 || !(it%(N/5)) || ptime > 10 || (shortcut && ptime>1) )  {
			const double  p = 1.*count/it,  err_est = 2*sqrt( p*(1.-p)/it );
			if( !shortcut || (shortcut && err_est < acc) )
				Rcout << setw(10) << it << setw(12) << p << setw(15) << err_est << endl; 
			Rflush();
			if( err_est < acc )  {it++; break;}
			pstart= pfinish;
			shortcut= false;
		}
	}

	PutRNGstate();

	Rcout << endl;
	it--;
	const double  sL = count*1./it;

	return sL; 
}





double Clmbr::sl_mc2(void)  const
// calculate significance level by clr, Monte Carlo evaluation method
// for changepoint = th0 and   alpha = alpha0 
{
	const double  acc= tol_sl_abs;
	double th0print = th0;   if(model_in < 0) th0print = -th0;
	Rcpp::Function Rflush("flush.console");
 
    Rcout << endl << _( "MC evaluation of conditional likelihood-ratio SL") << endl;
	Rcout << _("for ") << "(th0,a0)= (" << th0print << "," << alpha0 << "),  " 
		<< _("target accuracy =  ") << acc << ":" << endl << endl;
	Rcout << setw(10) << "iteration" << setw(14) << "est. SL" << setw(15) << "est. acc." << endl;
	Rflush();


	double  Fc=0, sL=0; 
	if (variance_unknown) {if (th0ex) Fc =F(m,-c); else  Fc =F(m-1,-c);} 
			else  Fc = Rf_pnorm5(-lambda*c ,0,1,1,0) ;


// generate mock results

	bool shortcut= true;
	const int  N = 10000000;
	int  it,  count= 0;
	const double  tstart = time( NULL );
	double  pstart= tstart,  sum=0.,  sumsqs=0.;

	GetRNGstate();

	for (it=1; it<N+1; it++)
	{
		const double xi = 2*c*(unif_rand()-0.5); 

		double z_tilde;
		if (th0ex)  z_tilde= 0.;  else  z_tilde = xi*c1 + c2; 

		const double  deltasq = lambdasq*(1-xi*xi) + z_tilde*z_tilde;

		double  z_;
		if (th0ex)  z_ = 0.;  else  { if (variance_unknown) z_ = z_tilde/sqrt(deltasq);  else  z_ = z_tilde; }

		double  wsq;
		if (variance_unknown)  wsq = 1 - omega/deltasq;  else  wsq = deltasq - omega;
		if (wsq<0.)  wsq= 0.;


		double  s0, testw;
		if (variance_unknown) { s0 = z_/sqrt(1-z_*z_); testw = wsq/(1-z_*z_); }
						else  { s0 = z_; testw = wsq; }


		Vector<double> s(m);
		if (th0ex)  s[0] = norm_rand();  else  s[0] = 0.;
		for (int i=1;i<m;i++) s[i] = norm_rand();
		if (variance_unknown) s =  1./sqrt(s*s) * s;
		if (!th0ex) s[0] = s0;

		const Vector<double>  &rs = s;
		if ( m_ge_w(testw, rs) )  {
			count++;
			double den;
			if (variance_unknown) {if (th0ex) den =fk(m,xi); else  den =fk(m-1,xi);} 
					else  den = Rf_dnorm4(lambda*xi, 0,1,0) ;
			sum += den;
			sumsqs += den*den;
		}


		double  pfinish=  time( NULL ),  ptime= pfinish - pstart;
		if ( it==N/10 || !(it%(N/5)) || ptime > 10 || (shortcut && ptime>1) )  {
			const double  p = 1.*sum/it,  sd = sqrt( (sumsqs/it - p*p)/it );
			double  err_est = 2*c*sd;
			if (!variance_unknown)  err_est *= lambda;
			if (variance_unknown) sL = 2*Fc + 2*c*sum/it; else sL = 2*Fc + lambda*2*c*sum/it;
			if( !shortcut || (shortcut && err_est < acc) )
				Rcout << setw(10) << it << setw(14) << sL << setw(15) << err_est << endl; 
			Rflush();
			if ( err_est < acc )  {it++; break;}
			pstart= pfinish;
			shortcut= false;
		}
	}

	PutRNGstate();

	it--;
	Rcout << endl;
	if (variance_unknown)  sL = 2*Fc + 2*c*sum/it;  else  sL = 2*Fc + lambda*2*c*sum/it;

	return sL;
}





bool Clmbr::m_ge_w( double wsq, const Vector<double> &s)  const
// Check whether  max(<gam(theta).s>^2)  >=  wsq  for some 'theta'.  Use pre-calculated  M*"stump of 1"  vectors 
// for the dot products of 'gamma' with 's',  by  gamma*u = gamma*(M-transpose*s) = (M*gamma)*s .  
{

	if( th0ex )  {		
// if th0ex= TRUE, vectors not pre-multiplied by 'M' matrix

		double  sf= (xs[ns-1]-xs[ns-2])*(s*pq1[ns-1]),  mk= sf*sf/qff[ns-1];
		if ( mk >= wsq )  return true; 

		for ( int k= ns-2; k > k1; k-- ) 
		{
			const double  s1 = s*pq1[k],  sx = sf + s1*xs[k];  
			if(k>0) sf = sx - s1*xs[k-1];

			const double  sa = sx*qx1[k] - s1*qxx[k],  sb = sx*q11[k] - s1*qx1[k],  thk = sa/sb;

			if( k > 0 )  {
				if( xs[k-1] < thk  &&  thk < xs[k] )  mk = (sx*sb - s1*sa)/ck[k];  else  mk = sf*sf/qff[k];
			}  else  {
				if( thk < xs[k] )  mk = (sx*sb - s1*sa)/ck[k];  else  mk= (s*(*pv1h)) * (s*(*pv1h)) ;		// lim sup
			}

			if ( mk >= wsq )  return true;
		}

	}  else  {

		double  sf= (xs[ns-1]-xs[ns-2])*(s*pmq1[ns-1]),  mk= sf*sf/qff[ns-1];
		if ( mk >= wsq )  return true; 

		for ( int k= ns-2; k > k1; k-- ) 
		{
			const double  s1 = s*pmq1[k],  sx = sf + s1*xs[k];  
			if(k>0) sf = sx - s1*xs[k-1];

			const double  sa = sx*qx1[k] - s1*qxx[k],  sb = sx*q11[k] - s1*qx1[k],  thk = sa/sb;

			if( k > 0 )  {
				if( xs[k-1] < thk  &&  thk < xs[k] )  mk = (sx*sb - s1*sa)/ck[k];  else  mk = sf*sf/qff[k];
			}  else  {
				if( thk < xs[k] )  mk = (sx*sb - s1*sa)/ck[k];  else  mk= (s*(*pm1h)) * (s*(*pm1h)) ;		// lim sup
			}

			if ( mk >= wsq )  return true;
		}

	}

	return false;
}



