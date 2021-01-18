//

#include "lmbr.h"




double Clmbr::sl( double th_0,  METHOD met,  bool verbose )
{
	double  sL,  err = 0.;

	if(trivial)  { 
		const double  thmle= mle(false);
		if(  ISNAN(thmle)  ||  th_0==thmle  || 
				(Model==M2 && thmle==xs[0] && th_0 <= thmle)  )   sL= 1.;  else  sL= 0.;

	}  else  {

		set_theta0( th_0, met );

		double *const per= &err;

		if ( fabs(z) >= w )  sL= 1.;  else  {

			if (met==AF) sL = sl_af();
			if (met==AF2) sL = sl_af(2);
			if (met==GEO) sL = sl_geo(per);
			if (met==GEO2) {
				const double  ah = ahigh(GEO,th_0);
				set_alpha0( ah, met );
				sL = sl_geo2(per);
			}
			if (met==MC) sL = sl_mc();
		}

	}


	if (verbose)  {
		const int  reflect = copysign( 1, model_in );
		Rcout << "  SL= " << sL << _("  for theta0 = ") << th_0*reflect ;
		if( !trivial )  {
			Rcout <<  _("  by method ");
			if (met==AF)  Rcout << "AF" ;  
			if (met==GEO)  Rcout << "CLR" ;  
			if (met==GEO && !th0ex)  Rcout << "  int.er.< " << err ; 
			if (met==MC)  Rcout << "CLR-MC" ;
		}
		Rcout << endl;
	}


	return sL ;
}




double Clmbr::sl( double th_0,  double a0,  METHOD met,  bool verbose )
{
	double  sL,  err = 0.;

	if(trivial) { 
		const double  thmle= mle(false);
		if( ISNAN(thmle)  ||  ( thmle==xs[0] && thmle >= th_0 )  )  {
			const double  
				slope =  ( (*py)[ is[1] ] - (*py)[ is[0] ] )/( xs[1]-xs[0]),  
				intercept =  (*py)[ is[0] ]  - slope*xs[0],  
				amle =  slope*th_0 + intercept;
			if( fabs(a0 - amle) < zero_eq )  sL= 1.;  else  sL= 0.;
		}  else  {
			if( lambdasq < zero_eq )  sL= 1.;  else  sL= 0.;
		}	

	}  else  {
		set_theta0(th_0, met);
		set_alpha0(a0, met);

		double *const  per= &err;

		if (met==AF) sL = sl_af2();
		if (met==GEO) sL = sl_geo2(per); 
		if (met==MC) sL = sl_mc2();
	}


	if (verbose)  {
		const int  reflect = copysign( 1, model_in );
		Rcout << "  SL= " << sL << _(" for (th0,a0)= ( ") << th_0*reflect << ", " << a0 << " )" ;
		if( !trivial )  {
			Rcout <<  _("  by method ");
			if (met==AF)  Rcout << "AF" ;  
			if (met==GEO)  Rcout << "CLR  int.er.< " << err ;  
			if (met==MC)  Rcout << "CLR-MC" ;
		}
		Rcout << endl;
	}


	return sL;
}


