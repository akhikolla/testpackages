//

#include "lmbr.h"




const int  bis_it_limit= 50;	// maximum number of iterations in bisection routines



double Clmbr::bisect( double a, double b, double (Clmbr::*fn)(double,int), int k, double value, double crit)
// find  x  such that  value < fn(x) < value + crit   if  crit > 0 ,   or   value - crit < fn(x) < value   if   crit < 0
{
	double  x1= a, x2= b, f1 = (this->*fn)(x1,k) - value,  f2 = (this->*fn)(x2,k) - value;
	if ( f1*f2>0 || f1==f2 || ISNAN(f1*f2) )
		stop( _("'bisect' cannot find interim point from starting values") );
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration < bis_it_limit) ) {
		const double  xmean = (x1+x2)/2,  fx = (this->*fn)(xmean,k)-value;
		if(f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; } 
		iteration++;
	}
	if(iteration==bis_it_limit)  Rf_warning( _("'bisect' failed to reach tolerance after maximum number of iterations") );
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}



double Clmbr::bisect( double a, double b, double (Clmbr::*fn)(double,int) const, int k, double value, double crit)  const
// "const" version of above routine
{
	double  x1= a, x2= b, f1 = (this->*fn)(x1,k) - value,  f2 = (this->*fn)(x2,k) - value;
	if ( f1*f2>0 || f1==f2 || ISNAN(f1*f2) ) 
		stop( _("'bisect' const  cannot find interim point from starting values") );
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration < bis_it_limit) ) {
		const double  xmean = (x1+x2)/2,  fx = (this->*fn)(xmean,k)-value;
		if(f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; } 
		iteration++;
	}
	if(iteration==bis_it_limit)  
		Rf_warning( _("'bisect' const  failed to reach tolerance after maximum number of iterations") );
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}




double Clmbr::bisect_sl( double a, double b, METHOD met, double crit)
// use bisection to find x such that   SL < sl(x) < SL + crit   if  crit > 0 ,
// or    SL - crit < sl(x) < SL   if  crit < 0 ,   used in 'ci' routine
{
	double  x1= a, x2= b, f1= sl(x1,met,false)-SL, f2= sl(x2,met,false)-SL;
	if( fabs(f1)<zero_eq && fabs(f1-f2)<zero_eq) return (x1+x2)/2.;
	const double p12 = f1*f2;
	if ( p12>0 || f1==f2 || fabs(p12)>1 || ISNAN(p12) ) 
		stop( _("'bisect_sl' cannot find interim point from starting values") );
	int iteration=0;
	while ( fabs(x1-x2) > fabs(crit)  && (iteration < bis_it_limit) ) {
		const double  xmean= (x1+x2)/2,  fx= sl(xmean,met,false) -SL;
		if (f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; }
		iteration++;
	}
	if(iteration==bis_it_limit)  
		Rf_warning( _("'bisect_sl' failed to reach tolerance after maximum number of iterations") );
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}



