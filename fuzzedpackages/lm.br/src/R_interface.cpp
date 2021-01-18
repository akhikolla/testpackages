//

#include "lmbr.h"

#define CLmsg			_( "confidence level must be between 0 and 1" )
#define methods2msg		_( "'method' must be 1 or 2" )
#define methods3msg		_( "'method' must be 1, 2 or 3" )
#define model_msg		_( "not applicable for this model" )





void  Clmbr::sl3R(  int met,  double tol,
	double theta_0 ) 
{ 
	METHOD MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			if(met==3)  MET=MC;  else  {
				stop( methods3msg );
			}
		}
	}

	const double tmp1 = tol_sl_abs, tmp2 = tol_sl_rel;
	tol_sl_abs= tol;
	tol_sl_rel= min(10*tol_sl_abs,0.01);

	if( model_in > 0 )  sl(theta_0, MET );  else  sl(-theta_0, MET );

	tol_sl_abs= tmp1;
	tol_sl_rel= tmp2;

	return;
}



void  Clmbr::sl4R( int met,  double tol,
	double theta_0,  double alpha_0 ) 
{ 
	if(Model==M3)  {
		Rcout << model_msg << endl << endl;
		return;
	}  

	METHOD MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			if(met==3)  MET=MC;  else  {
				stop( methods3msg );
			}
		}
	}

	const double tmp1 = tol_sl_abs, tmp2 = tol_sl_rel;
	tol_sl_abs= tol;
	tol_sl_rel= min(10*tol_sl_abs,0.01);

	if( model_in  > 0 ) 
		sl(theta_0, alpha_0, MET);
	else
		sl(-theta_0, alpha_0, MET);

	tol_sl_abs= tmp1;
	tol_sl_rel= tmp2;

	return;
}



double  Clmbr::sl5R( int met,  int verboseR,  int valueR,
	double tol,  double theta_0 ) 
{ 
	METHOD MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			if(met==3)  MET=MC;  else  {
				stop( methods3msg );
			}
		}
	}

	const bool  verbose = static_cast<bool>( verboseR );
	const bool  value = static_cast<bool>( valueR );
	if( !value )  stop( "dummy argument for dispatch, should be TRUE" );
	
	const double tmp1 = tol_sl_abs, tmp2 = tol_sl_rel;
	tol_sl_abs= tol;
	tol_sl_rel= min(10*tol_sl_abs,0.01);

	double result;
	if( model_in  > 0 ) 
		result= sl(theta_0, MET, verbose);
	else
		result= sl(-theta_0, MET, verbose);

	tol_sl_abs= tmp1;
	tol_sl_rel= tmp2;

	return  result;
}




double  Clmbr::sl6R( int met,  int verboseR,  int valueR,
	double tol,  double theta_0,  double alpha_0 ) 
{ 
	if(Model==M3)  {
		Rcout << model_msg << endl << endl;
		return  NA;
	}  

	METHOD MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			if(met==3)  MET=MC;  else  {
				stop( methods3msg );
			}
		}
	}

	const bool  verbose = static_cast<bool>( verboseR );
	const bool  value = static_cast<bool>( valueR );
	if( !value )  
		stop( "dummy argument for dispatch, should be TRUE" );

	const double tmp1 = tol_sl_abs, tmp2 = tol_sl_rel;
	tol_sl_abs= tol;
	tol_sl_rel= min(10*tol_sl_abs,0.01);

	double result;
	if( model_in  > 0 ) 
		result= sl(theta_0, alpha_0, MET, verbose);
	else
		result= sl(-theta_0, alpha_0, MET, verbose);

	tol_sl_abs= tmp1;
	tol_sl_rel= tmp2;

	return  result;
}



void  Clmbr::ciR( double CL,  int met) 
{ 
	if(CL <=0. || CL >=1.)  stop( CLmsg );

	METHOD  MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			stop( methods2msg );
		}
	}

	const double tmp = SL;
	set_SL(1.-CL);
	ci(MET); 
	set_SL(tmp);

	return; 
}




void  Clmbr::cr3R( double CL,  int met,  double incr) 
{ 
	if(Model==M3)  {
		Rcout << model_msg << endl << endl;
		return;
	}  

	if(CL <=0. || CL >=1.) stop( CLmsg );

	METHOD MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			stop( methods2msg );
		}
	}

	const double tmp = SL;
	set_SL(1.-CL);

	cr( MET, incr );

	set_SL(tmp);

	return;
}



NumericMatrix  Clmbr::cr4R( double CL,  int met, 
	 double incr,  int verboseR ) 
{ 
	if(Model==M3)  {
		Rcout << model_msg << endl << endl;
		return  NumericMatrix(0,0);
	}  

	if(CL <=0. || CL >=1.)  stop( CLmsg );
	const double tmp = SL;
	set_SL(1.-CL);

	METHOD MET;
	if(met==1)  MET=GEO;  else  {
		if(met==2)  MET=AF;  else  { 
			stop( methods2msg );
		}
	}

	double inc;
	if( incr == -1 )  inc= xinc;  else  inc= incr;

	const double  maxwidth = xs[ns-1] - xs[0] + 2;						
	const int  Nmax = maxwidth/inc + ns + 3;

	double*  Btmp= Calloc( Nmax*3, double );

	const bool  verbose = static_cast<bool>( verboseR );
	if( verbose )  
		stop( "dummy argument for dispatch, should be FALSE" );


	const int  nrows = cr( MET, incr, false, Btmp ); 


	set_SL(tmp);

	NumericMatrix  bds( nrows, 3 );
	for(int i=0;i<nrows;i++)  {
		bds(i,0) = *(Btmp + 0*nrows + i);
		bds(i,1) = *(Btmp + 1*nrows + i);
		bds(i,2) = *(Btmp + 2*nrows + i);
	}
	Free( Btmp );

	return bds;
}




void  Clmbr::MLE( void )  const
{ 
	mle(); 
	return; 
}



NumericVector  Clmbr::PARAM( void )  const
// function to pass parameter values to R-code
// internal, not meant for the user
{ 
	double  *const  pdummy =NULL,  *par= Calloc( 5, double );

	mle( false, pdummy, par ); 

	const double  th= par[0],  a= par[1],  b= par[2],  
			bp= par[3],  thfmin= par[4];

	Free( par );

	const double  syc = static_cast<double>( sety_called );

	return  NumericVector::create( th, a, b, bp, thfmin, syc ); 
}



void  Clmbr::SET_rWy( NumericVector rWy )  
{
	const int yn =rWy.size();
	if(yn!=n) stop( _("'rWy' vector has wrong dimension") );

	double*  Ytmp= Calloc( n, double );
	
	for (int i=0;i<n;i++) Ytmp[i] = rWy[i];

	set_sy( Ytmp, GEO2 );

	Free( Ytmp );

	sety_called = true;

	return;
}


