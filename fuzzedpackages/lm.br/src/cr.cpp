//

#include "lmbr.h"



int Clmbr::cr( METHOD met, double incr, bool verbose, double * bounds)
//
//      Checks whether { (theta,alpha)  such that  sig. level  > 'SL' }  is contiguous.
// Returns N = number of rows for the array 'bounds[N][3]' of (th, min_alpha, max_alpha) values.
// If (met==GEO) using Knowles, Siegmund and Zhang's geometric formula to calculate  sig. levels;
// if (met==AF) using Approximate-F method.
//
//	     Uses the 'ci' subroutine to check the 2-parameter significance levels along 
//	the  ( theta0, alphaMLE(theta0) )  ridge,  then finds the alpha-boundaries
//	at increments of theta0.  SL has a single alpha-maxima for a given theta0, 
//  in the broken line model.
//
{
	Rcpp::Function  getOption("getOption");
	Rcpp::Function  Rflush("flush.console");

	bool  verbose_progress;
	if( met==GEO  && !trivial )  verbose_progress = true;  else  verbose_progress = false;


// get theta-boundaries of confidence region(s)
	
	double inc;
	if( incr == -1 )  inc= xinc;  else  inc= incr;

	if(verbose_progress) { Rcout << "   " << _("getting theta-boundaries...   ");  Rflush(); }

	double*  tmp = Calloc( 2*ns, double );

	int numr;
	if(met==GEO)  numr = ci(GEO2,inc,false,tmp);  else  numr = ci(AF2,inc,false,tmp);

	double*  th_bds= Calloc( 2*numr, double );

	int i;
	for (i=0;i<2*numr;i++)  th_bds[i] = tmp[i];
	Free( tmp );


// get (theta,alpha)-boundaries of confidence region(s)
// store boundary values in an  N x 3  array

	const double  th_min= xs[0]-1,  th_max= xs[ns-1]+1;
	if(th_bds[0]== -Inf)  th_bds[0]= th_min;
	if(th_bds[2*numr-1]== +Inf)  th_bds[2*numr-1]= th_max;
	for(i=0;i<2*numr;i++) {
		if( fabs(th_bds[i]-xs[0]) < zero_eq )  th_bds[i]= xs[0];
		if( fabs(th_bds[i]-xs[ns-1]) < zero_eq )  th_bds[i]= xs[ns-1];
	}

	int  N=0;
	double width= 0;						
	for (i=0;i<numr;i++) width += th_bds[2*i+1] - th_bds[2*i];
	const int Nmax = width/inc + 1 + ns + 2*numr + 2;

	double*  bds= Calloc( 3*Nmax, double );


	const double thmle= mle(false);

	if( trivial && !ISNAN(thmle) && thmle!=xs[0] )  {

		const double amle= ahigh(AF,thmle);
		*(bds+N*3+0) = thmle; *(bds+N*3+1) = amle; *(bds+N*3+2) = amle; N++;


	}  else  {

		double  th= 0;

		int  width =0,  col =0;
		double  min_th=0,  max_th=0;
		double  tstart= time(NULL), tfinish= tstart, elapsed= 0;
		if(verbose_progress)  { 
			Rcout << endl << "   " << _("getting alpha-boundaries...   ");  Rflush();
			if(Model==M1)  min_th = max(th_bds[0],xs[0]);  else  min_th = th_bds[0];
			max_th = min(th_bds[2*numr-1],xs[ns-1]);
			Rcpp::IntegerVector  tw = getOption("width");
			width = tw[0];
			col = 33;
		}


		for (i=0;i<numr;i++) {

			const double  tha = th_bds[2*i], thb = th_bds[2*i+1];

// confidence region starts at "tha"  with an open-end if tha = x(1)-1,
// a vertical line if tha= x(n), a vertical line if  tha= x(1) in M1,
// otherwise an alpha-MLE point

			double high_a;
			if( tha == th_min ) {
				*(bds+N*3+0) =tha; *(bds+N*3+1) =a_sl(met,tha,-1); *(bds+N*3+2) =a_sl(met,tha,1); N++; 
			}
			else  if( tha == xs[ns-1] ) {
				th = xs[ns-1] + fabs(rel_print_eps*xs[ns-1]);
				high_a = ahigh(met,th);
				*(bds+N*3+0) =tha; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++; 
				*(bds+N*3+0) =th;  *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++; 
			}
			else  if( Model==M1  &&  tha == xs[0] ) {
				th = xs[0] + fabs(rel_print_eps*xs[0]);
				high_a = ahigh(met,th);
				*(bds+N*3+0) =tha; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
				*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++; 
			}
			else {
				high_a = ahigh(met,tha);
				*(bds+N*3+0) =tha; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
			}


// find (theta,alpha)-boundaries at increments of theta and at datapoints

			int ka= -1, kb= -1; 
			while (xs[ka+1]<=tha && ka<ns-1) ka++;
			while (xs[kb+1]<=thb && kb<ns-1) kb++;


			for (int k=ka;k<=kb;k++) {

				double th1i; 
				if( k < 0 ) th1i= tha; else th1i= max(tha, xs[k]);


// confidence region boundary is discontinuous at x(n), and at x(1) in M1,

				if( tha < th1i  &&  th1i < thb )  {

					if( th1i==xs[ns-1] ) {
						th = th1i - fabs(rel_print_eps*th1i);
						*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
					}

					*(bds+N*3+0) =th1i; *(bds+N*3+1) =a_sl(met,th1i,-1); *(bds+N*3+2) =a_sl(met,th1i,1); N++;

					if( Model==M1  &&  th1i==xs[0] ) {
						  th = th1i + fabs(rel_print_eps*th1i);
						  *(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
					}
				}


				double th1, th2;
				th1= floor( th1i );
				while ( th1 < th1i + tol_xb )  th1 +=inc;
				if(k< ns-1) th2= min(xs[k+1],thb); else th2= thb;


				for (th = th1; th < th2 - tol_xb; th += inc) {

					if( fabs(th)<zero_eq) th=0.;

					*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++; 

					if(verbose_progress)  { 
						tfinish = time( NULL );
						elapsed = tfinish - tstart;
						if( elapsed > 10 ) {
							if( col > width - 6 )  { Rcout << endl;  col= 0; } 
							const double  progress = floor( 100*(th-min_th)/(max_th-min_th) );
							Rcout << progress << "%...   ";   Rflush();
							tstart= tfinish;
							col += 9;
						}
					}

				}

			}


// confidence region ends at "thb" with an open end if  thb = x(n) + 1,
// a vertical line if thb= x(n), a vertical line if  thb= x(1) in M1,
// otherwise an 'alpha-MLE' point

			if( thb == th_max ) {
				*(bds+N*3+0) =thb; *(bds+N*3+1) =a_sl(met,thb,-1); *(bds+N*3+2) =a_sl(met,thb,1); N++; 
			}
			else  if( thb == xs[ns-1] ) {
				th = xs[ns-1] - fabs(rel_print_eps*xs[ns-1]);
				high_a = ahigh(met,th);
				*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
				*(bds+N*3+0) =thb; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
			}
			else  if( Model==M1  &&  thb == xs[0] ) {
				th = xs[0] - fabs(rel_print_eps*xs[0]);
				high_a = ahigh(met,th);
				*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
				*(bds+N*3+0) =thb; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
			}
			else {
				high_a = ahigh(met,thb);
				*(bds+N*3+0) =thb; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
			}

		}

	}


	if(verbose_progress)  Rcout << endl;

	if(verbose)  {
		Rcout << endl << 100*(1-SL) << _("-percent joint confidence region for  (theta, alpha)  by ");
		if(met==GEO)  Rcout << "CLR" << endl; else  Rcout << "AF" << endl;
		Rcout << endl << setw(10) << "theta" << setw(16) << "min. alpha" << setw(15) 
				<< "max. alpha" << endl << endl;
		if( model_in > 0 ) {
			for(i=0;i<N;i++)  {
				Rcout << setw(10) << *(bds+i*3+0) << "," << setw(14) << *(bds+i*3+1)
					<< "," << setw(14) << *(bds+i*3+2) << endl;
				for(int j=0;j<numr;j++)  if( *(bds+i*3+0) == th_bds[2*j+1] )  Rcout << endl;
			}
		}  else  {
			for(i=N-1;i>=0;i--)  {
				Rcout << setw(10) << *(bds+i*3+0)*(-1.) << "," << setw(14) << *(bds+i*3+1)
					<< "," << setw(14) << *(bds+i*3+2) << endl;
				for(int j=0;j<numr;j++)  if( *(bds+i*3+0) == th_bds[2*j] )  Rcout << endl;
			}
		}
		Rflush(); 
	}


	if(bounds!=0)  {
		for(i=0;i<N;i++) {
			if( model_in > 0 ) {
				for(int j=0;j<3;j++)  *(bounds+j*N + i) = *(bds+i*3+j);
			} else {
				*(bounds + i) = *(bds+(N-1-i)*3+0)*(-1.);
				for(int j=1;j<3;j++)  *(bounds+j*N + i) = *(bds+(N-1-i)*3+j);
			}
		}
	}

	Free( th_bds );  Free( bds );
	return N;
}




double Clmbr::a_sl( METHOD met, double th, int high_low)
// return alpha high or low value for a given SL value and a given theta0
// if high_low > 0 return high value, if high_low < 0 return low value
// if met=GEO use AF alpha as first guess, then grid search + bisection
{
	if(trivial)  return ahigh(met,th);
	if (met==AF)  {
		return a_af(th, high_low);  
	}  else  {
		if(th!=old_th) { ah =ahigh(met,th);  old_th =th; }
		if( sl(th,ah,met,false) < SL )  stop( _("'a_sl' initial point below critical SL") );
		const double incr = inc_y*high_low; 
		double guess = a_af(th, high_low);
		if ( ISNAN(guess) || (guess-ah)*high_low < zero_eq) guess = ah+incr;
		double guess2 = ah;
		while ( sl(th,guess,met,false) > SL) { 
			guess2 = guess;  
			guess += incr; 
		}
		const double a_geo =bisect(guess,guess2,&Clmbr::sl_a,1,SL,-tol_yb);
		return a_geo;
	}
}



double Clmbr::a_af( double th, int high_low )
// return alpha high or low value for a given SL and given theta0 by AF
// alpha-boundaries by AF are roots of a quadratic
// if 'high_low' > 0 return high value, if 'high_low' < 0 return low value
{

	if (th != prev_th) {
		prev_th = th;

		bool th_ex;
		if ( (Model==M1 && th<=xs[0]) || xs[ns-1]<=th) th_ex=true; else th_ex=false;

		double  a, bp, c, xc;
		if (th_ex) {
			if (variance_unknown) xc=x_vu_ex; else xc=x_vk_ex;
			if(Model==M1) {
				const double ff = sxx-2*sx1*th+s11*th*th, fy = yx-th*y1, f1 = sx1-th*s11;
				a =s11-f1*f1/ff; bp =y1-fy*f1/ff; c =sysq-fy*fy/ff -xc;
			} else {
				a= s11; bp= y1; c= sysq - xc;
			}
		} else {					
			if (variance_unknown) xc=x_vu; else xc=x_vk;
			int k=0; while (k<ns && xs[k]<th) k++;
			const double  fry = *psy*gfr(th,k), fr1 = *psig1*gfr(th,k);
			if(Model==M1) {
				const double smy = *psy*gsm(th,k), sm1=*psig1*gsm(th,k); 
				a =s11-sm1*sm1-fr1*fr1; bp =y1-smy*sm1-fry*fr1; c =sysq-smy*smy-fry*fry -xc;
			} else {
				a= s11-fr1*fr1; bp= y1-fry*fr1; c= sysq - fry*fry - xc;
			}
		}
		double rad = bp*bp-a*c, rrad; 
		if(fabs(rad)<zero_eq) rad=0.;
		if(rad<0) rrad=NaN; else rrad= sqrt(rad);
		a_low = (bp-rrad)/a;
		a_high = (bp+rrad)/a;
	}

	if (high_low <0) return a_low; else return a_high;
}



double Clmbr::ahigh( METHOD met, double th )
// return 'alpha' value that gives the highest significance level for a given theta
{
	if(trivial) {
		const double  slope= ( (*py)[ is[1] ] - (*py)[ is[0] ] )/( xs[1]-xs[0]),  intercept= (*py)[ is[0] ]  - slope*xs[0],
			amle = slope*th + intercept;
		return amle;
	}

	bool th_ex;
	if ( (Model==M1 && th<=xs[0]) || xs[ns-1]<=th) th_ex=true; else th_ex=false;

	double amle;
	if (th_ex) {
		if(Model==M1) {
			const double ff = sxx-2*sx1*th+s11*th*th, fy = yx-th*y1, f1 = sx1-th*s11;
			amle = y1/s11 - f1*(fy-f1*y1/s11)/(ff*s11-f1*f1);
		} else   amle = y1/s11;
	} else {
		int k=0; while (k<ns && xs[k]<th) k++;
		const double fr1 = *psig1*gfr(th,k), fry = *psy*gfr(th,k);
		if(Model==M1) {
			const double sm1 = *psig1*gsm(th,k), smy = *psy*gsm(th,k);
			amle = (y1-fry*fr1-smy*sm1)/(s11-fr1*fr1-sm1*sm1);
		} else   amle = (y1 - fr1*fry)/(s11-fr1*fr1);
	}
	double high_alpha = amle;


	if ( met==GEO  &&  !th_ex ) {

		if(th!=th0) set_theta0(th,GEO2);

// find maxima of sl_geo2 by grid search, starting from amle
		double incr=inc_y, ta=amle; 
		set_alpha0(ta);
		double tsl1=sl_geo2();
		ta +=incr;
		set_alpha0(ta);
		double tsl2=sl_geo2();

		if(tsl1 > tsl2) {
			incr = -incr;
			ta = amle;
			tsl2 = tsl1;
			tsl1 = tsl2-1.;
		}


		if(tsl2 < SL) {

			while( fabs(incr) > tol_yb/2 ) {

				while (tsl2 > tsl1) {
					tsl1 = tsl2;
					ta += incr;
					set_alpha0(ta);
					tsl2=sl_geo2();
				}
				incr = -incr/2;
				if(tsl1 < SL/64) break;		// far from boundary, not worthwhile
				tsl1 = tsl2 - 1.;
			}

			ta -= incr*10;
		}

		high_alpha = ta;
	}

	return high_alpha;
}



double Clmbr::sl_a( double alpha, int k )
// wrapper for use in the bisection routine
// return sl_geo2(th0,alpha)
// assume th0 preset and use default tolerance
{
	set_alpha0(alpha, GEO2);
	return sl_geo2();
}

