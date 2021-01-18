//
//  Newton_raphson.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 24/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//


//#include "Ulam.hpp"
#define MAXIT 200
#define ALF 1.0e-4 
#define TOLX 1.0e-7
#define MAXITS 200
#define EPS2 1.0e-4
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define STPMX 100.0

#include "Newton_raphson.h"
#include <iostream>
#include <math.h>
#include "Generic.h"
#include <cmath>
#include <float.h>
#include <cstring>
#include <math.h>
#include "Hamming.h"


#define FREERETURN {free_matrix(fjac,1,n,1,n);free_vector(fvec,1,n);\
free_vector(p,1,n);free_ivector(indx,1,n);return;}
/*
 * Estimates theta parameters by means of running newton iterative algorithm.
 */
double Newton_raphson::Newton_raphson_method( double dAvg_val, double initialGuess, int distanceModel_val, int model, int j_index, long double*count) {
    count_ = count;//count of num permus for ulam
    double xacc = 0.000001;
//    n_ = n_val;
    j_index_ = j_index; // for kendall GMM
    model_ =  model;
    dist_avg_=dAvg_val;
    UPPER_THETA=5;
    distance_id_=distanceModel_val;
    double theta;
    //TEST //for (double i = -5.1 ; i< 5.1 ; i++) cout<<"Theta "<<i<< " f: "<<f(i)<<" dev: "<<fdev(i)<<endl;
    theta= rtsafe(initialGuess, UPPER_THETA, xacc);
    //cout<<"newton "<<j_index_<<" "<<dist_avg_<<" "<<theta<<endl;
    return theta;
}

/*
 * Calculates Newton algorithms f and fdev functions
 */
void Newton_raphson::funcd(double theta, double *ff, double *ffdev) {
    //*ff = f(theta, n, Vjs);
    //*ffdev = fdev(theta, n);
    *ff=f(theta);
    *ffdev=fdev(theta);
}

/*
 * Newton - Rapshon execution algorithm.
 */
double Newton_raphson::rtsafe(double x1, double x2, double xacc) {
    int j;
    double dx, dxold;
    double  temp, xh, xl, rts;
    double f, df, fl, fh;
    
    funcd(x1, &fl, &df);//params,: theta, f, fdev
    funcd(x2, &fh, &df);
    //if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    //	cout<<"Root must be bracketed in rtsafe"<<endl;
    if (fl == 0.0)
        return x1;
    if (fh == 0.0)
        return x2;
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2);
    rts = x1; //<-fijamos un valor inicial para el theta.
    dxold = fabs(x2 - x1);
    dx = dxold;
    funcd(rts, &f, &df);
    for (j = 1; j <= MAXIT; j++) {
        //cout<<"rts: "<<rts<<". Function val: "<<f<<" f_dev "<<df<<endl;
        //Initialize the guess for root, the “stepsize before last,” and the last step.
        //Loop over allowed iterations.
        if ((((rts - xh) * df - f)*((rts - xl) * df - f) > 0.0)   || (fabs(2.0 * f) > fabs(dxold * df))) {
            //Bisect if Newton out of range, //or not decreasing fast enough.
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts)
                return rts;
        }
        else {
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts)
                return rts;
            //Change in root is negligible. Newton step acceptable. Take it.
        }
        //cout<<"DX: "<<dx<<endl;
        if (fabs(dx) < xacc) {
            return rts;
        }
        funcd(rts, &f, &df); //The one new function evaluation per iteration.
        //Orient the search so that f(xl) < 0.
        //Convergence criterion.
        
        if (f < 0.0) //Maintain the bracket on the root.
            xl = rts;
        else
            xh = rts;
        
        //cout<<"j_index: "<<j_index_<<" dist_avg_: "<<dist_avg_<<" rts: "<<rts<<" Function val: "<<f<<" f_dev "<<df<<endl;
    }
  //  cout << "Maximum number of iterations exceeded in rtsafe" << endl;
    
    return 0.0; //Never get here.
}


// * Theta parameter estimation function.
double Newton_raphson::f(double theta) {
    if(distance_id_== CAYLEY_DISTANCE){//Cayley distance
        double sum=0;
        for(int j = 1 ; j < n_ ; j++){
            double ex = exp(theta);
            double denom = j+ex;
            sum += (double)j/denom;
            //sum += (double) j /(double)(exp(theta)+(double)j);
        }
        return (double)(sum - (double)dist_avg_);
    }else if (distance_id_ == KENDALL_DISTANCE && model_ == MALLOWS_MODEL){//kendall MM
        double aux= 0;
        for(int j = 1 ; j < n_ ; j++){
            int k = n_ - j + 1;
            aux += (k * exp(-theta * k))/(1 - exp(-theta * k));
        }
        double aux2 = (n_-1) / (exp( theta ) - 1) - dist_avg_;
        //cout<<"trace theta :"<<theta<<" fval "<<aux2-aux<<endl;
        return aux2 - aux;
    }else if (distance_id_ == KENDALL_DISTANCE && model_ == GENERALIZED_MALLOWS_MODEL){//kendall GMM
        // oper = oper + oper1 + (pow(n - j + 1, 2) * exp(theta * (n - j + 1))) / (double) pow(exp(theta * (n - j + 1)) - 1, 2);
        int k = n_ - j_index_ + 1;
        double a = k /( exp(theta * k ) - 1 );
        double b = 1 / (exp( theta ) - 1 );
        //cout<<"trace theta :"<<theta<<" fval "<<dist_avg_<<" "<<a<<" "<<b<<endl;
        return dist_avg_ + a - b;
    }else if(distance_id_ == ULAM_DISTANCE){
        double numer = 0, denom = 0;
        for (int d = 0 ; d < n_ - 1; d++){
            double aux = count_[d ] * exp(-theta *d ) ;
            numer += aux * d;
            denom += aux;
        }
        return numer / denom - dist_avg_;
    }else if (distance_id_ == HAMMING_DISTANCE && model_ == MALLOWS_MODEL ){
        Generic gen;
        long double     x_j = 0 , sum_to_n = 0 , sum_to_n_1 = 0 , psi = 0 , psi_der = 0;
        for ( int j = 0 ; j <= n_ ; j++ ){
            x_j = (long double) pow(exp(theta)-1, j) / facts_[ j ];
            if ( j < n_ ) sum_to_n_1 += x_j;
            sum_to_n += x_j;
            if(sum_to_n> DBL_MAX ||sum_to_n != sum_to_n )
                return DBL_MAX;
        }
        psi = facts_[ n_ ] * exp(-theta * n_ ) * sum_to_n;
        psi_der = - n_ * psi + facts_[ n_ ] * exp( theta * ( 1 - n_ )) * sum_to_n_1;
        double f_fligner = (double)  (psi_der / psi + dist_avg_);
        //cout<<"theta "<<theta<<" dist "<<dist_avg_<<" fu "<<f_fligner<<endl;
        if(f_fligner != f_fligner )
            f_fligner=0;//trace
        return f_fligner;
    }else if (distance_id_ == HAMMING_DISTANCE && model_ == GENERALIZED_MALLOWS_MODEL ){
        //cout<<"Solve Weigthed Hamming Mallows with the multivariate newthon Raphson, mnewt "<<endl;
        //exit(1);
    }
    return 0;
}

// Theta parameter estimation function derivation.
double Newton_raphson::fdev(double theta) {
    if(distance_id_ == CAYLEY_DISTANCE ){
        double sum=0;
        for(int j = 1 ; j < n_ ; j++)
            sum += (double)( - j * exp( theta ))/pow(exp(theta) + j, 2);
        return sum;
    }else if ( distance_id_ == KENDALL_DISTANCE &&  model_ == MALLOWS_MODEL ){//kendall mm
        double aux= 0;
        for(int j = 1 ; j < n_ ; j++){
            int k = n_ - j + 1;
            aux += (k * k * exp( -theta * k ))/pow((1 - exp(-theta * k)) , 2 );
        }
        double aux2 = (- n_ + 1) * exp( theta ) / pow ((exp( theta ) - 1) , 2 );
        return aux2 + aux;
    }else if ( distance_id_ == KENDALL_DISTANCE && model_ == GENERALIZED_MALLOWS_MODEL){ //kendall gmm
        int k = n_ - j_index_ + 1;
        double a = - k * k * exp(theta * k )  /pow ((exp(theta * k ) - 1), 2);
        double b = exp(theta) / pow ( (exp( theta ) + 1 ), 2 );
        return a + b;
    }else if(distance_id_ ==  ULAM_DISTANCE){
        long double numer1 = 0, numer2 = 0, numer3 = 0, denom = 0;
        for (int d = 0 ; d < n_ - 1; d ++){
            long double aux = count_[d] * exp(-theta *d ) ;
            numer1 += aux * d * d ;
            numer2 += aux;
            numer3 += aux * d;
            denom += aux;
        }
        return (-numer1 * numer2 - numer3*numer3)/(denom*denom);
    }else if (distance_id_ == HAMMING_DISTANCE && model_ == MALLOWS_MODEL ){
        Generic gen;
        long double     x_j = 0 , sum_to_n = 0 , sum_to_n_1 = 0 , sum_to_n_2 = 0 , psi=0, psi_der = 0 , psi_der_2 = 0 ;

        for ( int j = 0 ; j <= n_ ; j++ ){
            x_j = (long double) pow(exp(theta)-1, j) / facts_[ j ];
            if ( j < n_ -1 ) sum_to_n_2 += x_j;
            if ( j < n_ ) sum_to_n_1 += x_j;
            sum_to_n += x_j;
            if(sum_to_n_1 > DBL_MAX || sum_to_n_1 != sum_to_n_1)
                return DBL_MAX-1;
        }
        psi = facts_[ n_ ] * exp(-theta * n_ ) * sum_to_n;
        psi_der = - n_ * psi + facts_[ n_ ] * exp(theta*( 1 - n_ )) * sum_to_n_1;
        psi_der_2 = -n_ * psi_der + facts_[ n_ ] * exp(theta *(1 - n_ )) * ((1 - n_ ) * sum_to_n_1 + sum_to_n_2 );
        double res = (double) ( - psi_der_2 * psi_der - psi_der * psi_der )/ (psi*psi);
        //cout<<"theta "<<theta<<" dist "<<dist_avg_<<"ps... "<<psi<<" "<<psi_der<<" "<<psi_der_2<<" fd "<<res<<endl;
        if(res > DBL_MAX)
            return DBL_MAX-1;
        return res;
    }else if (distance_id_ == HAMMING_DISTANCE && model_ == GENERALIZED_MALLOWS_MODEL ){
        //cout<<"Solve Weigthed Hamming Mallows with Newton_raphson::mle_theta_weighted_mallows_hamming "<<endl;
        //exit(1);
    }
    return 0;
}


void Newton_raphson::mle_theta_weighted_mallows_hamming(int m, double*h_avg , double * theta){
    //OJO : Given a starting point p[1..n]!!!
    //minimize a function on multiple variables
    Generic gen;
    m_          = m;
    h_avg_      = h_avg;
    int         iter  = 0 ;
    double      fmin  = 0 ;
    double      * point  = new double[ n_ + 1 ];
    double      * dpoint = new double[ n_ + 1 ];
    init_optim_wmh();
    for (int i = 0 ; i < n_; i++)point[ i + 1 ] = 0.2;
    for ( int j = 0 ; j < n_ ; j ++){
       // if (h_avg_[ j ] == 0) h_avg_[ j ] = 0.001;//los ceros dan problemas, txapu
       // if (h_avg_[ j ] >= 1) h_avg_[ j ] =  0.99;
    }
    frprmn(point, n_ , 0.0001, &iter , &fmin , &Newton_raphson::likeli_wmh, &Newton_raphson::dlikeli_wmh);
    for ( int j = 0 ; j < n_ ; j ++) theta[j] = (double) point [ j + 1 ];
    delete [] point;
    delete [] dpoint;
}

double Newton_raphson::likeli_wmh(double x[] ){
    //x is a vector from 1..n
    //both  likeli_wmh and dlikeli_wmh return the result *(-1) because frprmn is for minimization and we need maximization
    Generic gen;
    long double   like = 0, aux1 = 0 , aux2 = 0, sum_theta = 0;
    bool penalty = false;
    double  * do_x = new double[ n_ ];
    for (int i = 0 ; i < n_; i++) {
        do_x[i]  = (long double)x[i + 1];
        sum_theta += x[i + 1];
        if (do_x [ i ] < 0 || do_x [ i ] > 10  )
            penalty = true;
    }
    gen.elementary_symmetric_polynomial(do_x, n_,t_, aux_esp_, esp_ );
    for (int i = 0 ; i < n_ ; i++)   aux1 += do_x[ i ] * h_avg_[ i ];
    for (int k = 0 ; k <= n_ ; k++)  aux2 += facts_[ n_ - k ] * esp_[ k ];
    aux2 = aux2 * exp( -sum_theta ); //psi
    like = - m_ * ( aux1 + log( aux2 ));
    
    delete [] do_x;
    if (like != like || penalty ){//|| penalty != 0
        return DBL_MAX;
    }
    return - like ;//* (penalty+1);
}


void Newton_raphson::dlikeli_wmh(double x[] , double deriv[] ){
    //derivada de la verosimilitud
    //x y deriv [1..n]
    Generic gen;
    long double   psi = 0  ;
    double      * do_x      = new double[ n_ ];////el x es un puto vector de 1..n WTF!!!!!!!!!!
    long double * psi_der   = new long double [ n_ ];
    double        sum_theta = 0, aux = 0 ;
    for (int i = 0 ; i < n_; i++)  {
        do_x[i]  = ( double)x[i + 1];
        sum_theta += x[ i+1 ];
    }
    gen.elementary_symmetric_polynomial(do_x, n_, t_, aux_esp_,esp_ );
    gen.split_elementary_symmetric_polynomial (esp_, do_x , n_, esp_no_a_, esp_yes_a_);
    psi = 0 ;
    for (int k = 0 ; k <= n_ ; k++)
        psi += facts_[ n_ - k ] * esp_[ k ];
    psi = psi * exp( -sum_theta ); //psi
    for (int i = 0 ; i < n_ ; i++){
        psi_der[ i ] = 0;
        aux = 0 ;
        for (int k = 1 ; k <= n_ ; k++)
            aux += facts_[ n_ - k ] * esp_no_a_[ k - 1][ i ];
        aux = aux * exp( -sum_theta + do_x[ i ]); //psi
        psi_der[ i ] = - psi + aux;
        //both  likeli_wmh and dlikeli_wmh return the result *(-1) because frprmn is for minimization and we need maximization
        deriv[ i + 1 ] = ( -1 )  * (psi_der[ i ] / psi  + (long double) h_avg_[ i ]) ;
    }
    delete [] do_x;
//    delete [] removed;
    delete [] psi_der;
}


#include <math.h>
#define ITMAX 100
#define EPS1 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void Newton_raphson::frprmn(double p[], int n, double ftol, int *iter, double *fret, double (Newton_raphson::*func)(double []), void (Newton_raphson::*dfunc)(double [], double []))
{
    //Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func, using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is input as ftol. Returned quantities are p (the location of the minimum), iter (the number of iterations that were performed), and fret (the minimum value of the function). The routine linmin is called to perform line minimizations.
	int j,its;
	long double gg,gam,fp,dgg;
	double *g,*h,*xi;//,*vector();
	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*this.*func)(p);
	(*this.*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
//cout<<"p[i]: ";for (j=1;j<=n;j++) cout <<p[j]<<" ";cout<<" point from frprmn"<<endl;
//cout<<"h_avg: ";for (j=1;j<=n;j++) cout <<h_avg_[j-1]<<" ";cout<<" h from frprmn"<<endl;
		*iter=its;
        dlinmin(p,xi,n,fret,func, dfunc);
        //for (int i = 0 ; i < n ; i++) if (p[ i +1 ]<0) p[i+1] = 0 ;
        //linmin(p, xi, n, fret, func);
        if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS1)) {
			FREEALL
			return;
		}
		fp=(*this.*func)(p);
		(*this.*dfunc)(p,xi);
        dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
            //dgg += xi[j]*xi[j];	//or this or the next(see numerical receipes)
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
        for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
    //char str_msg[100];strcpy(str_msg, "Too many iterations in FRPRMN");nrerror(str_msg);
    //ERROR
}

#undef ITMAX 
#undef EPS1 
#undef FREEALL 


double   (Newton_raphson::*nrfunc)(double []);
void    (Newton_raphson::*nrdfun)(double [], double [] );

int ncom=0;
double *pcom=0,*xicom=0;
//void (*nrdfun)(double [], double [] );

#define TOL 2.0e-4


void Newton_raphson::dlinmin(double p[], double xi[], int n, double *fret,
                             double (Newton_raphson::*func) (double []),
                             void  (Newton_raphson::*dfunc)(double [], double []))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	nrdfun=dfunc;
	for (j=1;j<=n;j++) {
		pcom[j] =p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,&Newton_raphson::f1dim);
	*fret=dbrent(ax,xx,bx,&Newton_raphson::f1dim,&Newton_raphson::df1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
        //if (p[ j ] < 0 ) p[j] = 0 ;//test1
        //if (p[ j ] > 10 ) p[j] = 10 ;
	}
	free_vector(xicom,1,n);
	free_vector(pcom, 1,n);
}

#undef TOL

////
#define TOL 2.0e-4

//int ncom=0;	/* defining declarations */
//float *pcom=0,*xicom=0,(*nrfunc)();

//void linmin(p,xi,n,fret,func)
//float p[],xi[],*fret,(*func)();
//int n;
void    Newton_raphson::linmin  (double p [] ,double xi[],int n,double *fret,double (Newton_raphson::*func) (double []))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
//	float brent(),f1dim(),*vector();
//	void mnbrak(),free_vector();
    
	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,&Newton_raphson::f1dim);
	*fret=brent(ax,xx,bx,&Newton_raphson::f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}

#undef TOL

////

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void Newton_raphson::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (Newton_raphson::*func)(double))
{
	double ulim,u,r,q,fu,dum;
    
	*fa=(*this.*func)(*ax);
	*fb=(*this.*func)(*bx);
	/*if (*fb > *fa && *ax == 0 ) {
        *ax = *bx = *cx = 0.1 ;
        return ;
    }*/
    if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*this.*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
        (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*this.*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*this.*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*this.*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*this.*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*this.*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*this.*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT


#include <math.h>

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double Newton_raphson::brent(double ax,double bx,double cx,double (Newton_raphson::*f)(double),double tol,double *xmin)
//double ax,bx,cx,tol,*xmin;
//double (*f)();	/* ANSI: double (*f)(double); */
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;
	  
	a=((ax < cx) ? ax : cx);
	b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*this.*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
                q=fabs(q);
                etemp=e;
                e=d;
                if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                    d=CGOLD*(e=(x >= xm ? a-x : b-x));
                    else {
                        d=p/q;
                        u=x+d;
                        if (u-a < tol2 || b-u < tol2)
                            d=SIGN(tol1,xm-x);
                            }
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*this.*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
                SHFT(v,w,x,u)
                SHFT(fv,fw,fx,fu)
                } else {
                    if (u < x) a=u; else b=u;
                        if (fu <= fw || w == x) {
                            v=w;
                            w=u;
                            fv=fw;
                            fw=fu;
                        } else if (fu <= fv || v == x || v == w) {
                            v=u;
                            fv=fu;
                        }
                }
	}
    //char str_msg[100];strcpy(str_msg,"Too many iterations in BRENT");nrerror(str_msg);
//ERROR
	*xmin=x;
	return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN




#define ITMAX 100
#define ZEPS 1.0e-10 
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double Newton_raphson::dbrent(double ax, double bx, double cx, double (Newton_raphson::*f)(double), double (Newton_raphson::*df)(double), double tol, double *xmin)
{
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
    
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*this.*f)(x);//f1
	dw=dv=dx=(*this.*df)(x);//df1

    
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
       if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx)  d1=(w-x)*dx/(dx-dw);
                if (dv != dx)  d2=(v-x)*dx/(dx-dv);
                    u1=x+d1;
                    u2=x+d2;
                    ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
                    ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
                    olde=e;
                    e=d;
                    if (ok1 || ok2) {
                        if (ok1 && ok2)
                            d=(fabs(d1) < fabs(d2) ? d1 : d2);
                            else if (ok1)
                                d=d1;
                                else
                                    d=d2;
                                    if (fabs(d) <= fabs(0.5*olde)) {
                                        u=x+d;
                                        if (u-a < tol2 || b-u < tol2)
                                            d=SIGN(tol1,xm-x);
                                            } else {
                                                d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
                                            }
                    } else {
                        d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
                    }
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
      	if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*this.*f)(u);
		} else {
			u=x+SIGN(tol1,d);
			fu=(*this.*f)(u);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
      	du=(*this.*df)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
                MOV3(v,fv,dv, w,fw,dw)
                MOV3(w,fw,dw, x,fx,dx)
                MOV3(x,fx,dx, u,fu,du)
                } else {
                    if (u < x) a=u; else b=u;
                        if (fu <= fw || w == x) {
                            MOV3(v,fv,dv, w,fw,dw)
                            MOV3(w,fw,dw, u,fu,du)
                        } else if (fu < fv || v == x || v == w) {
                            MOV3(v,fv,dv, u,fu,du)
                        }
                }
	}
    //ERROR char str_msg[100];strcpy(str_msg, "Too many iterations in routine DBRENT");nrerror(str_msg);
    return 0;
}

#undef ITMAX
#undef ZEPS
#undef SIGN
#undef MOV3

double Newton_raphson::f1dim(double x)
{
	int j;
	double f,*xt;//,*vector();
    
	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    f=(*this.*nrfunc)(xt);
    free_vector(xt,1,ncom);
    return f;
}

double Newton_raphson::df1dim(double x)
{
	int j;
	double df1=0.0;
	double *xt,*df;//,*vector();
	//void free_vector();
    
	xt=vector(1,ncom);
	df=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    (*this.*nrdfun)(xt,df);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
    free_vector(df,1,ncom);
    free_vector(xt,1,ncom);

    return df1;
}

#define NR_END 1
#define FREE_ARG char*
void Newton_raphson::nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	//fprintf(stderr,"Numerical Recipes run-time error...\n");
	//fprintf(stderr,"%s\n",error_text);
	//fprintf(stderr,"...now exiting to system...\n");
	//exit(1);
}

double *Newton_raphson::vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
    
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) {
        //ERRORchar str_msg[100];strcpy(str_msg, "allocation failure in vector()");nrerror(str_msg);
    }
	return v-nl+NR_END;
}

void Newton_raphson::free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
