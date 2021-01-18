/* 
 *  random.cpp : Randomnumbergenerator for several distributions.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: random.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <math.h>
#ifdef _MSC_VER
	#include <ctime>
#else
	#include <sys/time.h>
#endif
#include "error.h"
#include "random.h"

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace std;

class _random_number_initializer
{
public:
	_random_number_initializer();
};

static _random_number_initializer __random_number_initializer;

_random_number_initializer::_random_number_initializer()
{
	// timeval t;	auskommentiert bei Sebastian Ruthe , weil nicht windows kompatibel
	//gettimeofday(&t,0); s.o. 
	//set_random_seed(t.tv_usec);    s.o.
	set_random_seed();
}

void set_random_seed()
{
	//srand48(s); // auskommentiert bei Sebastian Ruthe, weil nicht mit Windows kompatibel
#ifdef _MSC_VER
	srand(time(0));
#else
	GetRNGstate();
#endif
	//srand(int(s));
}

int random(int min,int max)
{
	errif(max < min,"random: illegal bounds " << min << " and " << max);
    //return random()%(max-min+1)+min;		// auskommentiert bei Sebastian Ruthe weil nicht mit Windows kompatibel 
    
#ifdef _MSC_VER
	return (rand() % (max - min + 1) + min);  //AP - debug under VC++
#else
	return ((int)(unif_rand()*32767)%(max-min+1)+min);
#endif
}


double Uniform(double a,double b)
{
	errif(a > b,"Uni: illegal bounds " << a << " and " << b);
	//return drand48()*(b-a)+a; 		// auskommentiert bei Sebastian Ruthe weil nicht mit Windows kompatibel
	return (((double)(unif_rand()*32767)) / 32767) * (b-a) + a; 
}

double Normal(double m,double s)
{
	errif(s <=0,"Normal: illegal variance " << s);
	
	if(m==0.0 && s==1.0)
	{
		double v1,v2,s;
		do
		{
			v1=2*Uniform(0,1)-1;
			v2=2*Uniform(0,1)-1;
			s=v1*v1+v2*v2;
		} while(s >= 1.0);

		return v1*sqrt(-2.0 * log(s) /s);
	}
	else
		return m+sqrt(s)*Normal(0.0,1.0);
}

double Exponential(double l)
{
	errif(l <=0,"Exponential: illegal parameter " << l);

	double u;

	do
	{
		u=Uniform(0,1);
	}
	while(u == 0);
	
	return -(1.0/l) * log(u);
}

double Gamma(double a)
{
	errif(a <=1,"Gamma: illegal parameter " << a);

	double y,x,v;

	do
	{
		y=tan(M_PI*Uniform(0,1));
		x=sqrt(2*a-1)*y+a-1;
		v=Uniform(0,1);
	}
	while(x <= 0.0 || v > (1+y*y) * exp( (a-1)*log(x/(a-1)) - sqrt(2*a-1)*y ) );

	return x;
}

double Chi2(int n)
{
	errif(n <=0,"Chi2: degrees of freedom " << n);

	return 2*Gamma(double(n)/2);
}

double Fischer(int m,int n)
{
	errif(m <=0||n <=0,"Fischer: illegal parameters " << m << "," << n);
	
	return Chi2(m)*double(n) / (Chi2(n)*double(m));
}

double gamma_2(int n)
{
	errif(n <= 0,"gamma_2: non-positive argument n=" << n);
	
	if(n % 2)
	{
		n=(n-1)/2;
		int prd=1;
		for(int i=1; i<=n; i++)
			prd*=(2*i-1);

		errif(prd < 0,"gamma_2: overflow");
		
		return double(prd)*sqrt(M_PI)/pow(2.0,n);
	}
	else
	{
		n=n/2 - 1;
		int prd=1;
		for(int i=2; i<=n; i++)
			prd*=i;

		errif(prd < 0,"gamma_2: overflow");

		return prd;
	}
}

double dfChi2(int n,double x)
{
	errif(x <= 0.0,"dfChi2: non-positive argument x=" << x);
	errif(n <= 0,"dfChi2: non-positive argument n=" << n);

	double n_2=double(n)/2.0;

 	return ( 1.0/(pow(2.0,n_2) * gamma_2(n)) )
 	  * pow(x,n_2-1.0)
 	  * exp(-x/2.0);
}

double limitChi2(int n,double p)
{
	errif(n <= 0,"limitChi2: non-positive argument n=" << n);
	errif(p <= 0.0 || p >= 1.0,"limitChi2: invalid percentage " << p);

	static double Chi2_table90[20]=
	{2.706,4.605,6.251,7.779,9.2236,
		 10.645,12.017,13.362,14.684,15.987,
		 17.275,18.549,19.218,21.064,22.307,
		 23.542,24.769,25.989,27.204,28.412};
	static double Chi2_table95[20]=
	{3.841,5.991,7.815,9.488,11.070,
		 12.592,14.067,15.507,16.919,18.307,
		 19.675,21.026,22.362,23.685,24.996,
		 26.296,27.587,28.869,30.144,31.410};
	static double Chi2_table99[20]=
	{6.635,9.210,11.345,13.277,15.086,
		 16.812,18.475,20.090,21.666,23.209,
		 24.725,26.217,27.688,29.141,30.578,
		 32.000,33.409,34.805,36.191,37.566};
	static double Chi2_table999[25]= // Hups
	{10.827,13.815,16.268,18.465,20.517,
	 22.457,24.322,26.125,27.877,29.588,
	 31.264,32.909,34.528,36.123,37.697,
	 39.252,40.790,42.312,43.280,45.315,
	 46.797,48.268,49.728,51.179,52.620};

	errif(n > 20,"limitChi2: unimplemented Chi2(" << n
	  << ") " << p*100 << "%");

	if(p==0.90)
		return Chi2_table90[n-1];
	else if(p==0.95)
		return Chi2_table95[n-1];
	else if(p==0.99)
		return Chi2_table99[n-1];
	else if(p==0.999)
		return Chi2_table999[n-1];

	errif(true,"limitChi2: unimplemented Chi2(" << n
	  << ") " << p*100 << "%");

	return 0.0;
}
