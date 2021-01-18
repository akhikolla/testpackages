
#include "doe_random.h"
#include <time.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static long seed=-19752;

void srandomt()
{
	if(seed>0) seed=-seed;
	srandom( -(long)time( NULL ) +seed);
}

void srandom(long newseed)
{
	if(newseed>0) seed=-newseed;
	else seed=newseed;
}

double Random()
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (seed <= 0 || !iy) {
		if (-(seed) < 1) seed=1;
		else seed = -(seed);
		for (j=NTAB+7;j>=0;j--) {
			k=(seed)/IQ;
			seed=IA*(seed-k*IQ)-IR*k;
			if (seed < 0) seed += IM;
			if (j < NTAB) iv[j] = seed;
		}
		iy=iv[0];
	}
	k=(seed)/IQ;
	seed=IA*(seed-k*IQ)-IR*k;
	if (seed < 0) seed += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = seed;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
