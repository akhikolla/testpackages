/***************************************************************************
                             SRC/mixmod/Utilities/Random.cpp  description
    copyright            : (C) MIXMOD Team - 2001-2016
    email                : contact@mixmod.org
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD
    
    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

    All informations available on : http://www.mixmod.org                                                                                               
***************************************************************************/
// Tiny Encryption Algorithm
//--------------------------

#include "mixmod/Utilities/Random.h" // prototypes
#include "mixmod/Utilities/Util.h"
#include <time.h>
#include <fstream>

namespace XEM {

// Algorithm implemented for portability 32-64 bits, ...

const double m = 4294967296.0;

const uint32_t d = 0x09E3779B9L;
const uint32_t k0 = 0x0C7D7A8B4L;
const uint32_t k1 = 0x09ABFB3B6L;
const uint32_t k2 = 0x073DC1683L;
const uint32_t k3 = 0x017B7BE43L;

uint32_t y = 123456789L;
uint32_t z = 987654321L;

// Return a value in [0...1[
double rnd() {
	uint32_t s = 0;
	uint32_t n = 8;
	while (n-- > 0) {
		s += d;
		y += (z << 4) + (k0^z) + (s^(z >> 5)) + k1;
		z += (y << 4) + (k2^y) + (s^(y >> 5)) + k3;
	}
	return (z + y / m) / m;
}

int64_t flip(double x) {
	return (rnd() <= x);
}

void initRandomize(int seed)
{
	seed < 0
		? randomize()
		: antiRandomize(seed);
}

// Randomly initialize seeds (non-deterministic runs)
void randomize() {
	struct timespec tp;

	clock_gettime(CLOCK_REALTIME, &tp);

	z = (uint32_t) tp.tv_nsec/1000000; //milliseconds
	y = (uint32_t) tp.tv_sec;
	rnd();
}

// Use user seed (deterministic runs)
void antiRandomize(int seed) {
	// bijection N-->N^2 by diagonals [TODO: check formulas]
	//1) look for the largest n such that n(n+1)/2 <= seed
	int n = (int)floor(0.5 * (-1.0 + sqrt(1.0 + 8.0 * seed)));
	//2) coordinates are seed-n(n+1)/2, (n+1)(n+2)/2-1-seed
	z = seed - (n*(n+1))/2;
	y = ((n+1)*(n+2))/2 - 1 - seed;
	rnd();
}

}
