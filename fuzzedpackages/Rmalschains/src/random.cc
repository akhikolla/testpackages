/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "random.h"
#include <iostream>
#include <math.h>
#include <cassert>

using namespace std;

Random::Random(IRealRandom *random) {
	assert(random != NULL);
	this->random = random;
}

Random::~Random(void) {
	delete this->random;
}


const double PI=3.1415926;

double Random::normal(double desv) {
  double u1, u2;
  double result;
  
  do {
    u1=Random::rand();
  } while (u1 == 0);

  u2=Random::rand();
  result = desv * sqrt (-2*log(u1)) * sin (2*PI*u2);
  
  return result;
}

void initSample(int *sample, int max) {
   int i;

    for (i = 0; i < max; i++) {
	sample[i] = i;
    }
}

int Random::getSample(int *sample, int *pmax) {
    int max = *pmax;
    int r, pos;

    assert(max >= 0);
    r = randint(0, max-1);
    pos = sample[r];
    sample[r] = sample[max-1];
    --max;

    *pmax = max;
    return pos;
}

unsigned long Random::getSeed(void) {
    return random->getSeed();
}

