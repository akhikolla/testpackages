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

#include "orandom.h"

#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9

void ORandom::setSeed(unsigned long seed) {
   m_seed = seed;
}

double  ORandom::rand(void) {
   m_seed = (m_seed * PRIME) & MASK;
   return (m_seed * SCALE);
}

unsigned long  ORandom::getSeed(void) {
   return m_seed;
}
