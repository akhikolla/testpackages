// Copyright (C) 2014 Mohammad H. Ferdosi
//
// HSPhase is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// HSPhase program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/*
 * diag.h
 *
 *  Created on: 17/01/2013
 *      Author: Mohammad H. Ferdosi
 */

#ifndef DIAG_H_
#define DIAG_H_
#include "memory.h"
#include "hsp.h"
void diagnostic(int const * matrix, int *nrow, int *ncol, int * result, vector<int> &hetSite,int* zeroFrq,  int* twoFrq);

#endif /* DIAG_H_ */
