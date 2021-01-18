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
 * File:   memory.h
 * Author: Mohammad H. Ferdosi
 *
 * Created on 27 October 2012, 3:25 PM
 */

#ifndef MEMORY_H
#define MEMORY_H

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "TypeConversion.h"

class memoryCLS
{
public:
    int *pMemory;
    int *pPreviousMemory;
    int *pValidity;

    memoryCLS(int* nrow = 0);
    memoryCLS(const memoryCLS& orig);

    int memoryMaker(int * col);
    int memoryMaker(int* col, int *validation);
    int evalution(int *swichLocations);
    int freeMemory();
    virtual ~memoryCLS();

private:
    int *pRows;
};

#endif	/* MEMORY_H */

