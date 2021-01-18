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
 * File:   swDetect.h
 * Author: Mohammad H. Ferdosi
 *
 * Created on 27 October 2012, 3:33 PM
 */

#ifndef SWDETECT_H
#define	SWDETECT_H

#include <iostream>

class swDetect
{
public:

    int **pRows;
    int col;
    int row;

    swDetect(int* R = 0, int* nrow = 0,  int* ncol = 0);


    virtual ~swDetect();
private:

};
int swFun( int *matrix1,  int *matrix2,  int *nrow,  int *ncol,  int *result);


#endif	/* SWDETECT_H */

