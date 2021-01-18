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
 * File:   swDetect.cpp
 * Author: Mohammad H. Ferdosi
 *
 * Created on 27 October 2012, 3:33 PM
 */

#include "swDetect.h"

using namespace std;

swDetect::swDetect(int *R /* = 0 */,  int *nrow /* = 0 */,  int *ncol /* = 0 */)
{

    col = (*ncol);
    row = (*nrow);
    pRows = new int*[*nrow];

    for (int i = 0; i < *nrow; i++)
    {
        pRows[i] = R + (*ncol)*i;
    }
}


swDetect::~swDetect()
{
    delete pRows;
}



int swFun( int *matrix1,  int *matrix2,  int *nrow,  int *ncol, int *result)
{
    swDetect mat1(matrix1, nrow, ncol);
    swDetect mat2(matrix2, nrow, ncol);

    int mySwitch = 0;

    for (int j = 0; j<*nrow; j++)
    {
        for (int i = 0; i<*ncol; i++)
        {
            if (mat1.pRows[j][i] != mat2.pRows[j][i] && mat1.pRows[j][i] != 0 && mat2.pRows[j][i]!=0)
            {
                mySwitch = mySwitch + 1;
                for (int k = i; k < *ncol; k++)
                {
                    if (mat2.pRows[j][k] == 3)
                        mat2.pRows[j][k] = 4;
                    else if (mat2.pRows[j][k] == 4)
                        mat2.pRows[j][k] = 3;
                }
            }
        }

    }
    *result = mySwitch;
    return (0);
}


