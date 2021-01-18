
/* CbinMMtable.h
 *
 * Copyright (C) 2014 Mahmoud Ghandi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __gkmsvmXC__CbinMMtable__
#define __gkmsvmXC__CbinMMtable__

#include <stdio.h>
#include "global.h"

class CbinMMtable
{
public:
    CbinMMtable();

    int **table;
    int *dat;
    int L, Dmax;
    int nrow;
    
    void deleteTable(); // deletes the table
    int createTable(int L, int Dmax); // creates a table with all l-mer with max Dmax 1s.
    
    ~CbinMMtable(void);
};

#endif /* defined(__gkmsvmXC__CbinMMtable__) */
