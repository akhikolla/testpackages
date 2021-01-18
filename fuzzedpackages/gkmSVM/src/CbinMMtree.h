
/* CbinMMtree.h
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

#ifndef __gkmsvmXC__CbinMMtree__
#define __gkmsvmXC__CbinMMtree__

#include <stdio.h>
#include "global.h"

class CbinMMtree
{
public:
    CbinMMtree();
    CbinMMtree *child0;     //match
    CbinMMtree *child1;     //mismatch
    int addSeq(int *seq, int n); //adds a binary sequence of length n to the tree
    int deleteTree(); // deletes the tree

    int addtree(int n0, int n1); //makes a tree consisting all the possible combinations of n0 zeros and n1 ones. returns number of leaves
    int addLDtree(int L, int Dmax); // makes a tree consisting of all L-mers with at most Dmax 1s. returns number of leaves
    int addTreeToTable(int ** table, int frompos, int n, int *tmpArray); // copies all the n-mers to a table. frompos should be 0 if called from outside 
    
    /*
    int makeLDtable(int **table, int L, int Dmax); //fills in the elements in a table all the possible combinations of n0 zeros and n1 ones. returns number of rows added
    int makeTable(int **table, int n0, int n1);    //fills in the elements in a table all the possible L-mers with at most Dmax 1s. returns number of rows added
     */
    double calcAddCost(int *lmer, double *w,  int L, double p); // calculates the cost of the additional edges // p = 1/b (prob of match)

    ~CbinMMtree(void);
};




#endif /* defined(__gkmsvmXC__CbinMMtree__) */
