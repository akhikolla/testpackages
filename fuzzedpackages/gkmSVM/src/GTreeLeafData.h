
/* GTreeLeafData.h
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

#ifndef __gkmsvmXC__GTreeLeafData__
#define __gkmsvmXC__GTreeLeafData__

#include <stdio.h>
#include "global.h"

class GTreeLeafData
{
public:

    int n; // number of L-mers in the list
    intintptr seqIDs_gbits; //if n==1, it is int and contains the ID, otherwise it is int* and is the array of IDs;
    //  LPTr Lmers; //pointer to the starts of the sequences
    int first_gbits; // gbits for the case n==1, otherwise, seqIDs and gapped_bits are both written in seqIDs_gbits (2 numbers for each L-mer)

    GTreeLeafData(void);
    ~GTreeLeafData(void);
    
    void add(int seqID, int gbits);
    void addLTreeSnodeData(LTreeSnodeData *nodeData, int curGapPosSeq);
    
    
    void process();// calculates the mismatch profiles
    int calcdist(int difx); // calculates number of mismatches
private:
    
};

extern	int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]


#endif /* defined(__gkmsvmXC__GTreeLeafData__) */
