
/* GTreeLeafData2.h
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


#ifndef __gkmsvmXC__GTreeLeafData2__
#define __gkmsvmXC__GTreeLeafData2__

#include <stdio.h>
#include "global.h"

class GTreeLeafData2  // handles addLTreeSnodeData instead of individual seqIDs , hence hopefully faster and less memory req.
{
public:
    
    int n; // number of L-mers in the list
    LTreeSnodeDataptr seqIDsets; //if n==1, it is LTreeSnodeData* and contains one seqIDset, otherwise it is LTreeSnodeData** and is the array of seqIDsets;  // each seqIDset is the list of seqIDs containing one l-mer
    intintptr gbits; // gbits for the case n==1, or array of gbits for n>1
    
    GTreeLeafData2(void);
    ~GTreeLeafData2(void);
    
    //void add(int seqID, int gbits);
    void addLTreeSnodeData(LTreeSnodeData *nodeData, int curGapPosSeq);
    
    
    void process();// calculates the mismatch profiles
    int calcdist(int difx); // calculates number of mismatches
private:
    
};

extern	int ***gMMProfile; //mismatchprofile[seqidi][mm][seqidj]


#endif /* defined(__gkmsvmXC__GTreeLeafData2__) */
