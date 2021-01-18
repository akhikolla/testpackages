
/* CiDLPasses.h
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

#ifndef __gkmsvmXC__CiDLPasses__
#define __gkmsvmXC__CiDLPasses__

#include <stdio.h>
#include "CbinMMtree.h"

class CiDLPasses
{
public:
    CiDLPasses();
    int L; // Length
    int M; // number of passes
    //int Dmax; // maximum number of mismatches
    int **passOrder; // order of the bases to be traversed. each of passOrder[0..M-1] is a permutation of 0..L-1
    CbinMMtree **passTrees; // for each of the M passes, there is a tree that holds all the mismatch (binary) L-mers that are checked in that pass. each L-mer should belong to one pass (tree) . the ordering is based on pass[i][] order
    
    void initPassOrderIDL(int L); // generates M=L circular pass orders
    void initPassOrderIDL(int L, int M, int Dmax); // generates M=L or 2L circular pass orders
    void initPassOrderAll(int L, int Dmax); // generates choose(L, Dmax) permutations
    void newPassOrderDesignCover(int L, int Dmax, int k);// generates M passes, that gaurantee that the first k places are matches (in all the trees)
        
    void newIDLPasses(int L, int Dmax); // generates M circular passes and assigns L-mers to the trees based on min iDL rule.
    void newGreedyIDLPasses(int L, int M, int Dmax, int *nodesAtDepthCnt, double p); // uses the greedty alg to assign to trees. p =1/b prob(match)
    void newGreedy2IDLPasses(int L, int M, int Dmax, int *nodesAtDepthCnt, double p); // uses the greedty 2alg to assign to trees. p =1/b prob(match)

    void deletePassOrder();
    double calcSlope(int *lmer, int *order, int L); // calculates slope (max #1s / total)
    double calcCost(int *lmer, int *order, double *w, double p, int L); // calcs a cost for traversing a lmer in specific order. for Greedy2alg
    
    int *reorder(int *lmer, int *order, int L, int *output); // reorders the L-mer
    
    int isCoprime(int a, int b); // checks if a and b are co-prime
    
    ~CiDLPasses(void);
};


#endif /* defined(__gkmsvmXC__CiDLPasses__) */
