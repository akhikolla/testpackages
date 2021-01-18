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
 * File:   block4Phase.h
 * Author: Mohammad H. Ferdosi
 *
 * Created on 27 October 2012, 3:34 PM
 */
#ifndef BLOCK4PHASE_H
#define	BLOCK4PHASE_H

#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>

#define uint unsigned int

using namespace std;


class SNP
{
public:

    vector<uint> strand1;
    vector<uint> strand2;


    SNP();
    virtual ~SNP();
    SNP recombination(unsigned int index);

};

class block4Phase
{
public:

    block4Phase(const uint *matrix, const uint *nrow, const uint *ncol, uint *result, const uint *siregenotype, const uint * str);
    block4Phase(const block4Phase& orig);
    virtual ~block4Phase();

    int sireStrdDetector(const SNP &sire,const SNP &halfsib);
    int blockMaker(SNP &sire, const SNP &halfsib, int * block = 0, const uint* ncol = 0);

private:
    int str_ ;

};


#endif	/* BLOCK4PHASE_H */

