
/* CbinMMtable.cpp
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

#include "CbinMMtable.h"
#include "CbinMMtree.h"



CbinMMtable::CbinMMtable(void){
    table=nullptr;
    dat=nullptr;
    L=Dmax=nrow=0;
}

int  CbinMMtable::createTable(int L, int Dmax){
    CbinMMtree *mmtree = new CbinMMtree();
    nrow = mmtree->addLDtree(L,Dmax);
    this->L = L;
    this->Dmax = Dmax;
    dat = new int[nrow*L];
    table = new int*[nrow];
    for(int i=0;i<nrow;i++){
        table[i]= dat+i*L;
    }
    int *tmpArray = new int[L];
    mmtree->addTreeToTable(table, 0, L, tmpArray);
    delete []tmpArray;
    mmtree->deleteTree();
    delete mmtree;
    return nrow;
}
    
void CbinMMtable::deleteTable(){
    if(dat!=nullptr){
        delete []dat;
        delete []table;
        table=nullptr;
        dat=nullptr;
        L=Dmax=nrow=0;
    }
}

CbinMMtable::~CbinMMtable(void){
    this->deleteTable(); 
}

/*
 class CbinMMtable tt;
 tt.createTable(4, 2);
 for(int i=0;i<tt.nrow;i++){
   printf("\n");
   for(int j=0;j<4;j++){
     printf("%d", tt.table[i][j]);
   }
 } 
 */
