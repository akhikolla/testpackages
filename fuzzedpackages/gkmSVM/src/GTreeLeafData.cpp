
/* GTreeLeafData.cpp
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

#include "GTreeLeafData.h"

GTreeLeafData::GTreeLeafData(void){
    this->n=0;
}

GTreeLeafData::~GTreeLeafData(void){
}

void GTreeLeafData::add(int seqID, int gbits){
 
    if(n==0){
        n = 1;
        this->seqIDs_gbits.i = seqID;
        this->first_gbits = gbits;
    }else{
        
        if (n==1){
            intintptr newseqids;
            newseqids.p= new int[2*2];
            newseqids.p[0]= seqIDs_gbits.i;
            newseqids.p[1]= first_gbits;

            newseqids.p[2]= seqID;
            newseqids.p[3]= gbits;
            
            this->seqIDs_gbits = newseqids;
            n = 2;
        }else{
            // n>1
            
            
            if ((n & (n-1))==0) // i.e. n is power of 2
            {
                // expand memory
                intintptr newseqids;
                newseqids.p= new int[n<<2];
                for(int j=0;j< (2*n); j++)
                {
                    newseqids.p[j]= seqIDs_gbits.p[j];
                }
                delete []seqIDs_gbits.p;
                seqIDs_gbits.p = newseqids.p;
            }
            seqIDs_gbits.p[2*n] = seqID;
            seqIDs_gbits.p[2*n+1] = gbits;
            n++;
        }
    }
    
}
 

void GTreeLeafData::addLTreeSnodeData(LTreeSnodeData *nodeData, int curGapPosSeq){
    // adds many sequences all at once
   /*
    if(n==0){
        n = nodeData->n;
        this->seqIDs_gbits.i = seqID;
        this->first_gbits = gbits;
    }else{
        
        if (n==1){
            intintptr newseqids;
            newseqids.p= new int[2*2];
            newseqids.p[0]= seqIDs_gbits.i;
            newseqids.p[1]= first_gbits;
            
            newseqids.p[2]= seqID;
            newseqids.p[3]= gbits;
            
            this->seqIDs_gbits = newseqids;
            n = 2;
        }else{
            // n>1
            
            
            if ((n & (n-1))==0) // i.e. n is power of 2
            {
                // expand memory
                intintptr newseqids;
                newseqids.p= new int[n<<2];
                for(int j=0;j< (2*n); j++)
                {
                    newseqids.p[j]= seqIDs_gbits.p[j];
                }
                delete []seqIDs_gbits.p;
                seqIDs_gbits.p = newseqids.p;
            }
            seqIDs_gbits.p[2*n] = seqID;
            seqIDs_gbits.p[2*n+1] = gbits;
            n++;
        }
    }
    
    */
}

void GTreeLeafData::process(){// calculates the mismatch profiles
    if(n==1){
        //    gMMProfile[seqIDs_gbits.i][0][seqIDs_gbits.i]++;  // it is not needed as it can be computed based on the sequence length !

    }else{
        
        int *iSidGb=seqIDs_gbits.p;
        for(int i=0;i<n;i++){
            int seqIDi =*iSidGb++;
            int gbitsi =*iSidGb++;
            
            int *jSidGb=seqIDs_gbits.p;
            
            int **gMMProfile_seqIDi = gMMProfile[seqIDi];
            //gMMProfile_seqIDi[0][seqIDi]++;//itself  // it is not needed as it can be computed based on the sequence length !
            for(int j=0;j<i;j++){
//              for(int j=0;j<n;j++){

                int seqIDj =*jSidGb++;
                int gbitsj =*jSidGb++;
                
                int dif=calcdist(gbitsj^gbitsi);
                gMMProfile_seqIDi[dif][seqIDj]++; //itself

            }
        }
        
    }

}

int GTreeLeafData::calcdist(int difx){
    int res = 0;
    while(difx>0){
        res+=((difx&NBITSONES)!=0);
        difx>>=NBITS;
    }
    return(res);
 }
