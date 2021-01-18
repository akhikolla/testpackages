#ifndef __paircompare_h__
#define __paircompare_h__

#include <vector>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

vector<long long> pairCompare(const NumericMatrix& X,
			const int curMarker, 
			const long start1, 
			const long end1, 
			const long start2, 
			const long end2, 
			const vector<long>& index);

/*helper function to conduct the comprision count for Mann-Whitney count*/
/*
	Variable definition: 
			int curMarker:  current working marker
			long  start1 :  the starting position of the genotype patient ID index of cu-
							rent type.
			long  end2   :  the ending pistion of genotype patient ID index of curent ty-
							pe.
			long  start2 :  the starting position of the genotype patient ID index of co-
							mparing range.
			long  end2   :  the ending pistion of genotype patient ID index of curent co-
							mparing range.

*/
vector<long long> pairCompare(const NumericMatrix& X, const int curMarker, 
			const long start1, 
			const long end1, 
			const long start2, 
			const long end2, 
			const vector<long>& index)
{
    vector<long long> res(2,0);
    long long nonEqCount = 0; 
    long long eqCount    = 0;

    /* receive data from matrix X in a new vector and sort them*/    
    vector<double> range1(end1-start1+1,0);
    vector<double> range2(end2-start2+1,0);
    for(unsigned long i = 0; i < range1.size(); i++)
        range1[i] = X(index[start1+i],curMarker);
    for(unsigned long i = 0; i < range2.size(); i++)
        range2[i] = X(index[start2+i],curMarker);
	
    sort(range1.begin(),range1.end());
    sort(range2.begin(),range2.end());
 
    /* running index for comparison*/   
    long p = 0, q = 0;
    long tmpP =0, tmpQ=0;

    long range1Size = range1.size();
    long range2Size = range2.size();
 
/* 
    searching for a number ties in two sorted arrays using two pointer moving 
	is nonmatch is found, and scan for repeated case when match is found.
*/
    for(;p < range1Size && q < range2Size;)
    {
        if(range1[p] < range2[q])
        {   p++;}  // array will smaller element value will move one step
        else if(range1[p] > range2[q])
        {   q++;}
        else 
        {   
            //found tie, then checking repeating elements
            tmpP = p;
            tmpQ = q;
            while(tmpP < range1Size - 1)
            {
                if(range1[tmpP] == range1[tmpP+1])
                    tmpP++;
                else
                    break;
            }
            while(tmpQ < range2Size - 1)
            {
                if(range2[tmpQ] == range2[tmpQ+1])
                    tmpQ++;
                else 
                    break;
            }
            eqCount += (tmpP-p+1)*(tmpQ-q+1);
            p = tmpP+1;
            q = tmpQ+1;
        }
    }

    /* pair comparison count for unmatchness (<) for two sorted array*/
    long m = 0, n = 0;
    long N = range2.size()-1; //store the number of array to an int

    while( m < range1Size && n < range2Size)
    {
        if(range1[m] < range2[n])
        {   
            m++;
            nonEqCount += (N-n+1);
        }else
        {
            n++;
        }
    }
    
    res[0] = nonEqCount;
    res[1] = eqCount;
    
    return res;
}

#endif

