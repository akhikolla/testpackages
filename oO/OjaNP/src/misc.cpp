/* 
 *  misc.cpp : Miscellanous subroutines.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: misc.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <limits.h>
#include <fstream>
#include "error.h"
#include "misc.h"


using namespace std;

//  MISC: fact, choices
//  ===================

int fact(int k)
{
   // static int factorials_max = 32;
    static int factorial[32];
    static int next_factorial=0;
	
    errif(k < 0,"fact: illegal parameter " << k);
  //  errif(k >= factorials_max,"fact: overflow");
	
    if(k < next_factorial)
		return factorial[0];  //Here was factorial[k];

    for(int i=next_factorial; i<=k; i++)
    {
		if(i==0)
			factorial[i]=1;
		else
		{
			errif(factorial[i-1] > INT_MAX/i,"fact: overflow");
			factorial[i]=factorial[i-1]*i;
		}
    }
	
    return factorial[k];
}

unsigned long choices(int n, int k)
{
	errif(n < 1,"choices: n is less than 1 (n = " << n << ")");
	errif(k < 1 || k > n,"choices: illegal value k (k = " << k 
 	  << ", n = " << n << ")");

	unsigned long r = n--; unsigned int d = 2;
	while (d <= k){ r *= n--; r /= d++; }
	return r;
}

double round(double x,double prec)
{
	return double(int((x/prec)+0.5))*prec;
}

bool is_file(const char* s)
{
	ifstream I(s);
	if(I)
		return true;

	return false;
}

void get_partitions(list<vector<int> >& P,int N)
{
	int n;
	vector<int> I,J;

	I=vector<int>(1);
	I[0]=N;

	P.clear();
	if(N < 1)
		return;

	P.push_back(I);
	if(N==1)
		return;
	
	for(;;)
	{
		I=P.back();
		n=I[I.size()-1];

		if(n != 1)
		{
			J=vector<int>(I.size()+1);
			for(size_t j=0; j<I.size()-1; j++)
				J[j]=I[j];
			J[I.size()-1]=n-1;
			J[I.size()]=1;
		}
		else
		{
			size_t j=0; // Ykk�st� edelt�v� alkio
			while(j<I.size() && I[j]!=1)
				j++;
			j--;
			
			int sum=0; // Alkioiden t�m�n hetkinen summa.
			
			J=vector<int>(j);
			
			for(size_t k=0; k<j; k++)
			{
				J[k]=I[k];
				sum+=I[k];
			}

			while(sum+I[j]-1 <= N)
			{
				J.push_back(I[j]-1);
				sum+=I[j]-1;
			}

			if(sum < N)
				J.push_back(N-sum);
		}

		P.push_back(J);
		if(J.size()==size_t(N))
			break;
	}
}

