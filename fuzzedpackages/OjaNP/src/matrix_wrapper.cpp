/* 
 *  matrix_wrapper.cpp : Interface to matrix class.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: matrix_wrapper.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "error.h"
#include "matrix_wrapper.h"

using namespace std; //df

#if 0
double cof(const matrix& M,int r,int c)
{
    errif(M.rows() == 0,"cof: matrix M has no elements");
    errif(M.rows() != M.columns(),"cof: matrix " << M.rows() << "x" << 
	  M.columns() << " is not quadratic");
	
    int n=M.rows();
	
    errif(r < 0 || r >= n,"cof: illegal row index " << r);
    errif(c < 0 || c >= n,"cof: illegal column index " << c);
    
    matrix C(n-1,n-1);
    
    for(int i1=0,i=0; i < n-1; i++,i1++)
    {
		if(i1==r)
			i1++;
		
		for(int j1=0,j=0; j < n-1; j++,j1++)
		{
			if(j1==c)
				j1++;
			
			C(i,j) = M.at(i1,j1);
		}
    }
    
    return ((r+c) % 2 ? -1.0 : 1.0) * C.det();
}
#else
double cof(const matrix& M,int r,int c)
{
	return M.Cofact(r,c);
}
#endif

matrix cof(const matrix& M)
{
    errif(M.rows() != M.columns(),"cof: matrix " << M.rows() << "x" << 
	  M.columns() << " is not quadratic");
	
    int n=M.rows();
    matrix C(n,n);
    
    for(int i=0; i<n;i++)
		for(int j=0; j<n;j++)
			C(i,j) = cof(M,i,j);
	
    return C;
}

valarray<double> to_vector(const matrix& M)
{
	errif(M.columns()!=1,"to_vector(const matrix&): invalid matrix " << M.rows() << "x" << M.columns());

	valarray<double> ret(M.rows());
	for(int i=0; i<M.rows(); i++)
		ret[i]=M.at(i,0);

	return ret;
}

matrix to_matrix(const valarray<double>& v)
{
	matrix M(v.size(),1);

	for(int i=0; i<int(v.size()); i++)
		M(i,0)=v[i];

	return M;
}

