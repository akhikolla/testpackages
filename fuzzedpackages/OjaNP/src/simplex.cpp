/* 
 *  simplex.cpp : Subroutines for d-dimensional simplex.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: simplex.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "simplex.h"
#include <math.h>

using namespace std; //df

//  CLASS: Simplex
//  ==============

Simplex::Simplex()
{
}

Simplex::~Simplex()
{
}

Simplex::Simplex(int d)
{
	M=matrix(d+1,d+1);
    for(int c=0; c<d+1; c++)
		M(0,c) = 1.0;
}

Simplex::Simplex(const Simplex& S)
{
	M=S.M;
}

Simplex& Simplex::operator=(const Simplex& S)
{
	M=S.M;
	return *this;
}

double Simplex::det() const
{
	return M.det();
}

double Simplex::size() const
{
    errif(dim()==0,"Simplex::size: no data");
    
    return fabs(det()) / fact(dim());
}

double Simplex::sign() const
{
    errif(dim()==0,"Simplex::sign: no data");
    
    double z = det();
    if(z < 0)
		return -1.0;
    else if (z > 0)
		return 1.0;
    else
		return 0.0;
}

void Simplex::get(const Data& D,const Index& I,const Point& V)
{
    if(dim() != D.dim())
	M = matrix(D.dim()+1,D.dim()+1);
    
    errif(V.dim() != D.dim(),"Simplex::get: Point '"
	  << V << "' has wrong dimension");
    errif(I.limit() > D.size() ||
	  I.dim() != D.dim(),"Simplex::get: index '"
	  << I << "' is not valid"); 
    
    for(int c=0; c<dim()+1; c++)
	M(0,c) = 1.0;
    
    for(int r=0; r<dim(); r++)
	M(r+1,dim()) = V.coord(r);
    
    for(int c=0; c<dim(); c++)
	for(int r=0; r<dim(); r++)
	    M(r+1,c) = D[I[c]][r];
}

void Simplex::get(const Data& D,const Index& I)
{
    Point Zero(D.dim());   
    get(D,I,Zero);
}

void Simplex::get(const vector<Point> points)
{
	errif(points.size() == 0, "Simplex::get: vector points has no elements")
		int d = points[0].dim();

	if (dim() != d)
		M = matrix(d + 1, d + 1);

	errif(d + 1 != points.size(), "Simplex::get: vector points size '" << points.size()
		<< "' must equal dimension + 1");

	for (int c = 0; c < dim() + 1; c++)
		M(0, c) = 1.0;

	for (int c = 0; c < dim() + 1; c++)
	for (int r = 0; r < dim(); r++)
		M(r + 1, c) = points[c].coord(r);
}

void Simplex::set_column(int col,const Point& v)
{
	errif(col < 0 || col >= dim()+1,"Simplex::get_column: invalid column " << col);
	errif(dim() != v.dim(),"Simplex::get_column: dimension mismatch");
	
	for(int r=0; r<dim(); r++)
	    M(r+1,col) = v.coord(r);
}

bool Simplex::contains(const Point& p)
{
	errif(dim() != p.dim(),"Simplex::contains: dimension mismatch");
	
	valarray<double> y(dim()+1);
	y[0]=1.0;
	for(int i=0; i<dim(); i++)
		y[i+1]=p.coord(i);

	valarray<double> z=M.solve(y);

	for(int i=0; i<dim(); i++)
		if(z[i] < 0 || z[i] > 1)
			return 0;
	
	return 1;
}

ostream& operator <<(ostream& os,const Simplex& S)
{
    os << S.M;

    return os;
}
