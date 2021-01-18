/* 
 *  line.cpp : Geometry of line.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: line.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "line.h"
#include "data.h"
#include "error.h"
#include "matrix_wrapper.h"
#include <math.h>
#include <valarray>

using namespace std; //df

//  CLASS: Line
//  ===========

Line::Line(const Point& p0,const Point& d)
{
    errif(p0.dim() != d.dim(),"Line::Line: dimensions do not match");
    start = p0;
    direction = d;
}

void Line::set_dim(int d)
{
	if(d != dim())
	{
		start=Point(d);
		direction=Point(d);
	}
}

Point Line::at(double t) const
{
    Point p(start + t * direction);
	
    return p;
}

double Line::x0(int i) const
{
    errif(i < 0 || i >= dim(),"Line::x0: illegal index " << i);
    return start.coord(i);
}

double Line::dx(int i) const
{
    errif(i < 0 || i >= dim(),"Line::dx: illegal index " << i);
    return direction.coord(i);
}

bool Line::is_nil() const
{
    for(int i=0; i<dim(); i++)
		if(direction.coord(i))
			return false;
	
    return true;
}

void Line::get(HyperplaneSet& H)
{
    errif(H.size() != H.dim()-1,"Line::get: incompatible size and dimension ("
	  << H.size() << " and " << H.dim() << ")");
	
    int n=H.dim();
    
	set_dim(n);
    
    n--;
    valarray<double> d0(n),d1(n);
	Point h0(n),h1(n);
    matrix D2(n,n);
	
    // K�yd��n l�pi kaikki dimensiot 'k' ja ratkaistaan muut muuttujat
    // muuttujan 'k' suhteen.
    for(int k=0; k<n+1; k++)
    {
		for(int i=0; i<n; i++)
		{
			d0[i] = H[i].c();
			d1[i] = H[i][k];
			for(int j=0; j<n; j++)
				D2(i,j) = H[i][j+(j >= k)];
		}
		
		if(D2.det()==0)
			continue;
		
		D2 = D2.inv();
		h0 = -(D2*Point(d0));
		h1 = -(D2*Point(d1));
		
		for(int i=0; i<n+1; i++)
		{
			if(i==k)
			{
				start[i] = 0;
				direction[i] = 1;
			}
			else
			{
				start[i] = h0[i-(i >= k)];
				direction[i] = h1[i-(i >= k)];
			}
		}
		
		return;
    }
    
    Point zero(dim());
    
    start=direction=zero;
}

void Line::get(const Data& D,const IndexSet& I)
{
    errif(D.size()==0,"OjaLine::get:no data");
    errif(I.indices()+1 != D.dim(),"OjaLine::get: wrong number of indices " 
	  << I.indices());
    errif(I.dim() != D.dim() ,"OjaLine::get: illegal index dimension " 
	  << I.dim());

    HyperplaneSet H;
    H.get(D,I);
    get(H);
}

double Line::angle(const Point& z) const
{
    errif(dim() != z.dim(),"Line::angle: dimensions do not match");
	
    return direction.angle(z);
}

Point Line::proj(const Point& z) const
{
    errif(dim() != z.dim(),"Line::proj: dimensions do not match");
    
    double d=(direction | z) / (direction.length() * direction.length());
    
    return d * direction;
}

void Line::get_through(const Point& p1,const Point& p2)
{
    if(p1==p2)
    {	
		Point zero(dim());
        start=direction=zero;
		return;
    }
    
    start=p1;
    direction=(p2-p1);
}

void Line::get_toward(const Point& p0,const Point& dir)
{
	start=p0;
	direction=dir;
}

ostream& operator <<(ostream& os,const Line& L)
{
    if(L.is_nil())
    {
		os << "(nil)";
		return os;
    }
    
    os << '(';   
    for(int i=0; i<L.dim(); i++)
    {
		if(i)
			os << ',';
		if(L.start.coord(i))
			os << L.start.coord(i);
		if(L.direction.coord(i) < 0)
			os << '-';
		else if(L.direction.coord(i) > 0)
		{
			if(L.start.coord(i))
				os << '+';
		}
		else
		{
			if(!L.start.coord(i))
				os << '0';
			
			continue;
		}
		
		if(fabs(L.direction.coord(i)) != 1.0)
			os << fabs(L.direction.coord(i));
		os << 't';
    }
    os << ')';
	
    return os;
}

