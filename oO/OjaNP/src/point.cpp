/* 
 *  point.cpp : Vector arithmetics.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: point.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "point.h"
#include "error.h"
#include "matrix_wrapper.h"
#include <math.h>
#include <list>

using namespace std;

//  CLASS: Point
//  ============
Point::Point(int n) : v(valarray<double>(n))
{
	for(int i=0; i<n; i++)
		v[i]=0.0;
}

Point& Point::operator=(const string& s0)
{
    list <double> pi;
    string s=s0;
	
    if(s=="nil" || s=="")
    {
		*this = Point();
		return *this;
    }
    
    for(;;)
    {
		int split;

		while(s[0]==' ')
			s=s.substr(1,s.length());
		
		if(s=="")
			break;
		
		pi.push_back(atof(s.c_str()));
		split=s.find(",");
		if(split < 0)
			split=s.find(" ");	
		if(split < 0)
			break;
		
		s=s.substr(split+1,s.length());
    }
    
    errif(pi.size()==0,"Point::operator=: string format '" << s0 << "' not valid");
    
    Point p(pi.size());

	list<double>::iterator i;
	int j=0;
	
	for(i=pi.begin(); i!=pi.end(); i++)
		p[j++]=*i;
	
    *this = p;
	
    return *this;
}

Point& Point::operator=(const Point& p)
{
	v.resize(p.dim());
	
	for(int i=0; i<dim(); i++)
		v[i]=p.coord(i);

	return *this;
}

double Point::length() const
{
	double sum=0.0;
	
	for(int i=0; i<dim(); i++)
		sum+=v[i]*v[i];

	return sqrt(sum);
}

bool operator==(const Point& x1,const Point& x2)
{
    errif(x1.dim() != x2.dim(),"Point::operator==: dimensions " << x1.dim() 
	  << " and " << x2.dim() << " do not match");
	
	for(int i=0; i<x1.dim(); i++)
		if(x1.coord(i) != x2.coord(i))
			return false;

	return true;
}

bool Point::in_box(const Point& min,const Point& max) const
{
    errif(min.dim() != max.dim(),"Point::in_box: dimensions " << min.dim() 
	  << " and " << max.dim() << " do not match");
	
    if(min.dim()==0)
		return false;
    errif(dim() != dim(),"Point::in_box: dimensions " << dim()
	  << " and " << max.dim() << " do not match");
	
    for(int i=0; i<dim(); i++)
		if(coord(i) < min.coord(i) || coord(i) > max.coord(i))
			return false;
    
    return true;
}

Point operator*(const matrix& M,const Point& x)
{
	Point y(x.dim());
	
	for(int i=0; i<x.dim(); i++)
		for(int j=0; j<x.dim(); j++)
			y[i]+=M.at(i,j) * x.coord(j);

	return y;
}

void Point::normalize()
{
    errif(is_nil(),"Point::normalize: point is nil");
    
    double n = length();
    // BUG: ep�tarkkaa vertailua
    errif(n==0.0,"Point::normalize: vector length is 0");
    if(n != 1.0)
		for(int i=0; i<dim(); i++)
			v[i] /= n;
}

matrix covariance(const Point& x,const Point& y)
{
	errif(x.dim() != y.dim(),"covariance: dimensions " << x.dim()
	  << " and " << y.dim() << " do not match");

	matrix M(x.dim(),y.dim());
	
	for(int i=0; i<x.dim(); i++)
		for(int j=0; j<y.dim(); j++)
			M(i,j)=x.coord(i)*y.coord(j);

	return M;
}

ostream& operator <<(ostream& os,const Point& p)
{
    if(p.is_nil())
    {
		return os;
    }
    
    for(int i=0; i<p.dim(); i++)
    {
		if(i)
			os << ' ';
		os << p.coord(i);
    }    
	
    return os;
}

double operator|(const Point& x1,const Point& x2)
{
    errif(x1.dim() != x2.dim(),"operator|: dimensions of ("
	  << x1 << ") and (" << x2 << ") do not match");
    double sum=0.0;
    for(int i=0; i<x1.dim(); i++)
		sum += x1.coord(i) * x2.coord(i);
    
    return sum;
}

double Point::angle(const Point& p2) const
{
    errif(dim() != p2.dim(),"operator|: dimensions of ("
	  << *this << ") and (" << p2 << ") do not match");
    
    double a=acos(fabs(*this | p2) / (length() * p2.length()));
	
    return a*180.0 / M_PI;
}

istream& operator >>(istream& is,Point& p)
{
	double t;
	list<double> x;

	while(is.peek() != '\n' && is.peek() != '\r')
	{
		is >> t;
		if(!is)
			break;
		
		x.push_back(t);
	}

	char c;
	while(is.peek() == '\n' || is.peek() == '\r')
		is.get(c);
	
	p=Point(x.size());
	list<double>::iterator i;
	int j=0;
	for(i=x.begin(); i!=x.end(); i++)
		p[j++]=*i;
	
    return is;
}


