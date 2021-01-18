/* 
 *  hyperplane.cpp : Hyperplane geometry.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: hyperplane.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "hyperplane.h"
#include "error.h"
#include "misc.h"
#include "data.h"
#include "simplex.h"
#include "line.h"
#include <math.h>
#include <valarray>

using namespace std;


//  CLASS: Hyperplane
//  =================

Hyperplane::Hyperplane()
{
    cofs = 0;
    cof = 0;
    isBound = false;
}

Hyperplane::~Hyperplane()
{
    if(cof)
		delete[] cof;
	cof = 0;
}

Hyperplane::Hyperplane(int dim)
{
    errif(dim < 1,"Hyperplane::Hyperplane: non-positive dimension " << dim);
    cofs = dim+1;
    cof = new double[cofs];
    isBound = false;
    
    errif(!cof,"Hyperplane::Hyperplane: out of memory");
}

Hyperplane::Hyperplane(const Hyperplane& H)
{
    cofs = 0;
    cof = 0;
    operator=(H);
}

Hyperplane& Hyperplane::operator=(const Hyperplane& H)
{
	if (H.cof == 0) return *this;	// uninitialized hyperplane
    set_dim(H.dim());
    for(int i=0; i<cofs; i++)
		cof[i] = H.cof[i];
    isBound = H.isBound;

    return *this;
}

void Hyperplane::set_dim(int n)
{
    if(cofs != n+1)
	{
		if(cof)
			delete[] cof;
		cofs = n+1;
		cof = new double[cofs];
		errif(!cof,"Hyperplane::set_dim: out of memory");
	}
}

double Hyperplane::operator[](int index) const 
{
    errif(index < 0 || index >= dim(),"Hyperplane::operator[]: illegal index "
	  << index);
    return cof[index+1];
}

double Hyperplane::c() const
{
    errif(cofs==0,"Hyperplane::c: no data");
    return cof[0];
}

Point Hyperplane::normal() const
{
    Point c(dim());
    
    for(int i=0; i<dim(); i++)
		c[i] = cof[i+1];
	
    return c;
}

Point Hyperplane::cof_at(const Point& x) const
{
	static double k = (1.0 / double(fact(dim())));

	if(side(x) < 0)
		return  -k * normal();
	else
		return k * normal();
}

double Hyperplane::cof0_at(const Point& x) const
{
	static double k = (1.0 / double(fact(dim())));
	if(side(x) < 0)
		return -k * c();
	else
		return k * c();
	  
}

double Hyperplane::operator|(const Point& x) const
{
    errif(dim() != x.dim(),"Hyperplabe::operator|: dimensions " << dim()
	  << " and " << x.dim() << " do not match");
    
    return (x | normal());
}

void Hyperplane::get(const Data& D,const Index& I)
{
    errif(D.dim() != I.dim(),"Hyperplane::get: dimensions " << D.dim() 
	  << " and " << I.dim() << " do not match");
	
    set_dim(D.dim());
	Simplex S(D.dim());
    S.get(D,I);
    for(int i=0; i<cofs; i++)
		cof[i] = S.row_cof(i);
}

void Hyperplane::get(vector<Point> points)
{
	errif(points.size() == 0, "Hyperplane::get: vector points has no elements")
	// This would be checked in Simplex.get(points) // errif(D.dim() != I.dim(), "Hyperplane::get: dimensions " << D.dim() << " and " << I.dim() << " do not match");

	set_dim(points[0].dim());
	points.push_back(Point(dim())/*Zero*/);
	Simplex S;
	S.get(points);
	for (int i = 0; i<cofs; i++)
		cof[i] = S.row_cof(i);
}

bool Hyperplane::degenerated() const
{
    for(int i=1; i<cofs; i++)
		if(cof[i] != 0.0)
			return false;
	
    return true;
}

bool Hyperplane::intersect(const Line& L,double& t) const
{
    errif(L.dim() != dim(),"Hyperplane::intersect: dimensions do not match");
    
    double a,b;
    
    a = cof[0];
    b = 0.0;
    for(int i=0; i<dim(); i++)
    {
		a += cof[i+1] * L.x0(i);
		b += cof[i+1] * L.dx(i);
    }
	
    t = 0.0;
	if (abs(b)<1.0e-10)
		return false;
    
    t = -a/b;
    return true;
}

Point Hyperplane::intersect(const Line& L) const
{
    errif(L.dim() != dim(),"Hyperplane::intersect: dimensions do not match");
	
    Point zero; 
    double t;      
    
    if(intersect(L,t))
		return L.at(t);
    else
		return zero;
}

double Hyperplane::side(const Point& x) const
{
	errif(x.dim() != dim(), "Hyperplane::side: dimensions " << dim()
		<< " and " << x.dim() << " do not match");
	double sg;
	sg = ((*this) | x) + c();
	if (sg < 0)
		return -1.0;
	else if (sg > 0)
		return 1.0;
	else
		return 0.0;
}

void Hyperplane::normalize()
{
    if(!degenerated())
    {
		double length;

		length=normal().length();
		for(int i=0; i < cofs; i++)
			cof[i] /= length;
    }
}

double Hyperplane::dist(const Point& x) const
{
	errif(dim() != x.dim(),"Hyperplane::dist: dimensions " << dim()
	  << " and " << x.dim() << " do not match");

	return fabs((*this | x) + c()) / normal().length();
}

istream& operator >>(istream& is,Hyperplane& H)
{
	err("operator >>(istream& is,Hyperplane& H): not implemented");
	return is;
}

ostream& operator <<(ostream& os,const Hyperplane& H)
{
    if(!H.cof)
		return os;
	
    if(H.degenerated())
    { 
		os << "degenerated";
		return os;
    }
	
    bool first=true;
	
    for(int i=1; i < H.cofs; i++)
    {
		double c=H.cof[i];
		if(c < 0)
			os << "- ";
		else if (c > 0 && !first)
			os << "+ ";
		if(c != 0)
		{
			if(c != -1.0 && c != 1.0)
				os << fabs(c) << " ";
			os << "x" << i << " ";
			first = false;
		}
    }
    
    os << "= " << -H.cof[0];
	
    return os;
}

//  CLASS: HyperplaneSet
//  ====================

HyperplaneSet::~HyperplaneSet()
{
    if(plane)
		delete[] plane;
}

void HyperplaneSet::resize(int n)
{
    if(planes != n)
    {
		if(plane)
			delete[] plane;
		
		plane = 0;
		planes = n;
		if(n)
		{
			plane = new Hyperplane[n];
			errif(!plane,"HyperplaneSet::resize: out of memory");
		}
    }
}

HyperplaneSet::HyperplaneSet(int n)
{
    errif(n < 1,"HyperplaneSet::HyperplaneSet: non-positive size " << n);
    planes = n;
    plane = new Hyperplane[planes];
    errif(!plane,"HyperplaneSet::HyperplaneSet: out of memory");
}

HyperplaneSet::HyperplaneSet(const HyperplaneSet& HS)
{
    planes = 0;
    plane = 0;
    operator=(HS);
}

HyperplaneSet& HyperplaneSet::operator=(const HyperplaneSet& HS)
{
    resize(HS.planes);
    
    for(int i=0; i<planes; i++)
		plane[i] = HS.plane[i];
	
    return *this;
}

Hyperplane& HyperplaneSet::operator[](int i) const
{
    errif(i<0 || i>=planes,"HyperplaneSet::operator[]:Illegal index " << i);
	
    return plane[i];
}

void HyperplaneSet::get(const Data& D,const IndexSet& I)
{
    resize(I.indices());
	
    for(int i=0; i<I.indices(); i++)
		plane[i].get(D,I[i]);
}

void HyperplaneSet::get_all(const Data& D)
{
    Index I(D.dim(),D.size());
	int combnum = I.combinations();
	resize(combnum + MAX_BOUNDS);
	planes = combnum;
	
    for(int i=0; I; i++,I++)
		plane[i].get(D,I);
}

void HyperplaneSet::add(const Hyperplane& H)
{
	errif (bounds >= MAX_BOUNDS, "HyperplaneSet::add: too much bounds")

	plane[planes++] = H;
	bounds++;
}

void HyperplaneSet::get(const vector<Hyperplane>& hyperplanes){
	resize(hyperplanes.size());
	for (int i = 0; i<hyperplanes.size(); i++)
		plane[i] = hyperplanes[i];
}

Point HyperplaneSet::crossing_point() const
{
    errif(!size(),"HyperplaneSet::crossing_point: no data");
    errif(dim() != size(),"HyperplaneSet::crossing_point: " << size() << 
	  " planes but dimension was " << dim());
    
    int n=dim();
    matrix A(n,n);
    valarray<double> x(n),c(n);
	
    for(int i=0; i<n; i++)
    {
		c[i] = -plane[i].c();
		for(int j=0; j<n; j++)
			A(i,j) = plane[i][j];
    }

    if(A.det())
		x = A.solve(c);
	else
		return Point();
	
    return Point(x);
}

Point HyperplaneSet::crossing_point(const Index& I) const
{
    for(int i=0; i<I.dim(); i++)
		errif(I[i] >= planes,"HyperplaneSet::crossing_point: index overflow I["
	      << i << "] = " << I[i]);
	
    int n=dim();
    errif(n != I.dim(),"HyperplaneSet::crossing_point: dimensions "
	  << n << " and " << I.dim() << " do not match");
    
    matrix A(n,n);
    valarray<double> x(n),c(n);
	
    for(int i=0; i<n; i++)
    {
		c[i] = -plane[I[i]].c();
		for(int j=0; j<n; j++)
			A(i,j) = plane[I[i]][j];
    }
    
    if(A.det())
		x = A.solve(c);
	else
		return Point();
    
    return Point(x);
}

double HyperplaneSet::oja(const Point& x) const
{
    double sum = 0.0;
    double k = (1.0 / double(fact(dim())));
    
    for(int i=0; i<planes; i++)
		sum += k * fabs((plane[i] | x) + plane[i].c());
    
	
    return sum;
}

Point HyperplaneSet::gradient(const Point& x) const
{
    double k = (1.0 / double(fact(dim())));
	double sgn;
	Point Gr(dim());

    for(int i=0; i<planes; i++)
	{
		sgn = plane[i].side(x);
		Gr += sgn * k * plane[i].normal();
	}
	
	return Gr;
}

Point HyperplaneSet::oja_rank(const Point& x) const
{
	double sgn;
	Point Gr(dim());

    for(int i=0; i<planes; i++)
	{
		sgn = plane[i].side(x);
		Gr += sgn * plane[i].normal();
	}
	
	for (int i = 0; i < dim(); i++)
		Gr[i] /= planes;

	return Gr;
}

void HyperplaneSet::oja_and_gradient(const Point& x,double& oja
  ,Point& grad) const
{
    double k = (1.0 / double(fact(dim())));
	double sgn;
	Point Gr(dim());

	oja = 0.0;
    for(int i=0; i<planes; i++)
	{
		oja += k * fabs((plane[i] | x) + plane[i].c());
		sgn = plane[i].side(x);
		Gr += sgn * k * plane[i].normal();
	}
	
	grad = Gr;
}

ostream& operator <<(ostream& os,const HyperplaneSet& H)
{
    if(H.size()==0)
    {
		os << "empty set" << endl;
		return os;
    }
	
    for(int i=0; i<H.size(); i++)
		os << H[i] << endl;
    
    return os;
}
