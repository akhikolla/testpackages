/* $Id: line.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef LINE_H
#define LINE_H
#include "point.h"
#include "index.h"
#include "hyperplane.h"
#include <iostream>

class Data;

using namespace std;

class Line
{
    friend ostream& operator <<(ostream& os,const Line& L);
    
    Point start,direction;
    
    void set_dim(int d);
	
  public:
    
    Line() {}
    Line(const Line& l2)
		{start=l2.start; direction=l2.direction;}
    Line& operator=(const Line& l2)
		{start=l2.start; direction=l2.direction; return *this;}
    Line(const Point& p0,const Point& d);
    Line(int d)
		{set_dim(d);}
	
    int dim() const {return start.dim();}
    void get(HyperplaneSet& H);
	void get(const Data& D,const IndexSet& I);
    void get_through(const Point& p1,const Point& p2); 
    void get_toward(const Point& p0,const Point& dir); 
    Point intersect(const Hyperplane& H) const
		{return H.intersect(*this);}
    bool intersect(const Hyperplane& H,double &l) const
		{return H.intersect(*this,l);}
    double x0(int i) const;
    double dx(int i) const;
    Point dir() const
		{return direction;} 
	Point origin() const
		{return start;}
    double angle(const Point& z) const;
    Point proj(const Point& z) const; 
    double proj_length(const Point& z) const
		{return proj(z).length();} 
    Point at(double t) const;
    bool is_nil() const;
};

ostream& operator <<(ostream& os,const Line& L);

#endif
