/* $Id: point.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef POINT_H
#define POINT_H

#include <valarray>
#include <string>

using namespace std;

class Point
{    
    valarray<double> v;
    	
  public:
	
    Point() {}
    Point(const Point& p) : v(p.v)
		{}
    explicit Point(const valarray<double>& v0) : v(v0)
		{}
    explicit Point(int n);
	
    Point(double x,double y) : v(valarray<double>(2))
		{v[0]=x; v[1]=y;}
    Point(double x,double y,double z) : v(valarray<double>(3))
		{v[0]=x; v[1]=y; v[2]=z;}
    Point& operator=(const Point& p);
    Point& operator=(const string& s);
    Point& operator=(const char* s0)
		{string s; s=s0; return operator=(s);}
	
    int dim() const
		{return v.size();}
    bool is_nil() const
		{return !dim();}
    bool in_box(const Point& min,const Point& max) const;
	double length() const;
    double dist(const Point& p2) const
		{return Point(p2-(*this)).length();}
	double angle(const Point& p2) const; /* New */
    void normalize(); 
    double& operator[](int i)
		{return v[i];}
    double coord(int i) const
		{return v[i];}
    
    friend Point operator+(const Point& x1,const Point& x2);
    friend Point operator-(const Point& x1,const Point& x2);
    friend double operator|(const Point& x1,const Point& x2);
    friend Point operator*(double t,const Point& x2);
    friend bool operator==(const Point& x1,const Point& x2);
    friend bool operator!=(const Point& x1,const Point& x2);
    friend Point operator-(const Point& x);
    Point& operator+=(const Point& x)
		{v+=x.v; return *this;}
    Point& operator-=(const Point& x)
		{v-=x.v; return *this;}
    Point& operator*=(double t)
		{v = v * t; return *this;}
};

inline Point operator+(const Point& x1,const Point& x2) 
{valarray<double> tmp(x1.v); tmp+=x2.v; return Point(tmp);}
inline Point operator-(const Point& x1,const Point& x2) 
{valarray<double> tmp(x1.v); tmp-=x2.v; return Point(tmp);}
inline Point operator*(double t,const Point& x2) 
{valarray<double> tmp(x2.v); tmp*=t; return Point(tmp);}
bool operator==(const Point& x1,const Point& x2);
inline bool operator!=(const Point& x1,const Point& x2)
{return !operator==(x1,x2);}
inline Point operator-(const Point& x) 
{return -1.0 * x;}

class matrix;

matrix covariance(const Point& x,const Point& y);
Point operator*(const matrix& M,const Point& x);

ostream& operator <<(ostream& os,const Point& p);
istream& operator >>(istream& is,Point& p);
#endif
