/* $Id: hyperplane.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef HYPERPLANE_H
#define HYPERPLANE_H
#include "index.h"
#include "point.h"

using namespace std; //df

class Line;
class Data;

class Hyperplane
{
    int cofs;
    double* cof;
    
    void set_dim(int d);
	
  public:

	bool isBound;
    
    Hyperplane();
    ~Hyperplane();
    Hyperplane(const Hyperplane& H);
    Hyperplane& operator=(const Hyperplane& H);
    Hyperplane(int dim);   
	
    void get(const Data& D,const Index& I);
	void get(Point& point, Point& grad); // create custom hyperplane through[point] with normal = [grad] 
	void get(vector<Point> points); // create custom hyperplane by [dim] points

    double operator[](int index) const;
    int dim() const
		{return cofs ? cofs-1 : 0;}
    bool degenerated() const;

	double angle(const Hyperplane& H) const /* New */
		{return normal().angle(H.normal());}
    Point intersect(const Line& L) const;
    bool intersect(const Line& L,double& l) const;
    double side(const Point& x) const;
	double dist(const Point& x) const; 
    double operator|(const Point& x) const;
    Point normal() const;
    double c() const;
	Point cof_at(const Point& x) const; /* New */
	double cof0_at(const Point& x) const; /* New */

    void normalize(); 
    
    friend ostream& operator <<(ostream& os,const Hyperplane& H);    
};

const int MAX_BOUNDS = 30;

class HyperplaneSet
{
    Hyperplane* plane;
    int planes;
	int bounds;
	
    void resize(int new_size);
    
  public:
	
    HyperplaneSet(const HyperplaneSet& HS);
    HyperplaneSet()
		{plane=0; planes=0;}
    HyperplaneSet(int size);
    ~HyperplaneSet();
    
    HyperplaneSet& operator=(const HyperplaneSet& HS);
    Hyperplane& operator[](int i) const;
	
    int dim() const
		{return planes ? plane[0].dim() : 0;}
    int size() const
		{return planes;}
    void get(const Data& D,const IndexSet& I);
    void get_all(const Data& D);  
	void get(const vector<Hyperplane>& hyperplanes);
	void add(const Hyperplane& H);
    Point crossing_point() const;
    Point crossing_point(const Index& I) const;
    double oja(const Point& x) const;
    Point gradient(const Point& x) const;
	Point oja_rank(const Point& x) const;
	void oja_and_gradient(const Point& x,double& oja,Point& grad) const;
};

ostream& operator <<(ostream& os,const Hyperplane& H);
istream& operator >>(istream& is,Hyperplane& H);
ostream& operator <<(ostream& os,const HyperplaneSet& H);

#endif
