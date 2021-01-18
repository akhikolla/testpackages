/* $Id: data.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#define R_NO_REMAP

#ifndef DATA_H
#define DATA_H
#include "point.h"
#include <vector>
#include <list>

using namespace std;

class Data
{
    int dimension;
    vector<Point> *data;
	
protected:
    
    void set_dim(int d);
    
public:
    
    Data();
    Data(const Data& S);
    Data& operator=(const Data& S);
    
    Data(int vecdim,int size);
    Data(const char* filename);
    ~Data();
	
    int size() const
		{return data==0 ? 0 : (*data).size();}
    int dim() const
		{return dimension;}
	
    Point& operator[](int index) const;
    Point average() const;
	matrix covariance() const;
    int center_index() const;
    Point center() const;   
    virtual Point min() const;
	virtual Point max() const;
	
    void enlarge(int n);
    void enlarge(list<Point> &Points);
	void enlarge(const vector<Point> &Points);
    void enlarge(const char* filename);
	void enlarge(const Point& point);
	
    void sort_by_distance(const Point& v);
};

// PNS sovitus s.e. y[i] = A x[i] + b
void linear_fit(matrix& A,Point& b,const Data& x,const Data& y);

ostream& operator <<(ostream& os,const Data& S);
istream& operator >>(istream& is,Data& S);
#endif
