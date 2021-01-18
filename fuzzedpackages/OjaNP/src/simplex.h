/* $Id: simplex.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */
#ifndef SIMPLEX_H
#define SIMPLEX_H
#include "data.h"
#include "index.h"
#include "misc.h"
#include "matrix_wrapper.h"

using namespace std; //df

class Simplex
{
    matrix M;
	
    friend ostream& operator <<(ostream& os,const Simplex& S);

	double det() const;
	
  public:
	
    Simplex();
    Simplex(int dim);
	~Simplex();
    Simplex(const Simplex& S);
	Simplex& operator=(const Simplex& S);
	
    void get(const Data& D,const Index& I,const Point& v);
    void get(const Data& D,const Index& I);
	void get(const vector<Point> points);
    void set_column(int col,const Point& v);
	
    int dim() const
		{return M.rows() ? M.rows()-1 : 0;}
    double size() const;
    double sign() const;
    double row_cof(int r) const
		{return cof(M,r,dim());}

	bool contains(const Point& p);
};

ostream& operator <<(ostream& os,const Simplex& S);

#endif
