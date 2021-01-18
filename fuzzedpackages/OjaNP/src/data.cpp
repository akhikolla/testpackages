/* 
 *  data.cpp : A collection of d-dimensional vectors.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: data.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "data.h"
#include "error.h"
#include "matrix_wrapper.h"
#include <math.h>
#include <fstream>
#include <algorithm>

using namespace std; //df

//  CLASS: Data
//  ===========

Data::Data()
{
    data = 0;
    dimension = 0;
}

Data::Data(const Data& S)
{
    data = 0;
    operator=(S);
}

Data& Data::operator=(const Data& S)
{
    if(data)
		delete data;
	
    dimension = S.dimension;
    
    if(S.data)
    {
		data = new vector<Point>(S.size());
		errif(!data,"Data::operator=: out of memory");
		for(int i=0; i<S.size(); i++)
			(*data)[i] = (*S.data)[i];
    }
    else
		data = 0;
    
    return *this;
}

Data::Data(int vecdim,int points)
{   
    errif(vecdim < 1,"Data::Data: dimension " << vecdim << " is less than 1");
    errif(points < 1,"Data::Data: number of points " 
	  << points << " is less than 1");
    
    Point v(vecdim);
	
    data = new vector<Point>(points);
    errif(!data,"Data::Data: out of memory");
    
    for(int i=0; i<points; i++)
		(*data)[i] = v;
    
    dimension = vecdim;
}

Data::Data(const char* filename)
{   
    data = 0;
    dimension = 0;
	
    enlarge(filename);
}

Data::~Data()
{
    if(data)
		delete data;
}

void Data::enlarge(int how_much)
{
    vector<Point> *newdata;
	
    errif(how_much < 0,"Data::enlarge: given value " << how_much << " was negative");
    errif(dimension==0,"Data::enlarge: dimension unknown");
    
    if(how_much==0)
		return;
    
    newdata = new vector<Point>(size()+how_much);
    errif(!newdata,"Data::enlarge: out of memory");
    Point v(dimension);
    
    for(int i=0; i<size(); i++)
		(*newdata)[i]=(*data)[i];
    for(int i=size(); i<size()+how_much; i++)
		(*newdata)[i]=v;
    
    if(data)
		delete data;
    
    data = newdata;
}

void Data::enlarge(const Point& p)
{
	if(size()==0)
		set_dim(p.dim());
	
	enlarge(1);
	(*this)[size()-1]=p;
}

void Data::enlarge(list<Point> &veclist)
{
    int oldsize;
    
    if(veclist.size() == 0)
		return;
    
    if(!dimension)
		dimension = veclist.front().dim();
	
    oldsize = size();
    enlarge(veclist.size());
	
    Point v;
    for(int i=0; veclist.size(); i++)
    {
		v = veclist.front();
		veclist.pop_front();
		errif(v.dim() != dimension,"Data::enlarge: vector has dimension "
	      << v.dim() << " but data set has " << dimension);
		
		(*this)[i+oldsize] = v;
    }
}
void Data::enlarge(const char* filename)
{
    ifstream F(filename);
    errif(!F,"Data::enlarge: file '" << filename << "' not found");
    F >> *this;
}

void Data::enlarge(const vector<Point> &Points){
	for (int i = 0; i<Points.size(); i++)
	{
		errif(Points[i].dim() != dimension, "Data::enlarge: vector has dimension "
			<< Points[i].dim() << " but data set has " << dimension);
	}
	data->insert(data->end(), Points.begin(), Points.end());
}

void Data::set_dim(int d)
{
    errif(dimension,"Data::set_dim: dimension already set");
    errif(d < 1,"Data::set_dim: illegal dimension " << d);
    
    dimension=d;
}

Point& Data::operator[](int index) const
{
    errif(size()==0,"Data::operator[]: no data");
    errif(index < 0 || index > size(),"Data::operator[]: illegal index " << index);
	
    return (*data)[index];
}

Point Data::average() const
{
    errif(size()==0,"Data::average: no data");
    
    Point sum(dim());

    for(int i=0; i<size(); i++)
		sum += (1.0 / double(size())) * (*data)[i];    
	
    return sum;
}

matrix Data::covariance() const
{
    errif(size()==0,"Data::covariance: no data");

	matrix C(dim(),dim());
	Point av=average();
	
	for(int k=0; k<size(); k++)
		for(int i=0; i<dim(); i++)
			for(int j=0; j<dim(); j++)
			{
				Point& x=operator[](k);
				
				C(i,j)+=(x.coord(i) - av[i])*(x.coord(j) - av[j]);
			}
	
	for(int i=0; i<dim(); i++)
		for(int j=0; j<dim(); j++)
			C(i,j)/=(size()-1);
	
	return C;
}

Point Data::min() const
{
    errif(!data,"Data::min: no data");
    
    Point mn=(*data)[0];
    
    for(int i=1; i<size(); i++)
		for(int j=0; j<mn.dim(); j++)
			if((*data)[i][j] < mn[j])
				mn[j] = (*data)[i][j];
	
    return mn;
}

Point Data::max() const
{
    errif(!data,"Data::max: no data");
    
    Point mx=(*data)[0];
	
    for(int i=1; i<size(); i++)
		for(int j=0; j<mx.dim(); j++)
			if((*data)[i][j] > mx[j])
				mx[j] = (*data)[i][j];
	
    return mx;
}

int Data::center_index() const
{
    errif(!data,"Data::center_index: no data");
    Point avg = average();
    double dmin,d;
    int index;
	
    dmin = avg.dist((*data)[0]);
    index = 0;
	
    for(int i=1; i<size(); i++)
    {
		d = avg.dist((*data)[i]);
		if(d < dmin)
		{
			dmin = d;
			index = i;
		}
    }
    
    return index;
}

Point Data::center() const
{
    errif(!data,"Data::center: no data");
	
    return (*data)[center_index()];
}

ostream& operator <<(ostream& os,const Data& S)
{
    for(int i=0; i<S.size(); i++)
    {
		for(int j=0; j<S.dim(); j++)
		{
			os << S[i][j];
			if(j < S.dim()-1)
				os << ' ';
		}
		os << endl;
    }
	
    return os;
}

istream& operator >>(istream& is,Data& S)
{
    Point p;
    list <Point> P;
	
    for(;;)
    {
		is >> p;
		if(p.is_nil())
			break;
		P.push_back(p);
    }
    S.enlarge(P);
    
    return is;
}

static Point v0;
int cmp_distance_from_origo(const Point& v1,const Point& v2)
{
    double d1,d2;
    
    d1 = (v1 - v0).length();
    d2 = (v2 - v0).length();
    if(d1 < d2)
		return -1;
    if(d1 > d2)
		return 1;
    return 0;
}

void Data::sort_by_distance(const Point& p)
{
    if(!size())
		return;

    v0 = p;
    sort(data->begin(),data->end(),&cmp_distance_from_origo);
}

void linear_fit(matrix& A,Point& z,const Data& X,const Data& Y)
{
	errif(X.size() != Y.size(),"linear_fit: data size mismatch");

	int d=X.dim(),N=X.size();

	z=Point(d);
 	A=matrix(d,d);

	Point Y0(N);
	Point Xy(d+1);
	matrix XX(d+1,d+1);

	for(int p=0; p<d; p++)
	{
		
		for(int i=0; i<N; i++)
			Y0[i]=Y[i][p];
		
		// Yksiulotteinen sovitus

		for(int i=0; i<d+1; i++)
			Xy[i]=0;

		for(int i=0; i<N; i++)
		{
			Xy[0]+=Y0[i];
			for(int j=0; j<d; j++)
				Xy[j+1]+=(X[i][j] * Y0[i]);
		}
				
		for(int i=0; i<d+1; i++)
			for(int j=0; j<d+1; j++)
				XX(i,j)=0.0;
		
		XX(0,0)=N;
		for(int i=0; i<d; i++)
		{
			for(int k=0; k<N; k++)
			{
				XX(0,i+1)+=X[k][i];
				XX(i+1,0)+=X[k][i];
			}
			
			for(int j=0; j<d; j++)
			{
				for(int k=0; k<N; k++)
				{
					XX(i+1,j+1)+=(X[k][i] * X[k][j]);
				}			
			}
		}

		Point bhat=XX.inv()*Xy;

		z[p]=bhat[0];
		for(int i=0; i<d; i++)
			A(p,i)=bhat[i+1];
	}
}
