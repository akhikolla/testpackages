/* 
 *  oja_data.cpp : Source data for computing oja median.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: oja_data.cpp,v 1.1 2008/01/25 11:47:50 ruthe Exp $
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

#include <math.h>
#include <fstream>
#include "oja_geometry.h"
#include "simplex.h"
#include "error.h"

using namespace std; //df

#define SCALE_TRESHOLD 0.00000001
/* Mik�li skaalattaessa pisteit� v�lille [0,1]^n jollakin koordinaatilla ei */
/* esiinny kynnyst� ylitt�v�� vaihtelua, kyseist� koordinaattia ei skaalata. */

//  CLASS: OjaData
//  ==============

OjaData::OjaData()
{
    method=EVAL_ALL_POINTS;
	lattice_measure=LM_DIAMETER;
    plane=0;
    planes=0;
    planeindex=0;

	epsilon=0.1;
	chi2_limit=0.0;
	set_size=20;
}

OjaData::OjaData(const char* filename) : Data(filename)
{
    method=EVAL_ALL_POINTS; 
	lattice_measure=LM_DIAMETER;
    plane=0;
    planes=0;
    planeindex=0;

	string f=filename;

	if(f.length() <= 5)
		f=f+".oja";
	else
	{
		if(f.substr(f.length()-5,5)==".data")
			f=f.substr(0,f.length()-5);
		f=f+".oja";
	}

	ifstream oja(f.c_str());

	if(oja)
		oja >> exact_median;

	original_size = size();
}

OjaData::OjaData(int vecdim,int size) : Data(vecdim,size)
{
    method=EVAL_ALL_POINTS; 
	lattice_measure=LM_DIAMETER;
    plane=0;
    planes=0;
    planeindex=0;
	original_size = size;
}

OjaData::~OjaData()
{
    if(plane)
		delete plane;
    if(planeindex)
		delete[] planeindex;
}

void OjaData::generate_hyperplanes()
{
    if(plane)
		return;
	
	regenerate_hyperplanes();
}

void OjaData::regenerate_hyperplanes()
{
    if(plane)
		delete plane;
    if(planeindex)
		delete[] planeindex;
	
    plane = new HyperplaneSet;
    errif(!plane,"OjaData::generate_hyperplanes: out of memory");
    planes=hyperplanes();
    plane->get_all(*this);
	errif(planes != plane->size(),"FATAL:OjaData::generate_hyperplanes: "
	  << "size calculation do not work");
	
    // BUG: Taulukoinnin sijasta pit�isi pysty� laskemaan muunnos
	//      int -> Index. Samalla kaikista funktioista voisi poistaa
	//      erillisen if(planes) testin ja siirt�� sen hyperplane():een.
	int hyps=hyperplanes();
	planeindex = new Index[hyps + MAX_BOUNDS];
    errif(!planeindex,"OjaData::generate_hyperplanes: out of memory");
    
    Index I(dim(),size());
    for(int i=0; i<planes; i++,I++)
		planeindex[i] = I;
}

void OjaData::get_oja_and_gradient(const Point& x,double& oja
  ,Point& grad) const
{
    errif(x.dim() != dim(),"OjaData::get_oja_and_gradient: point " << x
	  << " has illegal dimension");

 	if(plane)
 	{
 		plane->oja_and_gradient(x,oja,grad);
 		return;
 	}
			
    double sum = 0.0;
    double sgn;
    double k = 1.0 / double(fact(dim()));
    Simplex S;
    Point Gr(dim());
	
    for(Index I(dim(),size()); I; I++)
    {
		S.get(*this,I,x);
		sum += S.size();
		sgn = S.sign();
		for(int j=0; j<dim(); j++)
			Gr[j] +=  S.row_cof(j+1) * sgn * k;
    }
    
    oja = sum;
    grad = Gr;
}

Point OjaData::gradient(const Point& x) const
{
    errif(x.dim() != dim(),"OjaData::gradient: Illegal dimension on point " << x);

 	if(plane)
 		return plane->gradient(x);

    Point Gr(dim());
    Simplex S;
    double sgn;
    double k = 1.0 / double(fact(dim()));
	
    for(Index I(dim(),size()); I; I++)
    {
		S.get(*this,I,x);
		sgn =  S.sign();
		for(int j=0; j<dim(); j++)
			Gr[j] +=  S.row_cof(j+1) * sgn * k;
    }
    
    
    return Gr;
}

Point OjaData::oja_rank(const Point& x) const
{
	errif(x.dim() != dim(), "OjaData::oja_rank: Illegal dimension on point " << x);

	if (plane)
		return plane->oja_rank(x);

	Point Gr(dim());
	Simplex S;
	double sgn;
	double k = 0;

	for (Index I(dim(), size()); I; I++, k++)
	{
		S.get(*this, I, x);
		sgn = S.sign();
		for (int j = 0; j<dim(); j++)
			Gr[j] += S.row_cof(j + 1) * sgn;
	}

	for (int i = 0; i < dim(); i++)
		Gr[i] /= k;

	return Gr;
}

double OjaData::oja(const Point& x) const
{
	if (size() == 0) return -1;
    errif(x.dim() != dim(),"OjaData::oja: Illegal dimension on point " << x);
	
    if(plane)
		return plane->oja(x);
	
    double sum = 0.0;
    Simplex S;
	
    for(Index I(dim(),size()); I; I++)
    {
		S.get(*this,I,x);
		sum += S.size();
    }
    
    return sum;
}

bool OjaData::derivative(const Point& x,const Point& direction,double& Dpos
  ,double& Dneg) const
{
    Point d=direction;
    d.normalize();
    
	Index I(dim(),size());
	Hyperplane H;
    double D; // Apumuuttuja yhden termi laskemiseen
    double sum=0.0; // Summa positiiviseen suuntaan
	double sum_=0.0; // Summa negatiiviseen suuntaan
    double k=1.0 / double(fact(dim()));
    double sgn;

    for(int i=0; i<hyperplanes(); i++)
    {
		if(planes)
			H=hyperplane(i);
		else
		{
			H.get(*this,I);
			I++;
		}
		
		H.normalize();
		D = k * (H|d);       
		sgn = H.side(x);
		if(sgn > 0)              // Piste ei ole hypertasolla
			sum += D, sum_ -= D;
		else if (sgn < 0)
			sum -= D, sum_ += D;
		else                     // Piste on hypertasolla
			sum += fabs(D), sum_ += fabs(D);
    }

    Dpos = sum;
    Dneg = sum_;
    
    return Dpos >= 0.0 && Dneg >= 0.0;
}

const Hyperplane& OjaData::hyperplane(int i) const
{
    errif(!plane,"OjaData::hyperplane: hyperplanes not yet generated");
    errif(i < 0 || i >= hyperplanes(), "OjaData::hyperplane: illegal index "
	  << i);
    
    return (*plane)[i];
}

Index OjaData::hyperplaneindex(int i) const
{
    errif(!plane,"OjaData::hyperplaneindex: hyperplanes not yet generated");
    errif(i < 0 || i >= hyperplanes(), "OjaData::hyperplaneindex: illegal index "
	  << i);
	
    return planeindex[i];
}

const HyperplaneSet& OjaData::hyperplaneset() const
{
    errif(!plane,"OjaData::hyperplaneset: hyperplanes not yet generated");
	
    return *plane;
}

Point OjaData::scaled(Point x) const
{
    if(norm_offset.is_nil())
		return x;
    
    errif(x.dim() != dim(),"OjaData::scaled: dimensions do not match");
    
    x -= norm_offset;
    for(int i=0; i<x.dim(); i++)
		x[i] /= norm_scale.coord(i);    
	
    return x;
}

Point OjaData::descaled(Point x) const
{
    if(norm_offset.is_nil())
		return x;
    
    errif(x.dim() != dim(),"OjaData::descaled: dimensions do not match");
	
    for(int i=0; i<x.dim(); i++)
		x[i] *= norm_scale.coord(i);    
	
    x += norm_offset;
    
    return x;
}

void OjaData::scale()
{
    errif(size()==0,"OjaData::normalize: no data");
	
    norm_offset=min();
    norm_scale=(max()-norm_offset);
    for(int i=0; i<dim(); i++)
		if(fabs(norm_scale[i]) < SCALE_TRESHOLD)
			norm_scale[i] = 1.0;
    for(int i=0; i<size(); i++)
		operator[](i) =  scaled(operator[](i));
}

Point OjaData::min() const
{
	if (boundedMin.dim() != 0)
		return boundedMin;

	return (Data::min());
}

Point OjaData::max() const
{
	if (boundedMax.dim() != 0)
		return boundedMax;

	return (Data::max());
}

void OjaData::set_bounded_min_max(const Point bmin, const Point bmax)
{
	boundedMin = min();
	boundedMax = max();

	double enlarge = 1.0e-8;

	for (int j = 0; j < dim(); j++)
	{
		if (boundedMax[j] >= bmax.coord(j))
			boundedMax[j] = bmax.coord(j) + enlarge;
		if (boundedMin[j] <= bmin.coord(j))
			boundedMin[j] = bmin.coord(j) - enlarge;
	}
}

void OjaData::add_bound_points(const vector<Point> & crossing_points){
	enlarge(crossing_points);

	/*for (int i = 0; i < planes; i++){
		planeindex[i].set_limit(size());
	}*/
}

void OjaData::add_bound(const Hyperplane& bound, const set<int>& crossings){
	plane->add(bound);
	includedPlanes.insert(planes);
	planeindex[planes] = Index(dim(), size(), crossings);
	planes++;
}
