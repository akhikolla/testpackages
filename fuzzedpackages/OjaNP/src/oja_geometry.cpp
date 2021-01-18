/* 
 *  oja_geometry.cpp : Geometry with indices.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: oja_geometry.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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

#include "oja_geometry.h"
#include <math.h>
#include <algorithm>

using namespace std; //df

// Defines
// =======

#define OJA_FAVOUR 0.0
/* Testattaessa onko piste varmasti minimi sen kautta kulkevien suorien */
/* avualla, annetaan t�m�n verran etua nykyiselle pisteelle. Muutoin py�ristys- */
/* virhe voi hyl�t� oikean pisteen. */

/* Funktion minimointi suoralla (valitse yksi) */
/* ---------------------------- */

/*#define FNMIN_EVALUATE_ALL */
/* Muodostetaan kaikki pisteet suoralla ja evaluoidaan */

#define FNMIN_SWAP_FACTORS
/* K�ytet��n suoraa evaluointia vierekk�isten hypertasojen kertoimien */
/* vaihtotekniikkaa apuna k�ytt�en */

//  CLASS: OjaPoint
//  ===============

void OjaPoint::set_index(const IndexSet& I)
{
	errif(!data,"OjaPoint::set_index: no data");
    errif(I.indices() != I.dim(),"OjaPoint::set_index: cannot construct from " 
	  << I);
	errif(data->dim() != I.dim(),"OjaPoint::set_index: dimension " << I.dim()
	  << " do not match with data");
    errif(I.limit() > data->size(),"OjaPoint::set_index: maximum value "
	  << I.limit() << " is too big for data");
	     
    id = I;
}

void OjaPoint::set_index(const IndexSet& I,const Index& i)
{
	errif(!data,"OjaPoint::set_index: no data");
    errif(I.indices() != i.dim()-1,"OjaPoint::set_index: cannot construct"
	  "from " << I << " and " << i);
    errif(I.limit() != i.limit(),"OjaPoint::set_index: IndexSet and Index has"
	  " different maximum values");
    errif(I.dim() != i.dim(),"OjaPoint::set_index: dimensions " << I.dim()
	  << " and " << i.dim() << " do not match");
	errif(data->dim() != I.dim(),"OjaPoint::set_index: dimension " << I.dim()
	  << " do not match with data");
    errif(I.limit() > data->size(),"OjaPoint::set_index: maximum value "
	  << I.limit() << " is too big for data");
    
	int parts = I.indices();
	IndexSet new_id(parts+1, i.dim(), i.limit());
	for (int j = 0; j<parts; j++)
		new_id[j] = I[j];
	new_id[parts] = i;
    new_id.validate();
    
    id = new_id;
}

int OjaPoint::data_index() const
{
	errif(!data,"OjaPoint::data_index: no data");
    int idx=id.common_digit();
    errif(idx < 0,"OjaPoint::data_index: point " << *this << " is not"
	  " a data point");
    return idx;
}

bool OjaPoint::splits_line(int& n1,int &n2,Index& H) const
{
	errif(!data,"OjaPoint::splits_line: no data");
	n1 = -1;
	n2 = -1;

	if(data->dim() <= 2)
		return false;

	IndexSet I(data->dim()-1,data->dim(),data->size());

	for(int i=0; i<I.dim(); i++)
	{
		for(int j=0; j<I.indices(); j++)
			I[j] = id[j+(j >= i)];

		if(I.has_two_common_digits(n1,n2))
		{
			H = id[i];

			return true;
		}
	}
	
	return false;
}

bool OjaPoint::better_route(OjaLine& L,OjaPoint& P,double &o2)
{
	errif(!data,"OjaPoint::better_route: no data");
	double o1;

	o1 = data->oja(location())-OJA_FAVOUR;
	L = scan_all_routes(P,o2);

	return o2 < o1;
}

void OjaPoint::get(const IndexSet& I)
{
    errif(!data,"OjaPoint::get:no data");
    errif(I.indices() != data->dim(),"OjaLine::get: wrong number of indices " 
	  << I.indices());
    errif(I.dim() != data->dim() ,"OjaLine::get: illegal index dimension " 
	  << I.dim());

	HyperplaneSet H;
	
	H.get(*data,I);

	Point::operator=(H.crossing_point());
    id = I;
    id.validate();
}

OjaLine OjaPoint::scan_all_routes(OjaPoint& p,double& ojafn,double hi_score)
{
	OjaLine L(data);

	OjaLineIterator I(*this,false,false);
	double o;
	OjaPoint p0(data);

	ojafn=hi_score;

	while(I)
	{
		p0=I.line().min(o);
		if(!p0.is_nil())
		{
			if(ojafn==-1.0 || o < ojafn)
			{
				ojafn=o;
				p=p0;
				L=I.line();
				if(hi_score != -1.0)
					return L;
			}
		}
		
		I++;
	}
	
	return L;
}

ostream& operator <<(ostream& os,const OjaPoint& p)
{
	IndexIdentifier ID;

	ID.get(p.index());
	
	os << ID;
    return os;
}

// CLASS: DotSet
// =============

DotSet::DotSet(const OjaLine& L)
{
	line=&L;
	generate_dots();
}

DotSet::DotSet(const OjaLine* L)
{
	line=L;
	generate_dots();
}

DotSet::DotSet(const OjaLine* L, Point& h, double& h0)
{
	line = L;
	get_common_coefs(h, h0);
}

DotSet::~DotSet()
{
	if(sorted)
		delete dotarray;
}

void DotSet::generate_dots()
{
	if (Data()->get_includedPlanes().size() != 0){
		generate_dots_bounded();
		return;
	}

	sorted = false;

	h = Point(dim());
	h0 = 0;
	
    Point min=Data()->min();
    Point max=Data()->max();

#ifdef BOX_ENLARGEMENT
	// Eliminoidaan py�ristysvirheist� johtuva pisteiden katoaminen
	// suurentamalla laatikkoa
	for(int i=0; i<dim(); i++)
	{
		min[i] -= BOX_ENLARGEMENT;
		max[i] += BOX_ENLARGEMENT;
	}
#endif
		
	double t;
	pair<double,int> d0;
	Line L;
	Point p,x;

	L=line->line();

	// Etsit��n yksi piste 'x' laatikon sis�lt�.
	// BUG: t�m� on huono toetutus
    for(int i=0; i<Data()->hyperplanes(); i++)
		if(L.intersect(Data()->hyperplane(i),t))
		  {		  
			p = point_at(t);
			if(p.in_box(min,max))
			{
				x = p;
				break;
			}
		  }

	errif(x.is_nil(),"DotSet::generate_dots: no points at box");

    for(int i=0; i<Data()->hyperplanes(); i++)
    {
		if(L.intersect(Data()->hyperplane(i),t))
		{
			p = point_at(t);
			if(p.in_box(min,max))
			{
				d0.first = t;
				d0.second = i;
				dotlist.push_back(d0);
			}
			else
			{
				h += Data()->hyperplane(i).cof_at(x);
				h0 += Data()->hyperplane(i).cof0_at(x);
			}
		}
		else
		{
			h += Data()->hyperplane(i).cof_at(x);
			h0 += Data()->hyperplane(i).cof0_at(x);
		}
    }
}

/*
Searching the borders of the bounded region along a line
Copyright(C) 2015 Oleksii Pokotylo
*/
set<int> DotSet::find_valid_bounds(set<int>& includedPlanes, const OjaData* data, Point& x)
{
	double t;
	pair<double, int> d0;
	Line L = line->line();
	set<int> valid_bounds;
	x = Point();
	
    Point min=Data()->min();
    Point max=Data()->max();

	for (set<int>::reverse_iterator it = includedPlanes.rbegin(); it != includedPlanes.rend(); ++it){
		if (!data->hyperplane(*it).isBound) break;
		if (!L.intersect(data->hyperplane(*it), t)) continue;

		d0.first = t;
		d0.second = *it;
		dotlist.push_back(d0);
	}
	sort();
	if (dotarray->size() == 0){
		return set<int>();
	}
	else
	if (dotarray->size() == 2){
		valid_bounds.insert((*dotarray)[0].second);
		valid_bounds.insert((*dotarray)[1].second);
	}
	else {
		int i;
		Point ref = point(dot(0));
		for (i = 0; i < dotarray->size() - 1; i++){
			if (hyperplane(dot(i + 1)).side(ref) > 0){
				// the point on this bound is also on the positive side of the next bound
				if (!point(dot(i)).in_box(min, max) || !point(dot(i + 1)).in_box(min, max)){
					return set<int>();
				}
				valid_bounds.insert((*dotarray)[i].second);
				valid_bounds.insert((*dotarray)[i + 1].second);
				x = point(dot(i));
				break;
			}
		}
		for (i = i + 1; i < dotarray->size() - 1; i++){
			if (hyperplane(dot(i + 1)).side(ref) < 0){
				// the line lies outside the bounded region
				x = Point();
				valid_bounds.clear();
				break;
			}
		}
	}

	sorted = false;
	dotlist.clear();
	delete dotarray;

	return valid_bounds;
}

/*
Getting the coefficients' sums of the hyperplanes outside the bounded region
Copyright(C) 2015 Oleksii Pokotylo
*/
void DotSet::get_common_coefs(Point& h, double& h0)
{
	sorted = false;

	const OjaData* data = Data();
	h = Point(dim());
	h0 = 0;

	double t;
	pair<double, int> d0;
	Line L;
	Point x;

	L = line->line();

	set<int> includedPlanes = data->get_includedPlanes();
	int hyperplanes_count = data->hyperplanes();
	set<int> valid_bounds = find_valid_bounds(includedPlanes, data, x);

	if (valid_bounds.size() == 0){
		h = Point();
		return;	// the line lies out of bounds
	}

	errif(x.is_nil(), "DotSet::generate_dots: no points at box");

	for (int i = 0; i < hyperplanes_count; i++)
	{
		if (!includedPlanes.count(i))
		{
			const Hyperplane& hi = data->hyperplane(i);
			h += hi.cof_at(x);
			h0 += hi.cof0_at(x);
		}
	}
}

/*
The modification of the generating procedure for the bounded search
Copyright(C) 2015 Oleksii Pokotylo
*/
void DotSet::generate_dots_bounded()
{
	sorted = false;

	const OjaData* data = Data();
	Point min = data->min();
	Point max = data->max();
	h = data->h;
	h0 = data->h0;

#ifdef BOX_ENLARGEMENT
	// Eliminoidaan py�ristysvirheist� johtuva pisteiden katoaminen
	// suurentamalla laatikkoa
	for (int i = 0; i < dim(); i++)
	{
		min[i] -= BOX_ENLARGEMENT;
		max[i] += BOX_ENLARGEMENT;
	}
#endif

	double t;
	pair<double, int> d0;
	Line L;
	Point p, x;

	L = line->line();
	IndexSet lind = line->index();

	set<int> includedPlanes = data->get_includedPlanes();
	int hyperplanes_count = data->hyperplanes();
	set<int> valid_bounds = find_valid_bounds(includedPlanes, data, x);
	if (valid_bounds.size() == 0) return;


	errif(x.is_nil(), "DotSet::generate_dots: no points at box");

	for (set<int>::iterator it = includedPlanes.begin(); it != includedPlanes.end(); ++it)
	{
		int i = *it;
		const Hyperplane& hi = data->hyperplane(i);
		if (L.intersect(hi, t))
		{
			if (hi.isBound && !valid_bounds.count(i))
				continue;
			p = point_at(t);
			if (p.in_box(min, max) && !lind.has(data->hyperplaneindex(i)))
			{
				d0.first = t;
				d0.second = i;
				dotlist.push_back(d0);
				continue;
			}
		}
		if (!hi.isBound){
			h += hi.cof_at(x);
			h0 += hi.cof0_at(x);
		}
	}
}

void DotSet::sort()
{
	if(sorted)
		return;
	
	dotarray=new vector<pair<double,int> >(size());
	errif(!dotarray,"DotSet::sort: out of memory");
	
	int n=0;
	list<pair<double,int> >::iterator i;
	for(i=dotlist.begin(); i!=dotlist.end(); i++)
		(*dotarray)[n++] = *i;

	::sort(dotarray->begin(),dotarray->end());
	dotlist.clear();

	sorted = true;
}

pair<double,int>& DotSet::dot(int index) const
{
	errif(!sorted,"DotSet::dot: dot set is not sorted");
	errif(index < 0 || index >= size(),"DotSet::dot: illegal index " << index);
	
	return (*dotarray)[index];
}

double DotSet::oja(double at) const
{
	double sum=0.0;
    double k = (1.0 / double(fact(Data()->dim())));
	Point x = point_at(at);
	// BUG: Oikeastaan t�ss� pit�isi olla tarkistus:
	// errif(!x.in_box(Data()->min(),Data()->max()),"Err");

	if(sorted)
	{
		for(int i=0; i<size(); i++)
			sum += k * fabs((hyperplane(i)|x) + hyperplane(i).c());
	}
	else
	{
		list<pair<double,int> >::const_iterator i;
		for(i=dotlist.begin(); i!=dotlist.end(); i++)
			sum += k * fabs((hyperplane(*i)|x) + hyperplane(*i).c());
	}

	sum += (h | x) + h0;
	
	return sum;
}

OjaPoint DotSet::min_evaluate_all(double &ojafn)
{
	OjaPoint p(Data());

	if(size()==0)
		return p;

    double min=0.0,o;
 	pair<double,int> mt;
 	mt.second=-1;

	if(sorted)
	{
		for(int i=0; i<size(); i++)
		{
			o = oja(dot(i));
			
			if(mt.second==-1 || o < min)
			{
				mt = dot(i);
				min = o;
			}
		}
	}
	else
	{
		list<pair<double,int> >::const_iterator i;
		for(i=dotlist.begin(); i!=dotlist.end(); i++)
		{
			o = oja(*i);
			
			if(mt.second==-1 || o < min)
			{
				mt = *i;
				min = o;
			}
		}
	}
	
	p.set_index(line->index(),Data()->hyperplaneindex(mt.second));
	p.set_location(point(mt));
	ojafn=min;

	return p;
}

OjaPoint DotSet::min(double& ojafn)
{
#ifdef FNMIN_SWAP_FACTORS
	errif(size()==0,"DotSet::min: dot set empty");
	sort();

	// Hypertasojen m��r�
	int n=size(); 

	// Ensimm�ist� ja viimeist� pistett� k�ytet��n etumerkin
	// m��r��miseen siirrelt�ess� summasta hypertasoja.
	Point x0=point(0);
	Point xn=point(n-1);

	// K�sitelt�v� leikkauspiste.
	Point x=x0;


	// Kulloinenkin minimi ja funktion arvo minimiss�.
	int min;
	double omin;

	// Kumulatiivinen summa hypertasojen kertoimista.
	Point gi(dim());
	double g0=0.0;

	// Apumuuttuja.
	double o;

	int start = 1;

	if (Data()->get_includedPlanes().size() == 0){

	// Lasketaan ensin summa aloituspisteess� ja otetaan se
	// pohjatulokseksi.   
	for(int i=0; i<n; i++)
	{
		const Hyperplane& h = hyperplane(dot(i));
		if (h.isBound)
			continue;
		gi += h.cof_at(x0);
		g0 += h.cof0_at(x0);
	}

	}
	else {
		bool had_bound = false;
		x0 = xn;
	// Lasketaan ensin summa aloituspisteess� ja otetaan se
	// pohjatulokseksi.   
	for (int i = 0; i < n; i++)
	{
		const Hyperplane& hp = hyperplane(dot(i));
		if (hp.isBound){
			if (!had_bound){
				x0 = point(0);
				x = point(i);
				min = i;
				start = i + 1;
				had_bound = true;
			}
			continue;
		}
		gi += hp.cof_at(x0);
		g0 += hp.cof0_at(x0);
	}
	}


	omin=((h | x)+h0 + (gi | x)+g0);
#ifdef _MSC_VER 
#ifdef DEEPDEBUG
	double real = OjaData::S.oja(x);
#endif 
#endif
	
	bool f1 = false, f2 = false;
	// K�yd��n l�pi kaikki leikkauspisteet vaihtaen hypertasoja
	// summasta oikean merkkisiksi.
	for (int i = start; i<n; i++)
	{	
		const Hyperplane& hp = hyperplane(dot(i - 1));
		if (!hp.isBound)
		{
			gi -= hp.cof_at(x0);
			g0 -= hp.cof0_at(x0);

			gi += hp.cof_at(xn);
			g0 += hp.cof0_at(xn);
		}
		else {
#ifdef _MSC_VER 
#ifdef DEEPDEBUG
			real = OjaData::S.oja(x);
#endif 
#endif
			if (i != start) break; 
		}

		x=point(i);
		o=((h | x)+h0 + (gi | x)+g0);

		if(o<omin)
		{
			omin=o;
			min=i;
		}
	}

	OjaPoint p(Data());
	
	p.set_index(line->index(),Data()->hyperplaneindex(dot(min).second));
	p.set_location(point(dot(min)));
	ojafn=omin;
	//cout << p << "; " << line->index() << "; " << Data()->hyperplaneindex(dot(min).second) << endl;
	return p;
#endif
	
#ifdef FNMIN_EVALUATE_ALL
	return min_evaluate_all(ojafn);
#endif
}
	
//  CLASS: OjaLine
//  ==============

OjaLine::OjaLine() : Line()
{
	data=0;
}

OjaLine::OjaLine(const OjaLine& DL) : Line(DL)
{
    data = DL.data;
    id = DL.id;
}

OjaLine::OjaLine(const OjaData& D) 
{
    data = &D;
}

OjaLine::OjaLine(const OjaData* D) 
{
    data = D;
}

OjaLine& OjaLine::operator=(const OjaLine& DL)
{
    Line::operator=(DL);
    data = DL.data;
    id = DL.id;
    
    return *this;
}

OjaLine& OjaLine::operator=(const Line& DL)
{
    Line::operator=(DL);
    
    return *this;
}

int OjaLine::dim() const
{
	errif(!data,"OjaLine::dim: no data");
	return data->dim();
}

// BUG: t�m�n voisi tehd� ehk� suoremminkin
void OjaLine::get(const IndexSet& I)
{
    errif(!data,"OjaLine::get:no data");
    errif(I.indices()+1 != data->dim(),"OjaLine::get: wrong number of indices " 
	  << I.indices());
    errif(I.dim() != data->dim() ,"OjaLine::get: illegal index dimension " 
	  << I.dim());
	
	Line::get(*data,I);
    id = I;
    id.validate();
}

void OjaLine::set(const IndexSet& I,const Line& L)
{
	Line::operator=(L);
	id=I;
}

OjaPoint OjaLine::crossing_point(const Index& hyperplane) const
{
    errif(!data,"OjaLine::get:no data");

	OjaPoint p(data);

	p.set_index(id,hyperplane);
	p.get(p.index());

	return p;
}

void OjaLine::get_random_through(int index)
{ 
    errif(!data,"OjaLine::get_random_through:no data");
    errif(index < 0 || index >= data->size(),"OjaLine::get_random_through:"
	  " illegal index " << index);
	
    Line L;
    HyperplaneSet H(data->dim()-1);
    IndexSet Hindex(data->dim()-1,data->dim(),data->size());
    Index I(data->dim(),data->size());
	
    int last_try = 1000;
    do 
    {
		errif(last_try==0,"OjaData::randomline_through: giving up - not enough data points");
		last_try--;
		
		for(int i=0; i<data->dim()-1; i++)
		{
			I.random_with(index);
			Hindex[i]=I;
			H[i].get(*data,I);
		}
		L.get(H);
		
    } while(L.is_nil());    
    id = Hindex;
    id.validate();
    
    *this = L;
}

void OjaLine::get_random_through(int index1,int index2)
{
    errif(!data,"OjaLine::get_random_through:no data");
    errif(index1 < 0 || index1 >= data->size(),"OjaLine::get_random_through:"
	  " illegal index " << index1);
    errif(index2 < 0 || index2 >= data->size(),"OjaLine::get_random_through:"
	  " illegal index " << index2);
	
    Line L;
    HyperplaneSet H(data->dim()-1);
    IndexSet Hindex(data->dim()-1,data->dim(),data->size());
    Index I(data->dim(),data->size());
	
    int last_try = 1000;
    do 
    {
		errif(last_try==0,"OjaData::randomline_through: giving up - not enough data points");
		last_try--;
		
		for(int i=0; i<data->dim()-1; i++)
		{
			I.random_with(index1,index2);
			Hindex[i]=I;
			H[i].get(*data,I);
		}
		L.get(H);
		
    } while(L.is_nil());    
    id = Hindex;
    id.validate();
    
    *this = L;
}

OjaPoint OjaLine::min(double& ojafn) const
{
    errif(!data,"OjaLine::min:no data");

	DotSet DS(this);

	if (DS.size() == 0)
		return OjaPoint();

	double f;
	OjaPoint q=DS.min(f);

	ojafn=f;
    return q;
}

bool OjaLine::is_same(const OjaLine& L2)
{
	int n1,n2;
	if(index().has_two_common_digits(n1,n2))
	{
		int m1,m2;
		if(L2.index().has_two_common_digits(m1,m2))
		{
			if((n1==m1 && n2==m2) ||
			  (n2==m1 && n1==m2))
				return true;
		}
	}

	return false;
}

ostream& operator <<(ostream& os,const OjaLine& L)
{
// 	int n1,n2;

// 	if(L.is_nil())
// 		os << "{nil}";
// 	else if(L.dim() >= 3 && L.index().has_two_common_digits(n1,n2))
// 	{
// 		if(L.index().how_many_common_digits() > 2)
// 			os << "{N/A " << L.index() << '}';
// 		else
// 			os << '{' << n1 << '-' << n2 << '}';
// 	}
// 	else
// 		os << '{' << L.index() << '}';
	IndexIdentifier ID;

	ID.get(L.index());
	os << ID;
	
    return os;
}

istream& operator >>(istream& is,OjaLine& L)
{
    return is;
}

//  CLASS: OjaLineSet
//  =================

OjaLine& OjaLineSet::operator[](int i)
{
    errif(line.empty(),"OjaLineSet::next: list is empty");
    errif(i < 0 || i>=int(line.size()),"OjaLineSet::next: illegal index " << i);
    
    if(last_access_index==-1 && i==0)
    {
		last_access = line.begin();
		last_access_index = 0;
    }
    else if(i==last_access_index)
    {
    }
    else if(i==last_access_index+1)
    {
		last_access++;;
		last_access_index = i;
    }
    else
    {
		err("OjaLineSet::operator[]: only sequential access is allowes");
    }

	return *last_access;
}

// BUG: kaikista n�ist� validisuustesti monesta samasta indeksist�
// pit�isi siirt�� muualle
void OjaLineSet::make_all_combinations(int with)
{
	IndexSet I(data->dim()-1,data->dim()-1,data->size());
	IndexSet J(data->dim()-1,data->dim(),data->size());
	OjaLine L(*data);
	
	clear();
		
	while(I)
	{
		for(int i=0; i<data->dim()-1; i++)
		{
			J[i][0]=with;
			for(int j=0; j<data->dim()-1; j++)
				J[i][j+1]=I[i][j];
			if(J.validate() && J.how_many_common_digits() <= 2)
			{
				L.get(J);
				add(L);
			}
		}
		I++;
	}		
}

void OjaLineSet::make_data_combinations(int with)
{
    OjaLine L(*data);
    
    clear();

    for(int i=0; i<data->size(); i++)
		if(i != with)
		{
			L.get_random_through(with,i);
			add(L);
		}
}

void OjaLineSet::make_combinations(const IndexSet& I)
{
    errif(I.dim() != data->dim(),"OjaLineSet::make_combinations: illegal dimension");
    OjaLine L(*data);
    IndexSet J(I.indices()-1,I.dim(),I.limit()); 
    
    clear();

    for(int i=0; i<I.indices(); i++)
    {
		for(int j=0; j<J.indices(); j++)
			J[j] = I[j + (j >= i)];
		L.get(J);
		add(L);
    }
}

void OjaLineSet::make_split_line_combinations(int n1,int n2,const Index& H)
{
	errif(H.dim() != data->dim(),"OjaLineSet::make_split_line_combinations: illegal"
	  " dimension" << H.dim());
	errif(data->dim() < 3,"OjaLineSet::make_split_line_combinations: dimension "
	  << data->dim() << " too low for operation");
	errif(n1 < 0 || n1 >= data->size(),"OjaLineSet::make_split_line_combinations: n1="<<n1);
	errif(n1 < 0 || n2 >= data->size(),"OjaLineSet::make_split_line_combinations: n2="<<n2);
	
	OjaLine L(*data);
	clear();
	int d=data->dim();
	IndexSet I(d-1,d,data->size());
	IndexSet J(d-2,d-2,data->size());

	while(J)
	{
		for(int i=0; i<d-2; i++)
		{
			I[i][0]=n1;
			I[i][1]=n2;
			for(int j=0; j<d-2; j++)
				I[i][j+2] = J[i][j];

		}
		I[d - 2] = H;
		if(I.validate() && I.how_many_common_digits() <= 2)
		{
			L.get(I);
			add(L);
		}

		J++;
	}
	
// Lis�t��n viel� suora n1-n2

	L.get_random_through(n1,n2);
	add(L);
}

void OjaLineSet::make_full_combinations(const OjaPoint& p)
{
	int n1,n2;
	Index I;

	if(p.is_data())
		make_all_combinations(p.data_index());
	else if(p.splits_line(n1,n2,I))
		make_split_line_combinations(n1,n2,I);
	else
		make_combinations(p.index());
}

void OjaLineSet::make_combinations(const OjaPoint& p)
{
#ifdef STRONG_COMBINATORICS	
	make_full_combinations(p);
#endif

#ifdef NORMAL_COMBINATORICS
	int n1,n2;
	Index I;

	if(p.is_data())
		make_data_combinations(p.data_index());
	else if(p.splits_line(n1,n2,I))
		make_split_line_combinations(n1,n2,I);
	else
		make_combinations(p.index());
#endif

#ifdef WEAK_COMBINATORICS
	if(p.is_data())
		make_data_combinations(p.data_index());
	else
		make_combinations(p.index());
#endif
}

int OjaLineSet::scan_all_routes(OjaPoint& best,double& ojafn)
{
	errif(!size(),"OjaLineSet::scan_all_lines: lineset empty");
	
	int min_idx=-1;
	double o;
	OjaPoint p0(data);
	ojafn=-1.0;
	
	for(int i=0; i<size(); i++)
	{
		p0=(*this)[i].min();

		if(!p0.is_nil())
		{
			o=data->oja(p0.location());
			if(ojafn==-1.0 || o < ojafn)
			{
				ojafn = o;
				min_idx = i;
				best = p0;
			}
		}
	}

	return min_idx;
}

OjaLine OjaLineSet::best_at(const OjaPoint& p)
{
    errif(size()==0,"OjaLineSet::best: set is empty");

	Point z=p.location();
    int min_idx=-1;

#ifdef QUICK_CHECK
	double min,o;

	for(int i=0; i<size(); i++)
	{
		Point d;
		d=(*this)[i].line().dir();
		d.normalize();
		o=data->oja(z + QUICK_CHECK_DELTA * d);
		if(i == 0 || o < min)
		{
			min_idx = i;
			min = o;
		}
		o=data->oja(z - QUICK_CHECK_DELTA * d);
		if(i == 0 || o < min)
		{
			min_idx = i;
			min = o;
		}
	}
#define CHECK_ALL 
#endif
	
#ifdef CHECK_ALL
	if(min_idx == -1)
	{
		OjaPoint dummy1(data);
		double dummy2;
		min_idx = scan_all_routes(dummy1,dummy2);
	}
#endif

#ifdef DERIVATIVE
    double Dpos,Dneg,Dmin=0.0;
    for(int i=0; i<size(); i++)
    {
		data->derivative(z,(*this)[i].line(),Dpos,Dneg);
		if(Dpos < Dmin)
		{
			min_idx = i;
			Dmin = Dpos;
		}
		else if(Dneg < Dmin)
		{
			min_idx = i;
			Dmin = Dneg;
		}
    }
#endif
    
#ifdef GRADIENT
	Point gr;
	gr = -(data->gradient(z));
	double min_angle,a;
	for(int i=0; i<size(); i++)
	{
		a = (*this)[i].line().angle(gr);
		if(i==0 || a < min_angle)
		{
			min_idx = i;
			min_angle = a;
		}
	}
#endif
    
    OjaLine zero(*data);
    
    if(min_idx >= 0)
		return (*this)[min_idx];
    
    return zero;
}

