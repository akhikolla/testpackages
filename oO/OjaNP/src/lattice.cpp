/* 
 *  lattice.cpp : Lattice structure for approximation.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: lattice.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include "lattice.h"
#include "error.h"

using namespace std; //df

//
// Class: Lattice
// ==============

Lattice::Lattice(const Lattice& L)
{
	nodedata=0;
	child=0;
	parent=0;
	operator=(L);
}

Lattice& Lattice::operator=(const Lattice& L)
{
	errif(parent,"Lattice::operator=: sublattice substitution unimplemented");
	
	min=L.min;
	max=L.max;
	k=L.k;
	submin=L.submin;
	submax=L.submax;

	if(nodedata)
		delete[] nodedata;	
	nodedata=new NodeData[L.nodes()];
	errif(!nodedata,"Lattice::operator=: out of memory");
	
	for(int i=0; i<L.nodes(); i++)
		nodedata[i]=L.nodedata[i];

	if(child)
		delete child;
	
	if(L.has_sub())
	{
		child=new Lattice(*L.child);
		child->parent=this;
	}
	else
		child=0;

	return *this;
}

void Lattice::forget_old_levels()
{
	errif(!has_sub(),"Lattice::forget_old_levels: no sublattices");
	errif(parent,"Lattice::forget_old_levels: unimplemented operation");
	errif(child->has_sub(),"Lattice::forget_old_levels: unimplemented operation");

	// Irrotetaan alempi kerros pois
	Lattice* oldsub;
	oldsub=child;
	child->parent=0;
	child=0;

	// Kopioidaan ja tuhotaan vanha
	*this=*oldsub;

	delete oldsub;
}

void Lattice::initialize(Point mn,Point mx,SimpleIndex I)
{
	errif(mn.dim()==0 || mx.dim()==0,"Lattice::Lattice: zero dimension");
	errif(mn.dim() != mx.dim(),"Lattice::Lattice: dimensions " << mn.dim()
	  << " and " << mx.dim() << " do not match");
	errif(I.dim() != mn.dim(),"Lattice::Lattice: illegal index dimension"
	  << I.dim());	

	for(int i=0; i<I.dim(); i++)
	{
		errif(mn[i] > mx[i],"Lattice::Lattice: bad minimum..maximum "
		  << mn << ".." << mx);
	}
	
	min=mn;
	max=mx;
	k=I;
	child=0;
	parent=0;
	nodedata=new NodeData[nodes()];
	errif(!nodedata,"Lattice::Lattice: out of memory");
}

Lattice::Lattice(Point mn,Point mx,SimpleIndex I)
{
	initialize(mn,mx,I);
}

Lattice::Lattice(Point mn,Point mx,double h)
{
	errif(mn.dim() != mx.dim(),"Lattice::Lattice: dimensions " << mn.dim()
	  << " and " << mx.dim() << " do not match");
	errif(h <= 0.0,"Lattice::Lattice: invalid length " << h);

	int k0,kmx=1;
	for(int i=0; i<mn.dim(); i++)
	{
		errif(mn[i] > mx[i],"Lattice::Lattice: bad minimum..maximum "
		  << mn << ".." << mx);

		k0=int((mx[i]-mn[i])/h);
		if(k0 > kmx)
			kmx=k0;
	}

	SimpleIndex I(mn.dim(),kmx);
	
	for(int i=0; i<mn.dim(); i++)
		I[i]=int((mx[i]-mn[i])/h);

	initialize(mn,mx,I);	
}

Lattice::~Lattice()
{
	if(child)
		delete child;
	
	delete[] nodedata;
}

int Lattice::nodes() const
{
	int sz=1;
	for(int i=0; i<dim(); i++)
		sz *= (k[i]+1);

	return sz;
}

int Lattice::max_k() const
{
	int mx=0;

	for(int i=0; i<dim(); i++)
		if(k[i] > mx)
			mx=k[i];

	return mx;
}

double Lattice::step(int i) const
{
	errif(i < 0 || i >= dim(),"Lattice::coord: illegal index " << i);

	return k[i] ? (max.coord(i)-min.coord(i))/double(k[i]) : 0.0;
}

double Lattice::max_step() const
{
	double max_=0.0;
	
	for(int i=0; i<dim(); i++)
		if(step(i) > max_)
			max_=step(i);

	return max_;
}

int Lattice::steps(int i) const
{
	errif(i < 0 || i >= dim(),"Lattice::steps: illegal index " << i);

	return k[i]+1;
}

double Lattice::box_diameter() const
{
	double sum=0.0;
	for(int i=0; i<dim(); i++)
		sum += step(i)*step(i);
	
	return sqrt(sum);
}

double Lattice::box_average_edge_length() const
{
	double sum=0.0;
	for(int i=0; i<dim(); i++)
		sum += step(i);
	
	return sum / dim();
}

double Lattice::coord(const SimpleIndex& I,int idx) const
{
	errif(idx < 0 || idx >= dim(),"Lattice::coord: illegal index " << idx);

	return min.coord(idx)+step(idx)*I[idx];
}

Point Lattice::point(const SimpleIndex& I) const
{
	errif(I.dim() != dim(),"Lattice::point: dimensions " << I.dim()
	  << " and " << dim() << " do not match");
	
	Point p(dim());

	for(int i=0; i<dim(); i++)
		p[i]=coord(I,i);

	return p;
}

NodeData& Lattice::node(const SimpleIndex& I) const
{
	errif(I.dim() != dim(),"Lattice::node: wrong dimension " << I.dim());
		
 	// indeksi = X1 + X1*Sx1 + X2*Sx1*Sx2 + X3*Sx1*Sx2*Sx3 + ...
	int idx=0,sz=1;

	for(int i=0; i<dim(); i++)
	{
		errif(I[i] < 0 || I[i] > k[i],"Lattice::node: illegal index " << I);
		idx += sz*I[i];
		sz *= (k[i]+1);
	}

	errif(idx >= nodes(),"Lattice::node: index overflow " << idx << " >= " << nodes()
	  << " with " << *this);
	
	return nodedata[idx];
}

SimpleIndex Lattice::smallest_goodness() const
{
	LatticeIterator L(*this);
	double g,best=L.node().goodness;
	SimpleIndex ret=L.index();	

	L++;
	while(L)
	{
		g=L.node().goodness;
		if(g < best)
		{
			best=g;
			ret=L.index();
		}
		L++;
	}

	return ret;
}

Lattice& Lattice::sub() const
{
	errif(!child,"Lattice::sub: no sub-lattice");

	return *child;
}

Lattice* Lattice::top_layer()
{
	Lattice* l;

	for(l=this; l->parent; l=l->parent);

	return l;
}

Lattice* Lattice::bottom_layer()
{
	Lattice* l;

	for(l=this; l->child; l=l->child);

	return l;
}

bool Lattice::in_sub(const SimpleIndex& I) const
{
	errif(I.dim() != dim(),"Lattice::point: dimensions " << I.dim()
	  << " and " << dim() << " do not match");

	if(!has_sub())
		return false;
	
	for(int i=0; i<dim(); i++)
		if(I[i] < submin[i] || I[i] > submax[i])
			return false;

	return true;
}

bool Lattice::in_parent(const SimpleIndex& I) const
{
	errif(I.dim() != dim(),"Lattice::point: dimensions " << I.dim()
	  << " and " << dim() << " do not match");

	if(!has_parent())
		return false;

	for(int i=0; i<dim(); i++)
		if(I[i] % 2)
			return false;

	return true;
}

bool Lattice::in_lattice(const SimpleIndex& I) const
{
	errif(I.dim() != dim(),"Lattice::point: dimensions " << I.dim()
	  << " and " << dim() << " do not match");
	
	for(int i=0; i<dim(); i++)
		if(I[i] < 0 || I[i] > k[i])
			return false;
	
	return true;
}

SimpleIndex Lattice::parent_node_index(const SimpleIndex& I) const
{
	errif(I.dim() != dim(),"Lattice::point: dimensions " << I.dim()
	  << " and " << dim() << " do not match");
	errif(!has_parent(),"SimpleIndex::parent_node_index: no parents");
	errif(!in_parent(I),"SimpleIndex::parent_node_index: " << I << " not in parent");

	SimpleIndex J;

	J=parent->submin;
	for(int i=0; i<I.dim(); i++)
		J[i]+=I[i] / 2;

	return J;
}

bool Lattice::focus_on(const SimpleIndex& mn0,const SimpleIndex& mx0,bool expand)
{
	errif(child,"Lattice::focus_on: sub-lattice already exist");
	errif(dim() != mn0.dim(),"Lattice::focus_on: dimensions do not match");
	errif(dim() != mx0.dim(),"Lattice::focus_on: dimensions do not match");

	for(int j=0; j<dim(); j++)
		errif(mn0[j] > mx0[j] || mn0[j] > k[j] || mx0[j] > k[j]
		  || mn0[j] < 0 || mx0[j] < 0,"Lattice::focus_on: invalid"
		  " index pair " << mn0 << " and "<< mx0 << " with " << *this);

	// Luodaan kopiot mn0:sta ja mx0:sta, jotta niiden yl�- ja alarajat
	// eiv�t sotke laajennusta
	SimpleIndex mn(dim(),-1,max_k()+1),mx(dim(),-1,max_k()+1);
	mn.get_values(mn0);
	mx.get_values(mx0);
	
	if(expand)
	{
 		mn.dec_all();
 		mx.inc_all();
		for(int i=0; i<mx.dim(); i++)
		{
			if(mx[i] > k[i])
				mx[i]=k[i];
			if(mn[i] < 0)
				mn[i]=0;
		}
	}
	
	int n=0;
	for(int j=0; j<dim(); j++)
		if((mx[j]-mn[j])*2 > n)
			n=(mx[j]-mn[j])*2;

	SimpleIndex k0(dim(),n);
	for(int j=0; j<dim(); j++)
		k0[j]=(mx[j]-mn[j])*2;

	child=new Lattice(point(mn),point(mx),k0);
	errif(!child,"Lattice::focus_on: out of memory");

	child->parent=this;
	submin=mn;
	submax=mx;
	
	return true;
}
	
void Lattice::get_bounds_over(double limit,SimpleIndex& min_,SimpleIndex& max_)
{
	// M��ritell��n aluksi min ja max vastakkaisiin kulmiin
	SimpleIndex mn(k);
	SimpleIndex mx(dim(),max_k());

	// K�yd��n l�pi solmut ja etsit��n "minimi"- ja "maksimi"-
	// solmujen indeksit
	bool found=false;
	for(LatticeLevelIterator i(*this); i; i++)
	{
		if(i.node().goodness > limit)
		{
			found=true;
 			for(int j=0; j<dim(); j++)
 			{
				if(i[j] < mn[j])
					mn[j]=i[j];
				if(i[j] > mx[j])
					mx[j]=i[j];
			}
		}
	}


	if(found)
	{
		min_=mn;
		max_=mx;
	}
	else
	{
		min_=max_=SimpleIndex();
	}
}

void Lattice::get_bounds_under(double limit,SimpleIndex& min_,SimpleIndex& max_)
{
	// M��ritell��n aluksi min ja max vastakkaisiin kulmiin
	SimpleIndex mn(k);
	SimpleIndex mx(dim(),max_k());

	// K�yd��n l�pi solmut ja etsit��n "minimi"- ja "maksimi"-
	// solmujen indeksit
	bool found=false;
	for(LatticeLevelIterator i(*this); i; i++)
	{
		if(i.node().goodness < limit)
		{
			found=true;
 			for(int j=0; j<dim(); j++)
 			{
				if(i[j] < mn[j])
					mn[j]=i[j];
				if(i[j] > mx[j])
					mx[j]=i[j];
			}
		}
	}


	if(found)
	{
		min_=mn;
		max_=mx;
	}
	else
	{
		min_=max_=SimpleIndex();
	}
}

bool Lattice::focus_over(double limit,bool expand)
{
	SimpleIndex mn,mx;

	get_bounds_over(limit,mn,mx);
	  
	if(mn.dim()==0)
		return false;

	return focus_on(mn,mx,expand);
}

bool Lattice::focus_under(double limit,bool expand)
{
	SimpleIndex mn,mx;

	get_bounds_under(limit,mn,mx);
	  
	if(mn.dim()==0)
		return false;

	return focus_on(mn,mx,expand);
}

double Lattice::volume_diameter_under(double limit)
{
	SimpleIndex mn,mx,sz;
	double sum=0.0;

	get_bounds_under(limit,mn,mx);
	if(mn.dim()==0)
		return box_diameter();
	
	sz=mx;
	sz-=mn;
		
	for(int i=0; i<dim(); i++)
		sum += double(sz[i]+1)*step(i)*double(sz[i]+1)*step(i);
	
	return sqrt(sum);
}

double Lattice::average_edge_length_under(double limit)
{
	SimpleIndex mn,mx,sz;
	double sum=0.0;

	get_bounds_under(limit,mn,mx);
	if(mn.dim()==0)
	{
		for(int i=0; i<dim(); i++)
			sum += step(i);
	}
	else
	{
		sz=mx;
		sz-=mn;
		
		for(int i=0; i<dim(); i++)
			sum += double(sz[i]+1)*step(i);
	}
	
	return sum/double(dim());
}

double Lattice::max_edge_length_under(double limit)
{
	SimpleIndex mn,mx,sz;
	double max_=0.0;

	get_bounds_under(limit,mn,mx);
	if(mn.dim()==0)
	{
		for(int i=0; i<dim(); i++)
			if(step(i) > max_)
				max_=step(i);
	}
	else
	{
		sz=mx;
		sz-=mn;
		
		for(int i=0; i<dim(); i++)
			if(double(sz[i]+1)*step(i) > max_)
				max_=double(sz[i]+1)*step(i);
	}
	
	return max_;
}

void Lattice::update_from_parent()
{
	errif(!has_parent(),"SimpleIndex::update_from_parent: lattice has no parents");

	LatticeLevelIterator L(*this);

	// Nollataan aluksi kaikki hilapisteet.
	while(L)
	{
		L.node().goodness=0.0;
		L.node().gradient=Point(dim());
		L++;
	}

	// K�yd��n kerroksen kaikki pisteet l�pi. Jokaiselle "parilliselle" pisteelle
	// k�yd��n kaikki sen naapurit l�pi ja p�ivitet��n niiden gradienttia. 'goodness'
	// arvoa k�ytet��n keskiarvon laskemiseen siten, ett� se pit�� kirjaa kuhunkin
	// pisteeseen summatuista naapureiden m��r�st�.

	SimpleIndex I; // T�ss� s�ilytett��n alkuper�isen noden arvoa
	SimpleIndex N(dim(),-1,max_k()+1); // T�h�n lasketaan k�sittelyvuorossa oleva node
	SimpleIndex A(dim(),-1,1); // Naapurustoiteraattori k�y l�pi naapurioffsetit

 	for(L=0; L; L++)
 	{
 		if(in_parent(L.index()))
 		{
 			I=parent_node_index(L.index());
			for(A.set_to_min(); A; A++)
			{
				N.fill(0);
				N+=L.index();
				N+=A;
				if(in_lattice(N))
				{
					node(N).gradient+=parent->node(I).gradient;
					node(N).goodness++;
				}
			}
 		}
 	}

	// Lasketaan lopuksi keskiarvot
	for(L=0; L; L++)
		L.node().gradient=(1.0 / L.node().goodness) * L.node().gradient;

	// Huom! goodness arvo kertoo kuinka monesta naapurista
	// gradientti on keskiarvotettu.
}

ostream& operator<<(ostream& os,const Lattice& L)
{
	os << "(" << L.minimum() << ")..(" << L.maximum() << ") / ";
	for(int i=0; i<L.dim(); i++)
	{
		if(i)
			os << 'x';
		os << (L.steps(i));
	}
	
	if(L.has_sub())
		os << endl << "  Sub: " << L.sub();

	return os;
}

//
// Class: FreeLattice
// ==================

void FreeLattice::initialize_member_list()
{
	SimpleIndex I(k.dim(),k.max_digit());
	I.set_to_min();
	bool ok=true;

	member.clear();
	
	while(ok)
	{
		int i;
		member.push_back(I);
		for(i=dim()-1; i>=0; i--)
		{
			if(I[i] < k[i])
			{
				I[i]++;
				break;
			}
			else
				I[i]=0;
		}
		if(i==-1)
			ok=false;
	}
}

void FreeLattice::remove_node(const SimpleIndex& I)
{
	errif(I.dim() != dim(),"FreeLattice::remove_node: dimensions do not match");
	member.remove(I);
}

void FreeLattice::remove_nodes_over(double limit)
{
	for(LatticeLevelIterator I(*this); I; I++)
	{
		if(I.node().goodness > limit)
			remove_node(I.index());
	}
}

void FreeLattice::remove_nodes_under(double limit)
{
	for(LatticeLevelIterator I(*this); I; I++)
	{
		if(I.node().goodness < limit)
			remove_node(I.index());
	}
}

SimpleIndex FreeLattice::smallest_goodness() const
{
	errif(points()==0,"FreeLattice::smallest_goodness: lattice is empty");
	
	FreeLatticeIterator L(*this);
	double g,best=L.node().goodness;
	SimpleIndex ret=L.index();	

	L++;
	while(L)
	{
		g=L.node().goodness;
		if(g < best)
		{
			best=g;
			ret=L.index();
		}
		L++;
	}

	return ret;
}

bool FreeLattice::focus_on(const SimpleIndex& mn0,const SimpleIndex& mx0,bool expand)
{
	bool status;
	
	status=Lattice::focus_on(mn0,mx0,expand);
	
	sub().update_from_parent();
	forget_old_levels();
	initialize_member_list();

	return status;
}

//
// Class: LatticeLevelIterator
// ===========================

LatticeLevelIterator::LatticeLevelIterator(const Lattice& l)
{
	int max;

	max=l.max_k() ? l.max_k() : 1;
	
	L=&l;
	iter=SimpleIndex(l.dim(),max);
	overflow=false;
}

LatticeLevelIterator::LatticeLevelIterator(const LatticeLevelIterator& I)
{
	L=I.L;
	overflow=I.overflow;
	iter=I.iter;
}

void LatticeLevelIterator::increment()
{
	errif(overflow,"LatticeLevelIterator::increment: overflow");

	for(int i=L->dim()-1; i>=0; i--)
		if(iter[i] < L->k[i])
		{
			iter[i]++;
			return;
		}
		else
			iter[i]=0;

	overflow=true;
}

void LatticeLevelIterator::operator++(int)
{
	increment();
}

void LatticeLevelIterator::operator=(int i)
{
	errif(i,"LatticeLevelIterator::operator=: only allowed value is 0");
	
	overflow=false;
	iter.set_to_min();
}

int& LatticeLevelIterator::operator[](int i) const
{
	errif(i < 0 || i >= dim(),"LatticeLevelIterator::operator[]: illegal index "
	  << i);
	
	return iter[i];
}

ostream& operator<<(ostream& os,const LatticeLevelIterator& I)
{
	os << "[";
	for(int i=0; i<I.dim(); i++)
	{
		if(i)
			os << ' ';
		os << I[i];
	}
	os << "]";

	return os;
}

//
// Class: FreeLatticeIterator
// ==========================

void FreeLatticeIterator::operator=(int i)
{
	errif(i,"FreeLatticeIterator::operator=: only allowed value is 0");

	if(((FreeLattice*)L)->member.empty())
	{
		overflow=true;
		iter.set_to_min();
		next=((FreeLattice*)L)->member.end();
	}
	else
	{
		overflow=false;
		next=((FreeLattice*)L)->member.begin();
		iter=*next;
	}
}

void FreeLatticeIterator::operator++(int)
{
	next++;
	if(next!=((FreeLattice*)L)->member.end())
		iter=*next;
	else
		overflow=true;
}

//
// Class: LatticeBorderIterator
// ============================

// BUG: n�ille olisi tehokkaampikin toteutus
void LatticeBorderIterator::operator++(int)
{
	errif(overflow,"LatticeBorderIterator::operator++: overflow");

	increment();
	while(L->in_sub(iter) && !overflow)
		increment();
}

void LatticeBorderIterator::operator=(int i)
{
	errif(i,"LatticeBorderIterator::operator=: only allowed value is 0");
	
	overflow=false;
	iter.set_to_min();
	while(L->in_sub(iter) && !overflow)
		increment();
}

//
// Class: LatticeIterator
// ======================

void LatticeIterator::operator++(int)
{
	errif(overflow,"LatticeIterator::operator++: overflow");

	increment();
	while(L->in_sub(iter) && !overflow)
		increment();

	if(overflow)
	{
		if(L->has_sub())
			*this=LatticeIterator(L->sub());
	}
}

void LatticeIterator::operator=(int i)
{
	errif(i,"LatticeIterator::operator=: only allowed value is 0");
	
	overflow=false;
	iter.set_to_min();
	while(L->in_sub(iter) && !overflow)
		increment();

	if(overflow)
		if(L->has_sub())
			*this=LatticeIterator(L->sub());
}
