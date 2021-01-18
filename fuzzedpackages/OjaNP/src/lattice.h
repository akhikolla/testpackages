/* $Id: lattice.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef LATTICE_H
#define LATTICE_H
#include "point.h"
#include "index.h"
#include <list>

using namespace std;

// Lattices
// ========

struct NodeData
{
	Point gradient;
	double goodness;
};

class Lattice
{
  protected:
	
	Point min,max;
	SimpleIndex k; // Suurimmat indeksit kussakin dimensiossa
	Lattice *child,*parent;
	SimpleIndex submin,submax;
	NodeData* nodedata;

	void initialize(Point mn,Point mx,SimpleIndex I);

	friend class LatticeLevelIterator;
	friend class LatticeIterator;

	int nodes() const;
	
  public:

	Lattice(const Lattice& L);
	Lattice(Point mn,Point mx,SimpleIndex I);
	Lattice(Point mn,Point mx,double h);
	virtual ~Lattice();

	Lattice& operator=(const Lattice&);
	
	int dim() const
		{return min.dim();}
	virtual int points() const
		{return nodes();}
	
	int max_k() const;
	double step(int i) const;
	double max_step() const;
	int steps(int i) const;
	double box_diameter() const;
	double box_average_edge_length() const;
	
	Point minimum() const
		{return min;}
	Point maximum() const
		{return max;}
	bool has_sub() const
		{return child != 0;}
	bool has_parent() const
		{return parent != 0;}
	bool in_sub(const SimpleIndex& I) const;
	bool in_parent(const SimpleIndex& I) const;
	bool in_lattice(const SimpleIndex& I) const;
	SimpleIndex parent_node_index(const SimpleIndex& I) const;

	Lattice& sub() const;
	Lattice* top_layer();
	Lattice* bottom_layer();
	
	NodeData& node(const SimpleIndex&) const;
	double coord(const SimpleIndex&,int) const;
	Point point(const SimpleIndex&) const;
	virtual SimpleIndex smallest_goodness() const;

	void get_bounds_under(double limit,SimpleIndex& min,SimpleIndex& max);
	void get_bounds_over(double limit,SimpleIndex& min,SimpleIndex& max);
	
	virtual bool focus_on(const SimpleIndex& mn,const SimpleIndex& mx,bool expand=false);
	bool focus_over(double limit,bool expand=false);
	bool focus_under(double limit,bool expand=false);

	void forget_old_levels();
	
	double volume_diameter_under(double limit);
	double average_edge_length_under(double limit);
	double max_edge_length_under(double limit);
	
	void update_from_parent();
};

class FreeLattice:public Lattice
{
	list<SimpleIndex> member;
	
	friend class FreeLatticeIterator;
	
	void initialize_member_list();
	
  public:
	
	FreeLattice(const FreeLattice& L) : Lattice(L)
		{initialize_member_list(); member=L.member;}
	FreeLattice(Point mn,Point mx,SimpleIndex I) : Lattice(mn,mx,I)
		{initialize_member_list();}
	FreeLattice(Point mn,Point mx,double h) : Lattice(mn,mx,h)
		{initialize_member_list();}
	
	FreeLattice& operator=(const FreeLattice&L)
		{Lattice::operator=(L); member=L.member; return *this;}

	void remove_node(const SimpleIndex& I);
	void remove_nodes_over(double limit);
	void remove_nodes_under(double limit);
	
	int points() const
		{return member.size();}
	SimpleIndex smallest_goodness() const;
	SimpleIndex bigest_goodness() const;
	
	bool focus_on(const SimpleIndex& mn,const SimpleIndex& mx,bool expand=false);
};

// Lattice iterators
// =================

class LatticeLevelIterator
{
  protected:
	
	const Lattice* L;
	SimpleIndex iter;
	bool overflow;

	void increment();
	
  public:
	
	LatticeLevelIterator(const Lattice& l);
	LatticeLevelIterator(const LatticeLevelIterator& I);
	virtual ~LatticeLevelIterator() {}
	
	int dim() const
		{return L->dim();}
	
	virtual void operator=(int);
	virtual void operator++(int);
	
	int& operator[](int index) const;
	operator void*() const
		{return overflow ? (void*) 0 : (void *) -1;}
	
	SimpleIndex index()
		{return iter;}
	Point point() const
		{return L->point(iter);}
	double coord(int i) const
		{return L->coord(iter,i);}
	NodeData& node() const
		{return L->node(iter);}
};

class FreeLatticeIterator:public LatticeLevelIterator
{
	list<SimpleIndex>::iterator next;
	
  public:
	
	FreeLatticeIterator(const FreeLattice& L)
		: LatticeLevelIterator(L) {*this=0;}
	FreeLatticeIterator(const FreeLatticeIterator& I)
		: LatticeLevelIterator(I) {*this=0;}
	virtual ~FreeLatticeIterator() {}

	virtual void operator=(int);
	virtual void operator++(int);
};

class LatticeBorderIterator:public LatticeLevelIterator
{	
  public:

	LatticeBorderIterator(const Lattice& L)
		: LatticeLevelIterator(L) {*this=0;}
	LatticeBorderIterator(const LatticeBorderIterator& I)
		: LatticeLevelIterator(I) {*this=0;}
	virtual ~LatticeBorderIterator() {}
	
	virtual void operator=(int);
	virtual void operator++(int);
};

class LatticeIterator:public LatticeBorderIterator
{
  public:

	LatticeIterator(const Lattice& L0)
		: LatticeBorderIterator(L0) {*this=0;}
	LatticeIterator(const LatticeIterator& I)
		: LatticeBorderIterator(I) {}
	virtual ~LatticeIterator() {}

	virtual void operator=(int); 
	virtual void operator++(int);
};

ostream& operator<<(ostream&,const Lattice& L);
ostream& operator<<(ostream&,const LatticeLevelIterator& L);
#endif
