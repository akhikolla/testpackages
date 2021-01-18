/* $Id: index.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef INDEX_H
#define INDEX_H
#include <iostream>
#include <list>
#include <vector>
#include <set>

using namespace std;

class SimpleIndex
{
    bool overflow;
    int digits;
    int *digit;
    int min;
    int max;
    
    void initialize(int dim,int min_value,int max_value);
    
	friend int compare(const SimpleIndex&,const SimpleIndex&);

  public:
	
    SimpleIndex();
    SimpleIndex(int dim,int max_value)
		{initialize(dim,0,max_value-1);}
    SimpleIndex(int dim,int min_value,int max_value)
		{initialize(dim,min_value,max_value);}
    SimpleIndex(const SimpleIndex& i);
    virtual ~SimpleIndex();
    
    SimpleIndex& operator=(const SimpleIndex& i);
    SimpleIndex& operator++(int);
    SimpleIndex& operator--(int);
	SimpleIndex& operator+=(const SimpleIndex& i);
	SimpleIndex& operator-=(const SimpleIndex& i);
	
    int& operator[](int index) const;
    operator void*() const
		{return overflow ? (void*) 0 : (void *) -1;}
	
	int dim() const
		{return digits;}
	int min_digit() const
		{return min;}
	int max_digit() const
		{return max;}

	void fill(int n);
    void set_to_max()
		{fill(max);}
    void set_to_min()
		{fill(min);}
	void get_values(const SimpleIndex& I);

	void add_all(int how_much);
	void dec_all()
		{add_all(-1);}
	void inc_all()
		{add_all(+1);}

    void random(); 
};

int compare(const SimpleIndex&,const SimpleIndex&);
inline bool operator==(const SimpleIndex& I1,const SimpleIndex& I2)
{
	return compare(I1,I2)==0;
}
inline bool operator!=(const SimpleIndex& I1,const SimpleIndex& I2)
{
	return compare(I1,I2)!=0;
}
inline bool operator<(const SimpleIndex& I1,const SimpleIndex& I2)
{
	return compare(I1,I2)<0;
}
inline bool operator>(const SimpleIndex& I1,const SimpleIndex& I2)
{
	return compare(I1,I2)>0;
}
inline bool operator<=(const SimpleIndex& I1,const SimpleIndex& I2)
{
	return compare(I1,I2)<=0;
}
inline bool operator>=(const SimpleIndex& I1,const SimpleIndex& I2)
{
	return compare(I1,I2)>=0;
}

class Index
{
    bool overflow;
    int digits;
    int *digit;
    int max;
    
    bool inc_digit(int index);
    bool dec_digit(int index);

	friend int compare(const Index&,const Index&);
	
    void initialize(int dim,int max_value);
    
  public:
	
    Index();
    Index(int dim,int max_value)
		{initialize(dim,max_value);}
	Index(int dim, int max_value, set<int> digits);
    Index(const Index& i);
    virtual ~Index();
    
    Index& operator=(const Index& i);
    Index& operator++(int);
    Index& operator--(int);
    int& operator[](int index) const; 
    bool validate();
	bool is_valid() const;
	
    operator void*() const
		{return overflow ? (void*) 0 : (void *) -1;}
	
	int dim() const
		{return digits;}
	int limit() const
		{return max;}
	void set_limit(int lim)
		{ max = lim;}
	virtual int combinations() const;
	int has(int value) const; 
	bool has_sub_set(const Index& I) const;
	Index intersection(const Index& I) const;
	
	void random_with(int index);
    void random_with(int ind1,int ind2); 
    void random();
    void set_to_max();
    void set_to_min();
};

int ord(const Index& I);

int compare(const Index&,const Index&);
inline bool operator==(const Index& I1,const Index& I2)
{
	return compare(I1,I2)==0;
}
inline bool operator!=(const Index& I1,const Index& I2)
{
	return compare(I1,I2)!=0;
}
inline bool operator<(const Index& I1,const Index& I2)
{
	return compare(I1,I2)<0;
}
inline bool operator>(const Index& I1,const Index& I2)
{
	return compare(I1,I2)>0;
}
inline bool operator<=(const Index& I1,const Index& I2)
{
	return compare(I1,I2)<=0;
}
inline bool operator>=(const Index& I1,const Index& I2)
{
	return compare(I1,I2)>=0;
}
	

class IndexSet
{
  protected:
    
    bool overflow;
    int digits;
    Index *digit;
    Index *max;
	
    bool inc_digit(int idx);
	
    void initialize(int num_of_idxs,int dim,int max_value);
	bool search_for_common_digit(int* I) const;
	
	friend int compare(const IndexSet&,const IndexSet&);
	
  public:
	
    IndexSet();
    IndexSet(int num_of_idxs,int dim,int max_value)
		{initialize(num_of_idxs,dim,max_value);}
	IndexSet(int dim, int max_value, const vector<set<int> >&  indexes);
    IndexSet(const IndexSet& i);
    virtual ~IndexSet();
    IndexSet& operator=(const IndexSet& i);    
   
    virtual IndexSet& operator++(int);
    Index& operator[](int index) const;
    operator void*() const
		{return overflow ? (void*) 0 : (void *) -1;}
    int indices() const
		{return digits;}
    int dim() const
		{return indices() ? digit[0].dim() : 0;} 
    int limit() const
		{return indices() ? digit[0].limit() : 0;} 
    virtual int combinations() const; 

	int has(Index value) const;
	int has(int value) const;
	int common_digit() const; 
    bool has_common_digit() const
		{return common_digit() >= 0;} 
    bool has_two_common_digits() const; 
    bool has_two_common_digits(int& n1,int& n2) const; 
	int how_many_common_digits() const; 
	
    bool validate(); 
    void random(); 

	Index common_part() const;
	Index longest_common_subset(int& how_many) const;
	IndexSet sub_set_without(int number) const;
	IndexSet sub_set_with(int number) const;
	IndexSet sub_set(const Index& components) const;
	IndexSet complement_sub_set(const Index& components) const;
};

int compare(const IndexSet&,const IndexSet&);
inline bool operator==(const IndexSet& I1,const IndexSet& I2)
{
	return compare(I1,I2)==0;
}
inline bool operator!=(const IndexSet& I1,const IndexSet& I2)
{
	return compare(I1,I2)!=0;
}
inline bool operator<(const IndexSet& I1,const IndexSet& I2)
{
	return compare(I1,I2)<0;
}
inline bool operator>(const IndexSet& I1,const IndexSet& I2)
{
	return compare(I1,I2)>0;
}
inline bool operator<=(const IndexSet& I1,const IndexSet& I2)
{
	return compare(I1,I2)<=0;
}
inline bool operator>=(const IndexSet& I1,const IndexSet& I2)
{
	return compare(I1,I2)>=0;
}

ostream& operator <<(ostream& os,const SimpleIndex& I);
ostream& operator <<(ostream& os,const Index& I);
ostream& operator <<(ostream& os,const IndexSet& I);

class IndexType;
// BUG: Tehoton toteutus. Pit�isi olla Index** komponentteja sek�
// suorempi algoritmi.
class IndexIdentifier
{
  protected:

	int dimension; // Avaruuden dimensio
	int max_parts; // Alkuper�isen taulukon koko
	int parts; // Taulukon nykyinen koko
	Index* part;

	friend ostream& operator<<(ostream&,const IndexIdentifier&);
	friend int compare(const IndexIdentifier&,const IndexIdentifier&);

	
	void sort();
	void add_to_name(const Index& I);
	void delete_from_name(const Index& I);
	Index intersection_of(const Index& K);
	bool apply_k_intersect(int k);
	
  public:

 	IndexIdentifier();
 	IndexIdentifier(int dim);
 	virtual ~IndexIdentifier();
 	IndexIdentifier(const IndexIdentifier& Id);
	
	IndexIdentifier& operator=(const IndexIdentifier& Id);
	
	virtual void get(const IndexSet& I);
	void clear()
		{parts=0;}
	
	int dim() const; // Objektin dimensio
	int space_dim() const // Avaruuden
		{return dimension;}
	int limit() const
		{return parts ? part[0].limit() : 0;}
	int is_single_point() const;

    // Dimensioltaan yht� suurempien objektien m��r� (tai INT_MAX jos ylivuoto).
	int sup_objects() const;
	// Etsii kaikki ylemm�n dimension objektit.
	void put_sup_objects(set<IndexIdentifier>& S,int dim);	
	
	bool operator==(const IndexIdentifier& Id) const
	{return compare(*this,Id)==0;}
	bool operator!=(const IndexIdentifier& Id) const
	{return !(*this==Id);}
	bool operator<(const IndexIdentifier& Id) const
	{return compare(*this,Id)==-1;}

	void show_format(); // osien koot d1+d2+d3 notaationa
	IndexIdentifier format() const; // osien koot d1+d2+d3 indeksin�

	vector<int> partitions() const; // Palauttaa hypertasonippujen koot.
	int real_partitions(int max) const; // number of partitions without bounds
};

istream& operator>>(istream&,const IndexIdentifier&);
ostream& operator<<(ostream&,const IndexIdentifier&);

class JokerIdentifier : public IndexIdentifier
{
	friend ostream& operator<<(ostream&,const JokerIdentifier&);

	void joker_expand();

  public:

 	JokerIdentifier(): IndexIdentifier() {}
 	JokerIdentifier(const IndexIdentifier& Id) : IndexIdentifier(Id) {/*joker_expand(); is called from get()*/}
 	virtual ~JokerIdentifier() {}

	JokerIdentifier& operator=(const JokerIdentifier& J)
	{IndexIdentifier::operator=(J); return *this;}

	void get(const IndexSet& I)
		{IndexIdentifier::get(I); joker_expand();}
	void get(const IndexIdentifier& I)
		{*this=I; joker_expand();}
	
	bool put(IndexSet& I) const;

	int jokers() const;
	int dim() const
		{return parts ? dimension - parts: -1;}
	
	void put_sup_identifiers(set<JokerIdentifier>& L,int d) const;
	void convert_to_identifiers(set<IndexIdentifier>& N,int d=-1) const;

	IndexSet get_random_sup_object(int min_dim) const;
	int sup_combinations(int min_dim) const; // BUG: trivial limit only
};

inline int compare(const JokerIdentifier&J1,const JokerIdentifier&J2)
{return compare((IndexIdentifier)J1,(IndexIdentifier)J2);}


void inline remove_subset(set<int>& from, const set<int>& rem){
	for (set<int>::iterator i = rem.begin(); i != rem.end(); i++){
		int e = *i;
		from.erase(e);
	}
}

bool inline is_subset(set<int>& bigset, const set<int>& subset){
	for (set<int>::iterator i = subset.begin(); i != subset.end(); i++){
		int e = *i;
		if (!bigset.count(e))
			return false;
	}
	return true;
}

vector<set<int> > generate_unique_indexes(vector<set<int> >& bounds_crossing_indexes, set<int>& partset, int dim);

#endif
