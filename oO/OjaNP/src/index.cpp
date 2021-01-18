/* 
 *  index.cpp : Index geometry routines.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: index.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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
#include <limits.h>
#include <algorithm>
#include <iterator>
#include "random.h"
#include "error.h"
#include "misc.h"
#include "index.h" 

using namespace std; //df

//  CLASS: SimpleIndex
//  ==================

SimpleIndex::SimpleIndex()
{
    overflow = false;
    digits = 0;
    min = 0;
    max = 0;
    digit = 0;
}

void SimpleIndex::initialize(int dim,int min_value,int max_value)
{
    errif(dim < 1,"SimpleIndex::SimpleIndex: non-positive dimension " << dim);
	errif(min_value > max_value,"SimpleIndex::SimpleIndex: illegal range "
	  << min_value << ".." << max_value);
	
    overflow = false;
    digits = dim;
    min = min_value;
    max = max_value;
    digit = new int[digits];
    
    errif(!digit,"SimpleIndex::SimpleIndex: out of memory");
    set_to_min();
}

SimpleIndex::SimpleIndex(const SimpleIndex& i)
{
    digit = 0;
    operator=(i);
}

SimpleIndex& SimpleIndex::operator=(const SimpleIndex& i)
{
    overflow = i.overflow;
    digits = i.digits;
    min = i.min;
    max = i.max;
    if(digit)
		delete[] digit;
    digit = new int[digits];
    
    errif(!digit,"SimpleIndex::operator=: out of memory");
    for(int j=0; j<digits; j++)
		digit[j] = i.digit[j];
	
    return *this;
}

void SimpleIndex::get_values(const SimpleIndex& I)
{
	errif(dim() != I.dim(),"SimpleIndex::get_values: dimensions do not match");
	overflow=false;
	
    for(int j=0; j<digits; j++)
	{
		errif(I.digit[j] < min || I.digit[j] > max,"SimpleIndex::get_values: "
		  "value " << I.digit[j] << " out of range " << min << ".." << max);
		
		digit[j] = I.digit[j];
	}
}

SimpleIndex::~SimpleIndex()
{
    if(digit)
		delete[] digit;
	digit = 0;
}

SimpleIndex& SimpleIndex::operator++(int)
{
    errif(overflow,"SimpleIndex::operator++: overflow");
    
	for(int i=digits-1; i>=0; i--)
		if(digit[i] < max)
		{
			digit[i]++;
			return *this;
		}
		else
			digit[i]=min;

	overflow=true;
	
    return *this;
}

SimpleIndex& SimpleIndex::operator--(int)
{
    errif(overflow,"SimpleIndex::operator--: underflow");
	
	for(int i=digits-1; i>=0; i--)
		if(digit[i] > min)
		{
			digit[i]--;
			return *this;
		}
		else
			digit[i]=max;

	overflow=true;
    
    return *this;
}

int& SimpleIndex::operator[](int index) const
{
    errif(!digit,"SimpleIndex::operator[]: no digits allocated");
    errif(index < 0 || index >= digits,"SimpleIndex::operator[]: illegal index " << index);
	
    return digit[index];
}

void SimpleIndex::fill(int n)
{
	errif(n < min || n > max,"SimpleIndex::fill: value " << n
	  << " out of bounds " << min << ".." << max);
	
    overflow = false;
    for(int i=0; i<digits; i++)
		digit[i] = n;
}

SimpleIndex& SimpleIndex::operator+=(const SimpleIndex& I)
{
	errif(dim() != I.dim(),"SimpleIndex::operator+=: dimensions "
	  << dim() << " and " << I.dim() << " do not match");

	overflow=false;
	for(int i=0; i<dim(); i++)
	{
		digit[i]+=I.digit[i];
		if(digit[i] < min)
		{
			digit[i] = min;
			overflow = true;
		}
		if(digit[i] > max)
		{
			digit[i] = max;
			overflow = true;
		}
	}

	return *this;
}

SimpleIndex& SimpleIndex::operator-=(const SimpleIndex& I)
{
	errif(dim() != I.dim(),"SimpleIndex::operator-=: dimensions "
	  << dim() << " and " << I.dim() << " do not match");
	
	for(int i=0; i<dim(); i++)
	{
		digit[i]-=I.digit[i];
		if(digit[i] < min)
			digit[i] = min;
		if(digit[i] > max)
			digit[i] = max;
	}

	return *this;
}

void SimpleIndex::add_all(int how_much)
{
	for(int i=0; i<dim(); i++)
	{	
		digit[i]+=how_much;
		if(digit[i] > max)
		{
			digit[i] = max;
			overflow = true;
		}
		else if(digit[i] < min)
		{
			digit[i] = min;
			overflow = true;
		}
	}
}

int compare(const SimpleIndex& I1,const SimpleIndex& I2)
{
    errif(!I1.digit || !I2.digit,"compare(SimpleIndex,SimpleIndex): no digits allocated");
    errif(I1.dim() != I2.dim(),"compare(SimpleIndex,SimpleIndex): dimensions " << I1.dim()
	  << " and " << I2.dim() << " do not match");

	int diff;
	for(int i=0; i<I1.digits; i++)
	{
		diff=I1.digit[i] - I2.digit[i];
		if(diff)
			return diff;		
	}

	return 0;
}

void SimpleIndex::random()
{
	for(int i=0; i<digits; i++)
		digit[i] = ::random(min,max);
}

ostream& operator <<(ostream& os,const SimpleIndex& I)
{
    for(int i=0; i<I.dim()-1; i++)
		os << I[i] << ' ';
    os << I[I.dim()-1];
	
    return os;
}

//  CLASS: Index
//  ============

Index::Index()
{
    overflow = false;
    digits = 0;
    max = 0;
    digit = 0;
}

Index::Index(int dim, int max_value, set<int> digits)
{
	initialize(dim, max_value);

	errif(digits.size() != dim, "Index::Index: dim of digits differs ");

	int i = 0;
	for (set<int>::iterator si = digits.begin(); si != digits.end(); si++)
	{
		digit[i++] = *si;
	}
}

void Index::initialize(int dim, int max_value)
{
    errif(dim < 1,"Index::Index: non-positive dimension " << dim);
    errif(max_value < 1,"Index::Index: non-positive max_value " << max_value);
    errif(max_value < dim,"Index::Index: max_value " << max_value << " too small");
	
    overflow = false;
    digits = dim;
    max = max_value;
    digit = new int[digits];
    
    errif(!digit,"Index::Index: out of memory");
    set_to_min();
}

Index::Index(const Index& i)
{
    digit = 0;
    operator=(i);
}

Index& Index::operator=(const Index& i)
{
    overflow = i.overflow;
    digits = i.digits;
    max = i.max;
    if(digit)
		delete[] digit;
    digit = new int[digits];
    
    errif(!digit,"Index::operator=: out of memory");
    for(int j=0; j<digits; j++)
		digit[j] = i.digit[j];
	
    return *this;
}

Index::~Index()
{
    if(digit)
		delete[] digit;
	digit = 0;
}

bool Index::inc_digit(int i)
{   
    errif(!digit,"Index::inc_digit: no digits allocated");
    
    if (digit[i] < max - digits + i)
    {
		digit[i]++;
		return true;
    }
	
    if(i==0 || !inc_digit(i-1))
		return false;
    
    digit[i] = digit[i-1] + 1;
	
    return true;
}

bool Index::dec_digit(int i)
{   
    errif(!digit,"Index::dec_digit: no digits allocated");
	
    if(digit[i]==i)
		return false;
    
    digit[i]--;
    if(i==0)
		return true;
    
    if(digit[i] == digit[i-1])
    {
		dec_digit(i-1);
		digit[i] = max - digits + i;
    }
    
    return true;
}

Index& Index::operator++(int)
{
    errif(overflow,"Index::operator++: overflow");
    
    if(!inc_digit(digits-1))
		overflow = true;
    
    return *this;
}

Index& Index::operator--(int)
{
    errif(overflow,"Index::operator--: underflow");
	
    if(!dec_digit(digits-1))
		overflow = true;
    
    return *this;
}

int Index::combinations() const
{
	return choices(limit(),dim());
}

int& Index::operator[](int index) const
{
    errif(!digit,"Index::operator[]: no digits allocated");
    errif(index < 0 || index >= digits,"Index::operator[]: illegal index " << index
	  << " with '" << *this << "'");
	
    return digit[index];
}

int Index::has(int value) const
{
	// BUG: pit�isi olla puolitushaku
    errif(!digit,"Index::has: no digits allocated");
    int sum=0;
    for(int i=0; i<digits; i++)
		if (digit[i]==value)
			sum++;
    
    return sum;
}

bool Index::has_sub_set(const Index& I) const
{
	// BUG: tehoton
	for(int i=0; i<digits; i++)
		if(!I.has(digit[i]))
			return false;

	return true;
}

Index Index::intersection(const Index& I) const
{
	if(dim()==0 || I.dim()==0)
		return Index();
	
	errif(limit() != I.limit(),"Index::intersection: limit mismatch " << limit() << "!=" << I.limit());
	
	// BUG: tehoton
	int* tmp=new int[digits];
	int sz=0;
	
	for(int i=0; i<digits; i++)
		if(I.has(digit[i]))
			tmp[sz++]=digit[i];

	if(sz)
	{
		Index J(sz,I.limit());
		for(int i=0; i<sz; i++)
			J[i]=tmp[i];

		delete[] tmp;
		return J;
	}
	else
	{
		delete[] tmp;
		return Index();
	}
}

void Index::random()
{
	do
	{
		for(int i=0; i<digits; i++)
			digit[i] = ::random(0,max-1);
	}
	while(!validate());
}

void Index::random_with(int index)
{
    errif(index < 0 || index >= max,"Index::random_with: illegal index value "
	  << index);
	
    if(dim()==limit())
    {
		set_to_min();
		return;
    }

	do
	{
		digit[0]=index;
		for(int i=1; i<digits; i++)
			digit[i] = ::random(0,max-1);
	}
	while(!validate());	
}

void Index::random_with(int ind1,int ind2)
{
    errif(ind1 < 0 || ind1 >= max,"Index::random_with: illegal index value "
	  << ind1);
    errif(ind2 < 0 || ind2 >= max,"Index::random_with: illegal index value "
	  << ind2);
    errif(ind1==ind2,"Index::random_with: both numbers are the same");
    errif(dim() < 2,"Index::random_with: unable to fit two numbers in one");
    if(dim()==limit())
    {
		set_to_min();
		return;
    }
	
	do
	{
		digit[0]=ind1;
		digit[1]=ind2;
		for(int i=2; i<digits; i++)
			digit[i] = ::random(0,max-1);
	}
	while(!validate());	
}

void Index::set_to_max()
{
    overflow = false;
    for(int i=0; i<digits; i++)
		digit[i] = max-(dim()-i);
}

void Index::set_to_min()
{
    overflow = false;
    for(int i=0; i<digits; i++)
		digit[i] = i;
}

bool Index::validate()
{
    for(int i=0; i<digits-1; i++)
		for(int j=i+1; j<digits; j++)
			if(digit[i] > digit[j])
			{
				int tmp;
				tmp=digit[i];
				digit[i]=digit[j];
				digit[j]=tmp;
			}
    for(int i=1; i<dim(); i++)
		if(digit[i-1]==digit[i])
			return false;
    
    return true;
}

bool Index::is_valid() const
{
    for(int i=1; i<digits; i++)
			if(digit[i-1] >= digit[i])
				return false;
    
    return true;
}

int compare(const Index& I1,const Index& I2)
{
    errif(!I1.digit || !I2.digit,"compare(Index,Index): no digits allocated");
    errif(I1.dim() != I2.dim(),"compare(Index,Index): dimensions " << I1.dim()
	  << " and " << I2.dim() << " do not match");

	int diff;
	for(int i=0; i<I1.digits; i++)
	{
		diff=I1.digit[i] - I2.digit[i];
		if(diff)
			return diff;		
	}

	return 0;
}

ostream& operator <<(ostream& os,const Index& I)
{
	if(I.dim()==0)
		return os;
	
    for(int i=0; i<I.dim()-1; i++)
		os << I[i] << ' ';
    os << I[I.dim()-1];
	
    return os;
}

// Palauttaa ensimm�isen (d,n) indeksin, jonka ensimm�inen numero on k
//
// Huom. \sum_{k=i}^j \over{k}{d}= \over{j+1}{d+1} - \over{i}{d+1}
// 
static int _ord(int d,int n,int k)
{
	if(k==0)
		return 0;

	// lasketaan \sum_{k=n-k}^{n-1} \over{k}{d-1}
	
	return choices(n,d) - choices(n-k,d);
}

int ord(const Index& I)
{
	int sum=_ord(I.dim(),I.limit(),I[0]);
	for(int i=1; i<I.dim(); i++)
		sum += _ord(I.dim()-i,I.limit()-I[i-1]-1,I[i]-I[i-1]-1);
	return sum;
}


//  CLASS: IndexSet
//  ===============

IndexSet::IndexSet()
{
    overflow = false;
    digits = 0;
    digit = 0;
    max = 0;
}

void IndexSet::initialize(int num_of_idxs,int dim,int max_value)
{
    errif(num_of_idxs < 1,"IndexSet::IndexSet: non-positive index set size " << num_of_idxs);
	
    overflow = false;
    digits = num_of_idxs;
    
    digit = new Index[digits];
    max = new Index[digits];
    errif(!digit || !max,"IndexSet::IndexSet: out of memory");
    
    Index Imax(dim,max_value),Imin(dim,max_value);
    Imax.set_to_max();
    for(int i=0; i<digits; i++)
    {
		digit[i] = Imin;
		max[i] = Imax;
		for(int j=i; j; j--)
			digit[i]++;
		for(int j=digits-i-1; j; j--)
			max[i]--;
    }
}

IndexSet::IndexSet(int dim, int max_value, const vector<set<int> >& indexes)
{
	initialize(indexes.size(), dim, max_value);

	for (int i = 0; i < indexes.size(); i++)
	{
		digit[i] = Index(indexes[i].size(), max_value, indexes[i]);
	}
}

IndexSet::IndexSet(const IndexSet& i)
{
    digit = 0;
    operator=(i);
}

IndexSet& IndexSet::operator=(const IndexSet& I)
{
    overflow = false;
    digits = I.digits;
    if(digit)
    {
		delete[] digit;
		delete[] max;
    }
    
    digit = new Index[digits];
    max = new Index[digits];
    
    errif(!digit || !max,"IndexSet::operator=: out of memory");
    for(int j=0; j<digits; j++)
    {
		digit[j] = I.digit[j];
		max[j] = I.max[j];
    }
	
    return *this;
}

IndexSet::~IndexSet()
{
    if(digit)
    {
		delete[] max;
		delete[] digit;
		digit = 0;
    }
}

bool IndexSet::inc_digit(int i)
{   
    errif(!digit,"IndexSet::inc_digit: no digits allocated");
    
    if (!(digit[i] == max[i]))
    {
		digit[i]++;
		return true;
    }
	
    if(i==0 || !inc_digit(i-1))
		return false;
	
    digit[i] = digit[i-1];
    digit[i]++;
	
    return true;
}

IndexSet& IndexSet::operator++(int)
{
    errif(overflow,"IndexSet::operator++: overflow");
	
    if(!inc_digit(digits-1))
		overflow = true;
    
    return *this;
}

Index& IndexSet::operator[](int index) const
{
    errif(!digit,"IndexSet::operator[]: no digits allocated");
    errif(index < 0 || index >= digits,"IndexSet::operator[]: illegal index " << index
	  << " with '" << *this << "'");
	
    return digit[index];
}

int IndexSet::has(int value) const
{
    errif(!digit,"IndexSet::has: no indices");
    int sum=0;
    for(int i=0; i<digits; i++)
		sum+=digit[i].has(value);
    
    return sum;
}

int IndexSet::has(Index value) const
{
	errif(!digit, "IndexSet::has: no indices");
	for (int i = 0; i < digits; i++)
	if (digit[i] == value)
		return true;
	return false;
}

int IndexSet::combinations() const
{
    errif(1,"IndexSet::combinations: UNIMPLEMENTED");
    return 0;
}

bool IndexSet::search_for_common_digit(int* I) const
{
	/* K�yd��n l�pi kaikki viimeisen indeksin lukuarvot.
	   Kelataan kullekkin luvulle kaikkia muita indekseit� siten, ett� l�ytyy
	   v�hint��n yht�suuri. Jos kelaus pys�htyy samaan lukuun kaikissa
	   indekseiss�, yhteinen tekij� on l�ytynyt ja palautetaan TRUE.*/
	
	/* I = Laskuri kullekkin indeksille I_1, I_2, ... , I_{digits-1}.
	   Kukin laskuri k�y l�pi arvot 0,...,{dim-1}. */

	/* Alias viimeiselle indeksilaskurille I_{digits-1} */
	int& k = I[digits-1];

	/* Alias viimeiselle indeksille */
	Index& L = digit[digits-1];

	/* M��ritelm�n mukaan mitk��n indeksit eiv�t ole yksiulotteisessa
	   tapauksessa samoja */
    if(dim()==1)
		return -1;

    while(k<dim())
    {
		/* Kelaus ja ylivuototarkistus */
		for(int i=0; i<digits-1; i++)
			while(digit[i][I[i]] < L[k])
			{
				I[i]++;
				if(I[i] >= dim())
					return false;
			}
		
		int j=0;

		/* L�ytyik� samoja */
		for(j=0; j<digits-1; j++)
			if(digit[j][I[j]] != L[k])
				break;

		k++; /* Lis�t��n indeksi� valmiiksi mahdollista seuraavaa kutsu-
				kertaa varten */
		
		if(j==digits-1)
			return true;
	}
    
    return false;
}

int IndexSet::common_digit() const
{
	errif(!digit,"IndexSet::common_digit: no indices");
	errif(digits <= 1,"IndexSet::common_digit: all digits are common in "
	  << *this);

	int* I=new int[digits];
	   
    for(int i=0; i<digits; i++)
		I[i]=0;

	if(search_for_common_digit(I))
	{
		delete[] I;
		return digit[0][I[0]];
	}
	else
	{
		delete[] I;
		return -1;
	}
}

bool IndexSet::has_two_common_digits(int& n1,int& n2) const
{
    errif(!digit,"IndexSet::has_two_common_digit: no indices");
    errif(digits <= 1,
	  "IndexSet::has_two_common_digit: all digits are common in " << *this);

	int* I=new int[digits];
	   
    for(int i=0; i<digits; i++)
		I[i]=0;

	if(search_for_common_digit(I))
	{
		n1 = digit[0][I[0]];
		if(search_for_common_digit(I))
		{
			n2 = digit[0][I[0]];
			delete[] I;
			return true;
		}
	}	

	n1 = -1;
	n2 = -1;
	
	delete[] I;
	return false;
}

bool IndexSet::has_two_common_digits() const
{
	int dummy1,dummy2;

	return has_two_common_digits(dummy1,dummy2);
}

int IndexSet::how_many_common_digits() const
{
    errif(!digit,"IndexSet::how_many_common_digits: no indices");
    errif(digits <= 1,
	  "IndexSet::how_many_common_digits: all digits are common in " << *this);

	if(dim()==1)
		return 0;
	
	int* I=new int[digits];
	   
    for(int i=0; i<digits; i++)
		I[i]=0;

	int n;
	for(n=0; search_for_common_digit(I); n++);

	delete[] I;
	return n;
}

bool IndexSet::validate()
{
	for(int i=0; i<digits; i++)
		if(!digit[i].validate())
			return false;
	
    for(int i=0; i<digits-1; i++)
		for(int j=i+1; j<digits; j++)
			if(digit[i] > digit[j])
			{
				Index tmp;
				tmp=digit[i];
				digit[i]=digit[j];
				digit[j]=tmp;
			}
	
    for(int i=1; i<digits; i++)
		if(digit[i-1]==digit[i])
			return false;
	
    return true;
}

void IndexSet::random()
{
	do
	{
		for(int i=0; i<indices(); i++)
			digit[i].random();
	}
	while(!validate());
}

int compare(const IndexSet& I1,const IndexSet& I2)
{
    errif(!I1.digit || !I2.digit,"compare(IndexSet,IndexSet): no digits allocated");
    errif(I1.dim() != I2.dim(),"compare(IndexSet,IndexSet): dimensions " << I1.dim()
	  << " and " << I2.dim() << " do not match");
    errif(I1.indices() != I2.indices(),"compare(IndexSet,IndexSet): indices " << I1
	  << " and " << I2 << " do not match");

	int diff;
	for(int i=0; i<I1.digits; i++)
		for(int j=0; j<I1.dim(); j++)
		{
			
			diff=I1.digit[i][j] - I2.digit[i][j];
			if(diff)
				return diff;		
		}

	return 0;
}

ostream& operator<<(ostream& os,const IndexSet& I)
{
    if(I.indices()==0)
    {
		os << "empty";
		return os;
    }
    
    for(int i=0; i<I.indices()-1; i++)
		os << I[i] << "  ";
    
    os << I[I.indices()-1];
	
    return os;
}

Index IndexSet::longest_common_subset(int& how_many) const // BUG: tehoton toteutus
{
	/* K�sitell��n erikoistapaukset */
	if(indices()==0 || dim()==1)
		return Index();
	if(indices()==1)
		return digit[0];
		
	int hi=0; // Suurin yhteisten numeroiden m��r� kullakin hetkell�
	int nmbhi=0; // Kuinka monta t�h�n m��r��n ylt�v�� lukua on
	int* idx=new int[dim() * indices()]; // Taulukko kaikista numeroista
	int* nmb=new int[dim() * indices()]; // Kunkin numeron esiintymiskerrat
	int next=0; // Seuraava vapaa paikka edellisiss� taulukoissa
	
	for(int i=0; i<dim(); i++)
		for(int n=0; n<indices(); n++)
		{
			int k=digit[n][i];
			/* Lis�t��n yhdell� 'k'den m��r��. */

			/* Etsit�'n ensin 'k':n sijainti 'l':��n. */
			int l=next-1;
			while(l >= 0)
			{
				if(idx[l]==k)
					break;
				else
					l--;
			}
			
			/* Ei l�ytynyt. Luodaan uusi paikka. */
			if(l < 0)
			{
				l=next;
				idx[l]=k;
				next++;
				nmb[l]=0;
			}

			/* P�ivitet��n taulukko. */
			nmb[l]++;
			if(nmb[l] > hi)
			{
				hi=nmb[l];
				nmbhi=1;
			}
			else if(nmb[l] == hi)
				nmbhi++;
		}

	Index R(nmbhi,limit());

	int j=0;
	for(int i=0; i<next; i++)
		if(nmb[i]==hi)
			R[j++]=idx[i];

	how_many=hi;

	delete[] idx;
	delete[] nmb;
	
	return R;
}

IndexSet IndexSet::sub_set_without(int number) const
{
	Index** sub=new Index*[indices()];
	int subs=0;
	
	for(int i=0; i<indices(); i++)
		if(!digit[i].has(number))
			sub[subs++]=&digit[i];

	if(!subs)
	{
		delete[] sub;
		return IndexSet();
	}
	
	IndexSet R(subs,dim(),limit());

	for(int i=0; i<subs; i++)
		R[i]=*sub[i];

	delete[] sub;

	return R;
}

IndexSet IndexSet::sub_set_with(int number) const
{
	Index** sub=new Index*[indices()];
	int subs=0;
	
	for(int i=0; i<indices(); i++)
		if(digit[i].has(number))
			sub[subs++]=&digit[i];

	if(!subs)
	{
		delete[] sub;
		return IndexSet();
	}
	
	IndexSet R(subs,dim(),limit());

	for(int i=0; i<subs; i++)
		R[i]=*sub[i];

	delete[] sub;
	return R;
}

Index IndexSet::common_part() const
{
	Index C;
	int how_many;
	
	C=longest_common_subset(how_many);
	if(how_many < indices())
		return Index();
	else
		return C;
}

IndexSet IndexSet::sub_set(const Index& components) const
{
	if(components.dim()==0)
		return IndexSet();
	
	IndexSet R(components.dim(),dim(),limit());

	for(int i=0; i<components.dim(); i++)
		R[i]=digit[components[i]];

	return R;
}

IndexSet IndexSet::complement_sub_set(const Index& components) const
{
	// BUG: oletetaan liikoja
	if(components.dim()==indices())
		return IndexSet();

	IndexSet R(indices()-components.dim(),dim(),limit());

	// BUG: t�t� voisi marginaalisesti parantaa
	int j=0;
	for(int i=0; i<indices(); i++)
		if(!components.has(i))
			R[j++]=digit[i];

	return R;
}

// //  CLASS: CrpIndex
// //  ===============

// CrpIndex::CrpIndex(int num_of_idxs,int dim,int max_value)
// 	:IndexSet(num_of_idxs,dim,max_value)
// {
//     while(has_common_digit())
// 		IndexSet::operator++(0);
// }

// CrpIndex& CrpIndex::operator++(int)
// {
//     do
//     {
// 		IndexSet::operator++(0);
// 		if(overflow)
// 			return *this;
		
//     } while(has_common_digit());
    
//     return *this;
// }

// CLASS: IndexIdentifier
// ======================

IndexIdentifier::IndexIdentifier()
{
	max_parts=0;
	parts=0;
	part=0;
	dimension=-1;
}

IndexIdentifier::IndexIdentifier(int dim)
{
	errif(dim < 1,"IndexIdentifier::IndexIdentifier: illegal dim " << dim);

	max_parts=dim;
	part=new Index[max_parts];
	parts=0;
	dimension=dim;
}

IndexIdentifier::~IndexIdentifier()
{
	if(max_parts)
		delete[] part;
}

IndexIdentifier& IndexIdentifier::operator=(const IndexIdentifier& Id)
{
	if(this==&Id)
		return *this;
	
	if(max_parts)
		delete[] part;

	parts=Id.parts;
	max_parts=Id.max_parts;
	dimension=Id.dimension;
	part=new Index[max_parts];

	for(int i=0; i<parts; i++)
		part[i]=Id.part[i];
	
	return *this;
}

IndexIdentifier::IndexIdentifier(const IndexIdentifier& Id)
{
	max_parts=0;
	operator=(Id);
}

int compare(const IndexIdentifier&I,const IndexIdentifier&J)
{
	if(I.parts < J.parts)
		return -1;
	if(I.parts > J.parts)
		return 1;

	for(int i=0; i<I.parts; i++)
	{
		if(I.part[i].dim() < J.part[i].dim())
			return -1;
		if(I.part[i].dim() > J.part[i].dim())
			return 1;
		
		for(int n=0; n<I.part[i].dim(); n++)
		{
			if(I.part[i][n] < J.part[i][n])
				return -1;
			if(I.part[i][n] > J.part[i][n])
				return 1;
		}
	}


	return 0;
}

ostream& operator<<(ostream& os,const IndexIdentifier& Id)
{
	os << '{';
	if(Id.parts)
	{
		for(int i=0; i<Id.parts; i++)
		{
			if(i)
				os << ',';
			
			os << '{' << Id.part[i] << '}';
		}
	}
	os << '}';
	
	return os;
}

istream& operator>>(istream& is,const IndexIdentifier&)
{
	err("operator>>(istream&,const IndexIdentifier&): unimplemented");

	return is;
}

void IndexIdentifier::show_format()
{
	if(parts==0)
	{
		//cout << '0';
		return;
	}

	for(int i=0; i<parts; i++)
	{
		//if(i) cout << '+';
		//cout << part[i].dim();
	}
}

IndexIdentifier IndexIdentifier::format() const
{
	if(parts==0)
		return IndexIdentifier();

	Index J(1,limit());
	IndexIdentifier I(dimension);

	for(int i=0; i<parts; i++)
	{
		J[0]=part[i].dim();
		I.add_to_name(J);
	}

	return I;
}

// NIMENMUODOSTUSRUTIINIT
// ----------------------

void IndexIdentifier::sort()
{
	Index tmp;
	
	for(int i=0; i<parts-1; i++)
		for(int j=i+1; j<parts; j++)
			if(part[j].dim() < part[i].dim() || (part[j].dim() == part[i].dim() && part[j] < part[i]))
			{
				tmp=part[i];
				part[i]=part[j];
				part[j]=tmp;
			}
}

// Poistetaan komponentti I, jos sellainen l�ytyy
void IndexIdentifier::delete_from_name(const Index& I)
{
	for(int i=0; i<parts; i++)
		if(part[i].dim()==I.dim() && part[i]==I)
		{
			parts--;
			for(int j=i; j<parts; j++)
				part[j]=part[j+1];
			return;
		}
}

// Lis�t��n I komponentiksi
void IndexIdentifier::add_to_name(const Index& I)
{
	errif(parts==max_parts,"IndexIdentifier::add_to_name: table full (max " << max_parts << ")");
	part[parts++]=I;
}

// Palautetaan komponenttien part[K[0]],...,part[K[n]] leikkaus
Index IndexIdentifier::intersection_of(const Index& K)
{
	Index R=part[K[0]];
	
	for(int i=1; i<K.dim(); i++)
		R=R.intersection(part[K[i]]);

	return R;
}

// Etsi ja toteuta joku 'how_many' komponentin sallittu leikkaus.
bool IndexIdentifier::apply_k_intersect(int how_many)
{
	if(parts < 2)
		return false;

	Index C;
	for(Index K(how_many,parts); K; K++)
	{
		// Lasketaan leikkaus alkioista
		C=intersection_of(K);

		// Etsit��n teoreettinen alkuper�isten indeksien m��r� 'k'
		int k=how_many * (dimension + 1);
		for(int i=0; i<K.dim(); i++)
			k-=part[K[i]].dim();

		// Toteutetaan leikkaus, mik�li se on sallittu
 		if(C.dim() > dimension - k)
 		{
 			for(int i=K.dim()-1; i>=0; i--)
 				delete_from_name(part[K[i]]);

 			add_to_name(C);
 			return true;
 		}
	}

	return false;
}

void IndexIdentifier::get(const IndexSet& I0)
{
	// Tuhotaan vanha, jos sellainen on
	if(max_parts)
		delete[] part;

	// Varataan uusi taulukko
	dimension=I0.dim();

	max_parts=I0.indices();
	part=new Index[max_parts];
	parts=I0.indices();

	// Kopioidaan alkuper�iset
	for(int i=0; i<parts; i++)
		part[i]=I0[i];
	
	// Toistetaan leikkausoperaatioita kunnes tulos ei en�� muutu
	bool modified=true;
	while(modified)
	{
		modified=false;
		for(int k=2; k<=parts; k++)
			if(apply_k_intersect(k))
			{
				modified=true;
				break;
			}
	}
	
	sort();
}

// MUUT GEOMETRISET OPERAATIOT
// ---------------------------

int IndexIdentifier::dim() const
{
	if(parts)
	{
		int d=part[0].dim()-1;
		
		for(int i=1; i<parts; i++)
		{
			d=(d + part[i].dim() - dimension - 1);
		}

		if(d < 0)
 			d=-1;
		
		return d;
	}

	return -1;
}

void IndexIdentifier::put_sup_objects(set<IndexIdentifier>& S,int d)
{
	JokerIdentifier J;
	set<JokerIdentifier> L;
	set<JokerIdentifier>::iterator i;

	J.get(*this);
	J.put_sup_identifiers(L,1);

	for(i=L.begin(); i!=L.end(); i++)
		(*i).convert_to_identifiers(S,d);
}

vector<int> IndexIdentifier::partitions() const
{
	vector<int> ret(parts);

	for(int i=0; i<parts; i++)
	{
		ret[i]=dimension-part[i].dim()+1;
	}

	return ret;
}

int IndexIdentifier::is_single_point() const
{
	if (parts != 1 || part->dim() != 1) 
		return -1;
	return part[0][0];
}

int IndexIdentifier::sup_objects() const
{
	list<vector<int> > split;
	list<vector<int> >::const_iterator s;
	vector<int> P,Q,R;
	int total=0;

	P=partitions();
	for(size_t i=0; i<P.size(); i++)
	{
		Q=P;
		Q[i]--;
		
		if(Q[i]==0)
		{
			total++;
			continue;
		}
		
		get_partitions(split,Q[i]);
		for(s=split.begin(); s!=split.end(); s++)
		{
			int mul=1;
			R=*s;
			for(size_t j=0; j<R.size(); j++)
			{
				mul*=choices(limit() - (space_dim()-P[i]+1) , P[i]-R[j]);
				if(mul < 0)
					return INT_MAX;
			}
			total+=mul;
			if(total < 0)
				return INT_MAX;
		}
	}

	return total;
}

int IndexIdentifier::real_partitions(int omax) const
{
	int count = 0;
	for (int i = 0; i < parts; i++)
	if (part[i][0] < omax)
		count++;
	return count;
}

// class JokerIdentifier
// =====================

void JokerIdentifier::joker_expand()
{
	for(int i=0; i<parts; i++)
		if(part[i].dim() < dimension)
		{
			Index I(dimension,part[i].limit());
			for(int j=0; j<dimension; j++)
				if(j < part[i].dim())
					I[j]=part[i][j];
				else
					I[j]=-1;

			int k=dimension-part[i].dim()+1;
			delete_from_name(part[i]);
			while(k--)
				add_to_name(I);
			sort();
			
			joker_expand();
			return;
		}
}

ostream& operator<<(ostream& os,const JokerIdentifier& J)
{
	os << '{';
	if(J.parts)
	{
		for(int i=0; i<J.parts; i++)
		{
			if(i)
				os << ',';
			
			os << '{';
			for(int j=0; j<J.part[i].dim(); j++)
			{
				if(j) os << ' ';
				if(J.part[i][j]==-1)
					os << '*';
				else
					os << J.part[i][j];
			}
			os << '}';
		}
	}
	os << '}';
	
	return os;
}

// Lis�t��n listaan kaikki annetun dimensioiset objektit, joihin t�m� objekti sis�ltyy.
// Huom. dim parametri oikeastaan ei m��r�� dimensiota vaan on ainoastaan vihje.
//       objektin lopullinen dimensio voi olla my�s korkeampi
void JokerIdentifier::put_sup_identifiers(set<JokerIdentifier>& L,int d) const
{
	errif(d < 0 || d > dimension-1,"JokerIdentifier::put_sub_list_of_dim: illegal dimension " << d);

	Index K(dimension-d,parts);
	
	JokerIdentifier J(*this);
	
	while(K)
	{
		J.clear();

		for(int i=0; i<K.dim(); i++)
			J.add_to_name(part[K[i]]);
		J.sort();

		if(L.count(J)==0)
			L.insert(J);
		
		K++;
	}
}

int JokerIdentifier::jokers() const
{
	int sum=0;

	for(int i=0; i<parts; i++)
		for(int j=0; j<dimension; j++)
			if(part[i][j]==-1)
				sum++;

	return sum;
}

// T�ydennet��n jokerit numeroilla, jotka eiv�t ole esiintyneet indeksiss�.
// Palautetaan tyhj� indeksi, jos ei onnistu.
bool JokerIdentifier::put(IndexSet& I) const
{
	I=IndexSet(parts,dimension,limit());

	if(jokers()==0)
	{
		for(int i=0; i<parts; i++)
			I[i]=part[i];
		
		if(I.validate())
			return true;

		I=IndexSet();
		return false;
	}
	else // jokers()!=0
	{
		set<int> S;
		
		for(int i=0; i<parts; i++)
			for(int j=0; j<part[i].dim(); j++)
				if(part[i][j]!=-1)
					S.insert(part[i][j]);

		int next=0;

		for(int i=0; i<parts; i++)
		{
			for(int j=0; j<part[i].dim(); j++)
				if(part[i][j]==-1)
				{
					while(S.count(next) > 0)
					{
						errif(next>=limit(),"JokerIdentifier::put: overflow (limit " << limit() << " with " << *this); // BUG: ei tunnu toimivan
						next++;
					}
					I[i][j]=next;
					next++;
				}
				else
					I[i][j]=part[i][j];
		}
		
		if(I.validate())
			return true;
		
		I=IndexSet();
		return false;
	}
}

// Lis�t��n joukkoon kaikki mahdolliset jokerin laajennokset.
// Listaan lis�tt�vien tunnusten dimensio tarkistetaan, mik�li dimensio annettu.
void JokerIdentifier::convert_to_identifiers(set<IndexIdentifier>& N,int d) const
{
	JokerIdentifier J;
	IndexSet S;
		
	if(jokers()==0)
	{
		N.insert(*this);
		return;
	}

	SimpleIndex I(jokers(),0,limit()-1);
	
	for(;I;I++)
	{
		J=*this;
		int i;
		int k=0;

		// BUG: heikko toteutus (int-pointteri-taulukolla nopeampi)
		for(i=0; i<J.parts; i++)
		{
			for(int j=0; j<dimension; j++)
				if(J.part[i][j]==-1)
				{
					J.part[i][j]=I[k];
					k++;
				}
			if(!J.part[i].validate())
				break;
		}
		J.sort();

		if(i==J.parts && J.put(S))
		{
			IndexIdentifier Id;
			Id.get(S);
			if(N.count(Id)==0 && (d==-1 || Id.dim()==d))
				N.insert(Id);
		}
	}
}

IndexSet JokerIdentifier::get_random_sup_object(int d) const
{
	errif(d < 0 || d > dimension-1,"IndexIdentifier::get_random_sup_object: illegal dimension " << d);
	errif(dim() > d,"IndexIdentifier::get_random_sup_object: " << *this << " is not contained in any object of dimension " << d);

	
	Index K(dimension-d,parts);
	K.random();
	
	JokerIdentifier J(dimension);
	for(int i=0; i<K.dim(); i++)
		J.add_to_name(part[K[i]]);

	IndexSet I(dimension-d,dimension,limit());

	for(;;)
	{
		int i;
		// BUG: heikko toteutus (int-pointteri-taulukolla nopeampi)
		for(i=0; i<J.parts; i++)
		{
			for(int j=0; j<dimension; j++)
 				if(J.part[i][j]==-1)
					I[i][j]=random(0,limit()-1);
				else
					I[i][j]=J.part[i][j];
		}
		if(I.validate())
			break;
	}
		
	return I;
}

int JokerIdentifier::sup_combinations(int min_dim) const
{
	err("JokerIdentifier::sup_combinations: unimplemented");
	return 0;
}

