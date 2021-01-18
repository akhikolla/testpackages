/* 
 *  point.cpp : Latticeiterator.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: oja_iterator.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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

using namespace std; //df

/* Kombinatoriikkaa (valitse yksi) */
/* ---------------- */

/*#define STRONG_COMBINATORICS */
/* K�ytet��n t�ydellist� kombinatoriikkaa selvitett�ess� pisteen kautta */
/* kulkevia suoria. */

#define NORMAL_COMBINATORICS
/* Jos piste on havaintopiste, muodostetaan vain muiden havaintopisteiden */
/* suuntaan l�htev�t suorat. */

/*#define WEAK_COMBINATORICS */
/* Edellisen lis�ksi j�tet��n k�sittelem�tt� erikoistapaus, miss� hypertaso leikkaa kahta */
/* havaintopistett� yhdist�v�n suoran. */

// Class: OjaLineIndexIterator
// ===========================

/* Laskurin sis�isi� muuttujia k�ytet��n seuraavasti:

    int i,j;
	Index I;
	IndexSet S;

   moodi IT_NORMAL:

   S on alkuper�inen indeksijoukko, josta i:s indeksi poistetaan.

   moodi IT_SPLIT:

   Indeksi I lis�t��n aina sellaisenaan. Muihin indekseihin sis�llytet��n
   aina luvut i ja j, sek� S, joka on l�pik�yt�vien indeksien laskuri. Kun kaikki
   n�m� on k�yty l�pi, asetetaan i:=-1. Sen j�lkeen palautetaan viel� satun-
   nainen pisteidem i ja j l�pi kulkevan suoran indeksi.

   moodi IT_DATA:

   Indeksijoukok S k�y l�pi kaikki vapaiden indeksien vaihtoehdot, joista tulos
   saadaan lis��m�ll� kuhunkin indeksiin luku i.

   moodi IT_POINTOPOINT:

   K�yd��n l�pi kaikki pisteet j ja arvotaan satunnaisen suoran i--j indeksit.
   
*/
   
OjaLineIndexIterator::OjaLineIndexIterator(const OjaPoint& p, bool cheat_data
  ,bool cheat_split)
{
	errif(p.index().dim() < 2,"OjaLineIndexIterator::OjaLineIndexIterator: "
	  "invalid point '" << p << "'");
	
	overflow=false;
	current=IndexSet(p.index().indices()-1,p.index().dim(),p.index().limit());

	if(p.is_data())
	{
		if(cheat_data)
		{
			mode=IT_POINTOPOINT;
			i=p.data_index();
			j=0;
		}
		else
		{
			mode=IT_DATA;
			i=p.data_index();
			S=IndexSet(current.indices(),current.dim()-1,current.limit());
		}
	}
	else if(p.splits_line(i,j,I))
	{
		if(cheat_split)
		{
			mode=IT_NORMAL;
			i=0;
			S=p.index();
		}
		else
		{
			mode=IT_SPLIT;
			S=IndexSet(current.indices()-1,current.dim()-2,current.limit());
		}
	}
	else
	{
		mode=IT_NORMAL;
		i=0;
		S=p.index();
	}

	refresh_value();
}

void OjaLineIndexIterator::refresh_value()
{
	switch(mode)
	{
	  case IT_NORMAL:
		  for(int k=0; k<current.indices(); k++)
			  current[k] = S[k + (k >= i)];		  
		  break;

	  case IT_SPLIT:
	  {
		  int d=I.dim();

		  while(S)
		  {
			  for(int k=0; k<d-2; k++)
			  {
				  current[k][0]=i;
				  current[k][1]=j;
				  for(int n=0; n<d-2; n++)
					  current[k][n+2] = S[k][n];				  
			  }
			  
			  current[d - 2] = I;
			  
			  if(current.validate() && current.how_many_common_digits() <= 2)
				  break;
			  
			  S++;
		  }
		  
		  if(!S) // Kaikki muut paitsi suora i--j on k�yty l�pi
		  {
			  int k;
			  for(k=0; k<1000; k++)
			  {
				  for(int n=0; n<current.indices(); n++)
					  current[n].random_with(i,j);
				  
				  if(current.how_many_common_digits() <= 2)
					  break;
			  }
			  errif(k==1000,"OjaLineIndexIterator::refresh_value: "
				"giving up - not enough data points");

			  i=-1;
		  }
		  
		  break;
	  }
	  
	  case IT_DATA:

		  while(S)
		  {
			  for(int k=0; k<current.indices(); k++)
			  {
				  current[k][0]=i;
				  for(int n=0; n<current.dim()-1; n++)
					  current[k][n+1]=S[k][n];
			  }
				  
			  if(current.validate() &&
				(current.dim()<= 2 || current.how_many_common_digits() <= 2)
				  )
				  break;
			  
			  S++;
		  }				  
		  break;

	  case IT_POINTOPOINT:

		  if(i==j)
			  j++;
		  if(j < current.limit())
		  {
			  int k;
			  for(k=0; k<1000; k++)
			  {
				  for(int n=0; n<current.indices(); n++)
					  current[n].random_with(i,j);
				  
				  if(current.dim() <= 2 || current.how_many_common_digits() <= 2)
					  break;
			  }
			  errif(k==1000,"OjaLineIndexIterator::refresh_value: "
				"giving up - not enough data points");
		  }
		  break;
		  
	} // switch(mode)
}

OjaLineIndexIterator& OjaLineIndexIterator::operator++(int)
{
	errif(overflow,"OjaLineIndexIterator::operator++: overflow");
	
	switch(mode)
	{
	  case IT_NORMAL:
		  i++;
		  if(i >= S.indices())
			  overflow=true;
		  else
			  refresh_value();
		  break;

	  case IT_SPLIT:
		  if(i==-1)
			  overflow=true;
		  else
		  {
			  S++;
			  refresh_value();
		  }
		  break;

	  case IT_DATA:
		  S++;
		  if(S)
		  {
			  refresh_value();
			  if(!S)
				  overflow=true;
		  }
		  else
			  overflow=true;
		  break;

	  case IT_POINTOPOINT:
		  j++;
		  if(j>=current.limit())
			  overflow=true;
		  else
			  refresh_value();
		  break;
	}
	
	return *this;
}

// Class: OjaLineIterator
// ======================

OjaLineIterator::OjaLineIterator(const OjaPoint& p) :
#ifdef STRONG_COMBINATORICS
  OjaLineIndexIterator(p,false,false)
#endif
#ifdef NORMAL_COMBINATORICS
  OjaLineIndexIterator(p,true,false)
#endif
#ifdef WEAK_COMBINATORICS
  OjaLineIndexIterator(p,true,true)
#endif
{
	L=OjaLine(p.data_set());
}

OjaLineIterator::OjaLineIterator(const OjaPoint& p,bool chdt,bool chsp) :
  OjaLineIndexIterator(p,chdt,chsp)
{
	L=OjaLine(p.data_set());
}

const OjaLine& OjaLineIterator::line()
{
	if(L.index().dim()==0 || L.index() != index())
		L.get(index());

	return L;
}
