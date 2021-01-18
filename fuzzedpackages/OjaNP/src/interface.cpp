/* 
 *  interface.cpp: Interface to R.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: sinterface.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
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

#include <R.h>
#include "interface.h"
#include "oja_geometry.h"
#include "matrix_wrapper.h"
#include "global.h"
#include "random.h"



extern "C"
{
	//
	//  This is a interface to oja median.
	//
	//  rows - number of input data points
	//  cols - dimension of the data
	//  data - data itself
	//  vec_out - output vector depending on called function (median etc.)
	//  mat_out - output matrix depending on called function (covariance etc.)
	//  param1,param2,param3,param4 - miscellanous parameters to called function
    //  func - function to call:
	//         1 = calculate oja median (vec_out) using exact algorithm
	//             param1 = maximum number of lines to follow.
	//         2 = calculate oja median (vec_out) using lattice approximation
	//             param1 = maximum size of the confidence area (diameter)
	//             param2 = chi2 limit value
	//             param3 = number of hyperplanes to sample between lattice updates
	//         3 = calculate value of the oja object function at point param1 (vector)
	//             return value is first component of the output vector
	//         4 = calculate bootstrap estimate of covariance
	//             param1, param2, param3 = same as in case 2
	//             param4 = number of bootstrap samples
	//  dbg - if non-zero, show debugging information
   
	//XXXvoid r_oja(long* rows,long* cols,double* data,double* vec_out,double* mat_out,long* func,double* param1, double* param2, long* param3, long* param4,long* dbg, long* rSeed)
	void r_oja(long* rows,long* cols,double* data,double* vec_out,double* mat_out,long* func,double* param1, double* param2, long* param3, long* param4,long* dbg)
	{
		int dim=(int)*cols;
		int size=(int)*rows;
		double* d=data;
		OjaPoint v_output;

		for(int i=0; i<dim; i++)
			vec_out[i]=0;

		OjaData D(dim,size);
		for(int i=0; i<dim; i++)
			for(int j=0; j<size; j++)
				D[j][i]=*d++;

		debug=*dbg;
		verbose=(*dbg != 0);
		
		set_random_seed();
		
		switch((int)*func)
		{
		  case 1:
		  {
			  D.set_median_method(FOLLOW_INTERSECTION_LINES);
			  D.set_max_searchlines(int(*param1));
			  v_output=D.median();
			  break;
		  }

		  case 6:		//AP method
		  {
		    D.set_median_method(FOLLOW_INTERSECTION_LINES_BOUNDED_APPROX);  // This is a test, if this return fixed the valgrind message on this line! Check with Oleksii if the return is correct!!! (DF)
			  if (int(*param3)==1)	D.set_median_method(FOLLOW_INTERSECTION_LINES_BOUNDED);
			  D.set_max_searchlines(int(*param1));
			  D.set_volume(double(*param2));
			  v_output = D.median();
			  break;
		  }
			  
		  case 2:
		  {      //XXX if (*rSeed!=0){ srand(*rSeed);}
		  	  //GetRNGstate();
		  	  
			  D.set_median_method(LATTICE_APPROX3);
			  D.set_lattice_measure(LM_DIAMETER);
			  D.set_epsilon(*param1);
			  if(*param2 < 0)
				  D.set_chi2_limit(- *param2);
			  else
				  D.set_chi2_limit(limitChi2(D.dim(),*param2));
			  D.set_set_size(int(*param3));
			  v_output=D.median();
			  break;
		  }
		  
		  case 3:
		  {
			  Point x(dim);
			  for(int i=0; i<dim; i++)
				  x[i]=*param1++;
			  vec_out[0]=D.oja(x);
			  for(int i=1; i<dim; i++)
				  vec_out[i]=0.0;
			  return;
		  }

		  case 4:
		  {
			  D.set_median_method(BOOTSTRAP);
			  D.set_lattice_measure(LM_DIAMETER);
			  D.set_epsilon(*param1);
			  if(*param2 < 0)
				  D.set_chi2_limit(- *param2);
			  else
				  D.set_chi2_limit(limitChi2(D.dim(),*param2));
			  D.set_set_size(int(*param3));

			  matrix cov(dim,dim);
			  v_output = D.medianBootstrap(cov,int(*param4));
			  for(int i=0; i<dim; i++)
				  for(int j=0; j<dim; j++)
					  *mat_out++=cov(j,i);
			  break;
		  }

		  default:
			//  cerr << "oja: unsupported function " << *func << endl;
			  return;
		}

		for(int i=0; i<dim; i++)
			vec_out[i]=v_output.location()[i];
	}
}
