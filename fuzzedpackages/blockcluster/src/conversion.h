
/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
*/

/*
 * Project:  cocluster
 * created on: Jan 10, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file conversion.h
 *  @brief 
 **/


#ifndef CONVERSION_H_
#define CONVERSION_H_
#include <RTKpp.h>

template<class out,class inp>
inline out convertMatrix(const inp& matrixinput)
{
  int rows = matrixinput.sizeRows();
  int cols = matrixinput.sizeCols();
  out matrixOutput(rows,cols);
   for(int i=0;i<rows;i++)
   {
     for(int j=0;j<cols;j++)
     {
       matrixOutput(i,j) = matrixinput(i,j);
     }
   }
   return matrixOutput;
}

template<class inp,class out>
inline void convertMatrix(const inp& matrixinput,out& matrixOutput)
{
  int rows = matrixinput.rows();
  int cols = matrixinput.cols();
  matrixOutput = out(rows,cols);
   for(int i=0;i<rows;i++)
   {
     for(int j=0;j<cols;j++)
     { matrixOutput(i,j) = matrixinput(i,j);}
   }
}

template<class out,class inp>
inline out convertvector(const inp& vectorinput)
{
  int len = vectorinput.size();
  out vectoroutput(len);
  for (int i = 0; i < len; ++i)
  { vectoroutput[i] = vectorinput[i];}
  return vectoroutput;
}

template<class inp,class out>
inline void convertvector(const inp& vectorinput,out& vectoroutput)
{
  int len = vectorinput.size();
  vectoroutput = out(len);
  for (int i = 0; i < len; ++i)
  { vectoroutput[i] = vectorinput[i];}
}

#endif /* CONVERSION_H_ */
