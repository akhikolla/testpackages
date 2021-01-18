#ifndef RUN_H_
#define RUN_H_

#include "RankCluster.h"
#include <RcppEigen.h>

RcppExport SEXP semR(SEXP X,SEXP m,SEXP K,SEXP Qsem,SEXP Bsem,SEXP Ql,SEXP Bl,SEXP RjSE,SEXP RjM,SEXP maxTry,SEXP run,SEXP detail);


template<class out,class inp>
inline out convertMatrix(const inp& matrixinput)
{
  int rows = matrixinput.rows();
  int cols = matrixinput.cols();
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
     {
       matrixOutput(i,j) = matrixinput(i,j);
     }
   }
}

template<class out,class inp>
inline out convertvector(const inp& vectorinput)
{
    int len = vectorinput.size();
    out vectoroutput(len);
    for (int i = 0; i < len; ++i) {
      vectoroutput(i) = vectorinput(i);
    }
    return vectoroutput;
}

template<class inp,class out>
inline void convertvector(const inp& vectorinput,out& vectoroutput)
{
    int len = vectorinput.size();
    vectoroutput = out(len);
    for (int i = 0; i < len; ++i) {
      vectoroutput(i) = vectorinput(i);
    }
}


#endif /* RUN_H_ */
