#ifndef CONVERSION_H_
#define CONVERSION_H_
#include <RcppEigen.h>

template<class inp,class out>
inline void convertMatrix(const inp& matrixInput,out& matrixOutput){
  int rows = matrixInput.rows();
  int cols = matrixInput.cols();
  matrixOutput = out(rows,cols);
   for(int i=0;i<rows;i++)
   {
     for(int j=0;j<cols;j++)
     {
       matrixOutput(i,j) = matrixInput(i,j);
     }
   }
}
template<class inp,class out>
inline void convertVector(const inp& vectorInput,out& vectorOutput)
{
    int len = vectorInput.size();
    vectorOutput = out(len);
    for (int i = 0; i < len; ++i) {
      vectorOutput(i) = vectorInput(i);
    }
}

// Output conversion...thanks to Parmeet
template<class out,class inp>
inline out outMatrix(const inp& matrixinput)
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
template<class out,class inp>
inline out outVector(const inp& vectorinput)
{
    int len = vectorinput.size();
    out vectoroutput(len);
    for (int i = 0; i < len; ++i) {
      vectoroutput(i) = vectorinput(i);
    }
    return vectoroutput;
}

#endif /* CONVERSION_H_ */
