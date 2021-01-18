#ifndef STCXEMnonspatial_H
#define STCXEMnonspatial_H
#include "STCXEM.h"


class STCXEMnonspatial : public STCXEM{
  public:
  STCXEMnonspatial(){};
  ~STCXEMnonspatial(){};  
  STCXEMnonspatial(const S4 &, const List &, const NumericMatrix &);

  virtual double ComputeLogLike();
  virtual void Mstep();
  virtual void NewtonLogitWeighted(const int);
};
#endif
