#ifndef STCXEMspatial_H
#define STCXEMspatial_H
#include "STCXEM.h"


class STCXEMspatial : public STCXEM{
  public:
  
  STCXEMspatial(){};
  ~STCXEMspatial(){};  
  STCXEMspatial(const S4 &, const List &, const NumericMatrix &);

  virtual double ComputeLogLike();
  virtual void Mstep();
  virtual void NewtonLogitWeighted(const int);
};
#endif
