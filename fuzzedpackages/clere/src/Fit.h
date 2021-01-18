#ifndef FIT_H
#define FIT_H
#include "IO.h"
#include "Model.h"
using namespace std;
using namespace Eigen;
class Fit{
 private:
  IO *io;
  Model model;
  MatrixXd theta;
 public:
  Fit(IO *io);
  ~Fit(){};
  void fitModel();
  void output();
};
#endif
