#ifndef MODELSELECT_H
#define MODELSELECT_H
#include "IO.h"
#include "Model.h"

using namespace std;
using namespace Eigen;

class ModelSelect{
 private:
  int g_best;
  IO *io;
  Model *models;
  MatrixXd *thetas;
 public:
  ModelSelect(IO *io);
  ~ModelSelect(){};
  void fitAllModels();
  void findBestModel();
  void output();
};
#endif
