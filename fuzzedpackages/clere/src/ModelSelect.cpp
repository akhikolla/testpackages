#include "ModelSelect.h"
ModelSelect::ModelSelect(IO *io){
  int g;
  this->io          = io;
  this->models      = new Model[io->g];
  this->thetas      = new MatrixXd[io->g];
  for(g=1;g<=io->g;g++){
    thetas[g-1].resize(io->nItEM,2*g+4);
    models[g-1] = Model(io->p,g,io->nsample);
  }
};
void ModelSelect::fitAllModels(){
  int g;
  for(g=1;g<=io->g;g++){
    if(io->algorithm=="MCEM"){
      models[g-1].fitMCEM(io,thetas[g-1]);
    }else{
      models[g-1].fitSEM(io,thetas[g-1]);
    }
  }
};
void ModelSelect::findBestModel(){
  int g,df;
  double penalty;
  // The Model Selection Criterion is minimized
  double tmp,criterion = INFINITY;
  for(g=1;g<=io->g;g++){
    if(io->sparse){
      df = 2*g + 1;
    }else{
      df = 2*g + 2;
    }
    penalty = log(io->n) * df;
    if(io->analysis=="aic"){
      penalty = 2.0 * df;
    }
    if(io->analysis=="icl"){
      penalty += models[g-1].get_entropy();
    }
    tmp = -2*models[g-1].get_likelihood() + penalty;
    if( tmp < criterion ){
      criterion = tmp;
      g_best = g;
    }
  }
};
void ModelSelect::output(){
  int ig_best = g_best-1;
  io->b  = models[ig_best].get_b();
  io->pi = models[ig_best].get_pi();
  io->intercept = models[ig_best].get_intercept();
  io->sigma2 = models[ig_best].get_sigma2() ;
  io->gamma2 = models[ig_best].get_gamma2() ;
  io->likelihood = models[ig_best].get_likelihood();
  io->entropy = models[ig_best].get_entropy();
  io->theta   = thetas[ig_best];
  io->P       = models[ig_best].get_P();
  io->Zw      = models[ig_best].get_Zw();
  io->Bw      = models[ig_best].get_Bw();
};
