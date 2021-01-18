#include "Fit.h"
Fit::Fit(IO *io){
  this->io = io;
  this->model = Model(io->p,io->g,io->nsample);
  this->theta.resize(io->nItEM,2*io->g+4);
};
void Fit::fitModel(){
  if(io->algorithm=="MCEM"){
    model.fitMCEM(io,theta);
  }else{
    model.fitSEM(io,theta);
  }
};
void Fit::output(){
    io->b  = model.get_b();
    io->pi = model.get_pi();
    io->intercept = model.get_intercept();
    io->sigma2 = model.get_sigma2() ;
    io->gamma2 = model.get_gamma2() ;
    io->likelihood = model.get_likelihood();
    io->entropy = model.get_entropy();
    io->theta   = theta;
    io->P       = model.get_P();
    io->Zw      = model.get_Zw();
    io->Bw      = model.get_Bw();
};
