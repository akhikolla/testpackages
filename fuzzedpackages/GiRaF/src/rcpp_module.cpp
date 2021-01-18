#include <GiRaF.hpp>

using namespace Rcpp;

RCPP_MODULE(modGiRaF){
  
  class_<Border>("Border")
  .constructor< unsigned, unsigned, unsigned, arma::vec >()
  .method("setBorders", &Border::set_borders)
  ;
  
  class_<Lattice>("Lattice")
    .constructor< unsigned, unsigned, unsigned, unsigned,
  arma::vec >()
    .constructor< unsigned, unsigned, unsigned, unsigned,
  arma::vec, arma::vec >()
    .constructor< unsigned, unsigned, unsigned, unsigned,
  arma::vec, arma::vec, std::vector<unsigned> >()
    .method("vertices", &Lattice::view_vertices)
    .method("GibbsSampler", &Lattice::GibbsSampler)
    .method("GibbsSamplerCond", &Lattice::GibbsSamplerCond)
    .method("SWSampler", &Lattice::SWSampler)
    .method("SWSamplerCond", &Lattice::SWSamplerCond)
  ;
  
  class_<Block>("Block")
    .derives<Lattice>("Lattice")
    .constructor< unsigned, unsigned, unsigned, unsigned,
  arma::vec >()
    .constructor< unsigned, unsigned, unsigned, unsigned,
  arma::vec, arma::vec >()
    .constructor< unsigned, unsigned, unsigned, unsigned,
  arma::vec, arma::vec, std::vector<unsigned> >()
    
    .method("initFactor", &Block::initFactor)
    .method("correctFactor", &Block::correctFactor)
    .method("recursion", &Block::recursion)
    .method("recursion.mem", &Block::recursion_mem)
    .method("recursion.cond", &Block::recursion_cond)
    .method("recursion.cond.mem", &Block::recursion_cond_mem)
    .method("sample", &Block::exact_sample)
    .method("sampleCond", &Block::exact_sample_cond)
    ;
}

