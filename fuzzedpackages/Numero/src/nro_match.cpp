/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_match(SEXP codebook_R, SEXP data_R) {
  mdreal rlnan = medusa::rnan();

  /* Check prototypes. */
  vector<vector<mdreal> > protos = nro::matrix2reals(codebook_R, 0.0);
  if(protos.size() < 1) return CharacterVector("Empty codebook.");

  /* Check data. */
  vector<vector<mdreal> > vectors = nro::matrix2reals(data_R, 0.0);
  if(vectors.size() < 1) return CharacterVector("Too few data.");

  /* Set map topology. */
  punos::Topology topo(protos.size());
  
  /* Estimate coverage of valid data. */
  vector<mdreal> covers;
  mdsize ncols = vectors[0].size();
  for(mdsize i = 0; i < vectors.size(); i++) {
    mdreal nv = abacus::statistic(vectors[i], "number");
    covers.push_back(nv/ncols);
  }

  /* Create a self-organizing map. */
  koho::Model model(topo, vectors.size(), 0.0); string err;
  for(mdsize k = 0; k < protos.size(); k++) {
    err = model.configure(k, protos[k]);
    if(err.size() > 0) return CharacterVector(err);
  }

  /* Transfer data into the model. */
  for(mdsize i = 0; i < vectors.size(); i++) {
    string key = medusa::long2string(i); /* temporary identifier */
    string err = model.insert(key, vectors[i]);
    if(err.size() > 0) return CharacterVector(err);
    vectors[i].clear(); /* reduce memory footprint */
  }
  
  /* Find best matching units. */
  vector<mdsize> bmus;
  vector<mdreal> dist;  
  for(mdsize i = 0; i < vectors.size(); i++) {
    string key = medusa::long2string(i);
    vector<mdreal> delta = model.distance(key);
    mdsize bmu = abacus::extrema(delta).first;    
    if(bmu < topo.size()) {
      bmus.push_back(bmu + 1); /* R-style indexing */
      dist.push_back(delta[bmu]);
    }
    else {
      bmus.push_back(0);
      dist.push_back(rlnan);
    }
  }
  
  /* Return results. */
  List res;
  res.push_back(bmus, "DISTRICT");
  res.push_back(covers, "COVERAGE");
  res.push_back(nro::reals2vector(dist), "RESIDUAL");
  return res;
}
