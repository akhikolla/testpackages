/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_train(SEXP topo_R, SEXP sigma_R, SEXP codebook_R, SEXP data_R,
	  SEXP nsub_R, SEXP eq_R, SEXP lag_R) {
  mdreal rlnan = medusa::rnan();
  mdreal sigma = as<mdreal>(sigma_R);
  time_t stamp = time(NULL);
  string err;
  
  /* Input parameters. */
  mdsize nsub = as<mdsize>(nsub_R);
  mdreal eq = as<mdreal>(eq_R);
  mdreal lag = as<mdreal>(lag_R);
  
  /* Determine map topology. */
  vector<vector<mdreal> > topodata = nro::matrix2reals(topo_R, 0.0);
  punos::Topology topo = nro::reals2topology(topodata, sigma);
  if(topo.size() < 1) return CharacterVector("Unusable topology.");
  
  /* Check prototypes. */
  vector<vector<mdreal> > protos = nro::matrix2reals(codebook_R, 0.0);
  if(protos.size() != topodata.size())
    return CharacterVector("Incompatible codebook.");

  /* Check if K-means clustering. */
  if(topodata.size() == 1) {
    if(topodata[0].size() == 1)
      topo = Topology((mdsize)(topodata[0][0] + 0.5));
    topodata.clear();
    protos.clear();
  }
  
  /* Add pseudo-random jitter to ensure convergence. */
  vector<vector<mdreal> > vectors = nro::matrix2reals(data_R, 0.001);
  if(vectors.size() < 1) return CharacterVector("Too few data.");
  
  /* Check the number of columns. */
  mdsize ncols = vectors[0].size();
  if(ncols < 3) return CharacterVector("Too few dimensions.");
  
  /* Create a self-organizing map. */
  koho::Model model(topo, nsub, eq);
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
  
  /* Fit the map to the training data. */
  vector<mdreal> history;
  vector<koho::Resident> layout;
  time_t reset = stamp;
  while(true) {
    vector<mdreal> trace;
    if(lag >= 0.0) err = model.train(layout, trace, lag);
    else err = model.train(layout, trace, rlnan);
    history.insert(history.end(), trace.begin(), trace.end());
    if(err.size() > 0) return CharacterVector(err);
    if(trace.size() < 1) break;

    /* Progress message. */
    if(lag < 0.0) continue;
    if(difftime(time(NULL), reset) < lag) continue;
    string dt = medusa::time2text(difftime(time(NULL), stamp));
    Rprintf("%d cycles in %s\n", history.size(), dt.c_str());
    reset = time(NULL);
  }
  
  /* Update prototypes. */
  protos.resize(topo.size());
  for(mdsize i = 0; i < topo.size(); i++) {
    protos[i] = model.prototype(i);
    protos[i].resize(ncols, rlnan);
  }
  
  /* Recode centroid labels. */
  vector<mdsize> labels(layout.size(), 0);
  vector<mdreal> residuals(layout.size(), 0.0);
  for(mdsize i = 0; i < layout.size(); i++) {
    mdsize key = string2size(layout[i].identity);
    mdsize label = layout[i].district;
    if(label >= protos.size()) labels[key] = 0;
    else labels[key] = (label + 1); /* R-style indexing */
    residuals[key] = layout[i].residual;
  }
  
  /* Return results. */
  List res;
  res.push_back(labels, "layout");
  res.push_back(nro::reals2vector(residuals), "residuals");
  res.push_back(nro::reals2matrix(protos), "centroids");
  res.push_back(nro::reals2vector(history), "history");
  return res;
}
