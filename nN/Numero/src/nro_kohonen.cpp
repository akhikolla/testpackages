/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_kohonen(SEXP seeds_R, SEXP rho_R, SEXP sigma_R) {
  vector<vector<mdreal> > seeds = nro::matrix2reals(seeds_R, 0.0);
  mdreal rho = as<mdreal>(rho_R); 
  mdreal sigma = as<mdreal>(sigma_R);
  
  /* Check inputs. */
  if(seeds.size() < 3)
    return CharacterVector("Too few seeds.");
  if(seeds[0].size() < 3)
    return CharacterVector("Too few data columns.");
  for(mdsize i = 0; i < seeds.size(); i++) {
    mdreal s = abacus::statistic(seeds[i], "number");
    if(s < 1) return CharacterVector("Empty seed.");
  }
  if(rho < 2.0) return CharacterVector("Too small map radius.");

  /* Create map topology. */
  vector<mdreal> epochs(1, 0.0);
  punos::Topology topo(epochs, (mdsize)(rho + 0.5));
  mdsize ndistricts = topo.size();
  if(ndistricts < 1) return CharacterVector("Cannot create topology.");

  /* Set neighborhood network. */
  if(topo.rewire(sigma) == false)
      return CharacterVector("Topology failed.");
  
  /* Interpolate component planes. */
  vector<vector<mdreal> > protos = topo.interpolate(seeds);
  if(protos.size() != ndistricts)
    return CharacterVector("Prototype interpolation failed.");

  /* Collect coordinate data. */
  vector<vector<mdreal> > coord(ndistricts);
  for(mdsize i = 0; i < ndistricts; i++) {
    punos::District u = topo[i];
    vector<mdreal>& c = coord[i];
    c.push_back(u.x);
    c.push_back(u.y);
    c.push_back(u.radii.first);
    c.push_back(u.radii.second);
    c.push_back(u.angles.first);
    c.push_back(u.angles.second);
  }
  
  /* Return results. */
  List res;
  res.push_back(reals2matrix(protos), "centroids");
  res.push_back(reals2matrix(coord), "topology");
  return res;
}
