/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

#define N_TOPO_COLS 6

/*
 *
 */
punos::Topology
nro::reals2topology(const vector<vector<mdreal> >& vectors,
		    const mdreal sigma) {
  mdreal rlnan = medusa::rnan();

  /* Disconnected topology. */
  if(vectors.size() == 1)
    if(vectors[0].size() == 1)
      return Topology(vectors[0][0]);
  
  /* Check district coordinates. */
  mdsize ndistricts = vectors.size();
  vector<punos::District> districts(ndistricts);
  for(mdsize i = 0; i < ndistricts; i++) {
    const vector<mdreal>& x = vectors[i];
    if(x.size() < N_TOPO_COLS) return Topology();
    for(mdsize j = 0; j < N_TOPO_COLS; j++)
      if(x[j] == rlnan) return Topology();
    districts[i].x = x[0];
    districts[i].y = x[1];
    districts[i].radii.first = x[2];
    districts[i].radii.second = x[3];
    districts[i].angles.first = x[4];
    districts[i].angles.second = x[5];
  }

  /* Create map topology.*/
  vector<mdreal> epochs(1, 0.0);
  punos::Topology topo(epochs, districts);
  if(topo.rewire(sigma) == false) return Topology();
  return topo;
}
