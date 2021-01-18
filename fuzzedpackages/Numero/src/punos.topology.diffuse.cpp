/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
vector<mdreal>
Topology::diffuse(const vector<mdsize>& bmus,
		  const vector<mdreal>& vals) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  mdsize nelem = bmus.size();
  mdsize nunits = (p->coord).size();
  mdreal rlnan = medusa::rnan();
  if(vals.size() != nelem)
    panic("Incompatible inputs.", __FILE__, __LINE__);
  
  /* Calculate running sums. */
  vector<mdreal> xsums(nunits, 0.0);
  for(mdsize i = 0; i < nelem; i++) {
    mdsize bmu = bmus[i];
    if(bmu >= nunits) continue;
    mdreal x = vals[i];
    if(x == rlnan) continue;
    xsums[bmu] += x;
  }

  /* Return smoothed sums. */
  return punos_local::smoothen(xsums, p->network);
}
  
/*
 *
 */
vector<vector<mdreal> >
Topology::diffuse(const vector<Site>& layers, const vector<mdsize>& bmus,
		  const vector<mdreal>& vals) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  mdsize nelem = bmus.size();
  mdsize nunits = (p->coord).size();
  mdsize nlevels = (p->levels).size();
  mdreal rlnan = medusa::rnan();
  if(layers.size() != nelem)
    panic("Incompatible inputs.", __FILE__, __LINE__);
  if(vals.size() != nelem)
    panic("Incompatible inputs.", __FILE__, __LINE__);

  /* Check if anything to do. */
  if(p->maxradius <= 0.0) return vector<vector<mdreal> >();
  if(nlevels < 1) return vector<vector<mdreal> >();
  
  /* Allocate temporary arrays. */
  vector<mdreal> empty(nunits, 0.0);
  vector<vector<mdreal> > xsums(nlevels, empty);

  /* Calculate running sums. */
  for(mdsize i = 0; i < nelem; i++) {
    const Site& layer = layers[i];
    mdsize bmu = bmus[i];
    if(bmu >= nunits) continue;
    mdreal x = vals[i];
    if(x == rlnan) continue;
    mdsize a = layer.bounds.first;
    mdsize b = layer.bounds.second;
    mdreal wa = layer.weights.first;
    mdreal wb = layer.weights.second;
    if(a >= nlevels) continue;
    if(b >= nlevels) continue;
    if(wa == rlnan) continue;
    if(wb == rlnan) continue;
    if(wa != 0.0) xsums[a][bmu] += wa*x;
    if(wb != 0.0) xsums[b][bmu] += wb*x;
  }

  /* Return smoothed sums. */
  vector<vector<mdreal> > planes(nlevels);
  for(mdsize k = 0; k < nlevels; k++)
    planes[k] = punos_local::smoothen(xsums[k], p->network);
  return planes;
}
