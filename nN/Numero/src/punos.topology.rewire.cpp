/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

#define R_CUTOFF 2.5
#define NLINK_CAP 2097152

/*
 *
 */
static unordered_set<mdsize>
find_neighbors(const map<mdreal, vector<mdsize> >& table,
	       const mdreal x, const mdreal sigma) {
  unordered_set<mdsize> neighbors;
	       
  /* Find self. */
  map<mdreal, vector<mdsize> >::const_iterator it;
  if((it = table.find(x)) == table.end())
    medusa::panic("Bad input.", __FILE__, __LINE__);
  
  /* Find neighbors with smaller values. */
  while(it != table.begin()) { it--;
    mdreal dx = (x - it->first)/(sigma + 1e-9);
    if(dx > R_CUTOFF) break;
    neighbors.insert((it->second).begin(), (it->second).end());
  }

  /* Find neighbors with equal or larger values. */
  for(it = table.find(x); it != table.end(); it++) {
    mdreal dx = (it->first - x)/(sigma + 1e-9);
    if(dx > R_CUTOFF) break;
    neighbors.insert((it->second).begin(), (it->second).end());
  }
  return neighbors;
}

/*
 *
 */
bool
Topology::rewire(const mdreal s) {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  
  /* Check if anything to do. */
  if(p->maxradius <= 0.0) return false;

  /* Check input. */
  if((s < 0.0) || (s >= p->maxradius)) {
    worry("Unusable input.", __FILE__);
    return false;
  }

  /* Check if changes needed. */
  vector<LinkMap>& network = p->network;
  if(network.size() < 1) p->sigma = 0.0;
  if(s == p->sigma) return true;
  network.clear();

  /* No neighbors. */
  if(s <= 0.0) {
    p->sigma = 0.0;
    return true;
  }
  
  /* Organize coordinate data. */
  map<mdreal, vector<mdsize> > xtable;
  map<mdreal, vector<mdsize> > ytable;
  const vector<District>& coord = p->coord;
  for(mdsize i = 0; i < coord.size(); i++) {
    xtable[coord[i].x].push_back(i);
    ytable[coord[i].y].push_back(i);
  }
  
  /* Calculate (non-normalized) Gaussian link weights. */
  mdbyte btmin = medusa::bmin();
  mdbyte btmax = medusa::bmax();
  for(mdsize i = 0; i < coord.size(); i++) {
    mdreal x = coord[i].x;
    mdreal y = coord[i].y;
    
    /* Find neighbors. */
    unordered_set<mdsize> xneigh = find_neighbors(xtable, x, s); 
    unordered_set<mdsize> yneigh = find_neighbors(ytable, y, s); 
    
    /* Set neighborhood topology. */
    LinkMap links;
    for(unordered_set<mdsize>::iterator it = xneigh.begin();
	it != xneigh.end(); it++) {
      if(yneigh.count(*it) < 1) continue;
      const District& district = coord[*it];
      mdreal dx = (district.x - x)/(s + 1e-9);
      mdreal dy = (district.y - y)/(s + 1e-9);
      mdreal r = sqrt(dx*dx + dy*dy);
      if(r > R_CUTOFF) continue;
      long w = (long)(exp(-0.5*r*r)*btmax + 0.5);
      if(w > btmax) w = btmax;
      if(w < btmin) continue;
      links[*it] = w;
    }
    network.push_back(links);

    /* Protect against memory crunch. */
    if(network.size() <= NLINK_CAP) continue;
    medusa::panic("Capacity exceeded.", __FILE__, __LINE__);
  }
  
  /* Save new scale factor. */
  p->sigma = s;
  return (network.size() > 0);
}
