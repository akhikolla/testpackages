/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
Topology::Topology() {
  buffer = new TopologyBuffer();
}

/*
 *
 */
Topology::Topology(const mdsize k) {
  TopologyBuffer* p = new TopologyBuffer();

  /* Unconnected centroids. */
  vector<District>& coord = p->coord; coord.resize(k);
  for(mdsize i = 0; i < k; i++) {
    coord[i].x = 0.0;
    coord[i].y = 0.0;
    coord[i].radii.first = 0.0;
    coord[i].radii.second = 0.0;
    coord[i].angles.first = 0.0;
    coord[i].angles.second = 0.0;
  }

  /* Set default parameters. */
  p->maxradius = 0.0;
  p->sigma = 0.0;
  this->buffer = p;
}

/*
 *
 */
Topology::Topology(const vector<mdreal>& zpos,
		   const mdsize ncircles) {
  TopologyBuffer* p = new TopologyBuffer();
  mdreal rlnan = medusa::rnan();

  /* Check level positions. */
  if(zpos.size() < 1) return;
  if(zpos[0] == rlnan) return;
  for(mdsize k = 1; k < zpos.size(); k++) {
    if(zpos[k] <= zpos[k-1]) return;
    if(zpos[k] == rlnan) return;
  }

  /* First district. */
  District district;
  double rad = 1.05/sqrt(M_PI);
  district.x = 0.0;
  district.y = 0.0;
  district.radii.first = 0.0;
  district.radii.second = rad;
  district.angles.first = 0.0;
  district.angles.second = 360;
  vector<District> districts(1, district);
  
  /* Create a concentric circular lattice. */
  while(rad < ncircles) {
    double rA = rad;
    double rB = (rA + 1.0);
    mdsize n_div = (mdsize)(M_PI*(rB*rB - rA*rA) + 0.5);
    district.radii.first = rA;
    district.radii.second = rB;
    for(mdsize i = 0; i < n_div; i++) {
      double phiA = (360.0*i/n_div - (11*n_div)%180);
      double phiB = (360.0*(i + 1)/n_div - (11*n_div)%180);
      district.x = 0.5*(rA + rB)*cos(M_PI*(phiA + phiB)/360.0);
      district.y = 0.5*(rA + rB)*sin(M_PI*(phiA + phiB)/360.0);
      district.angles.first = phiA;
      district.angles.second = phiB;
      districts.push_back(district);
    }
    rad = rB;
    if(districts.size() >= USHRT_MAX)
      panic("Too many map districts.", __FILE__, __LINE__);
  }

  /* Fine-tune coordinates to ensure correct surface area. */
  mdsize ndistricts = districts.size();
  double scale = sqrt(ndistricts/M_PI)/rad;
  for(mdsize i = 0; i < ndistricts; i++) {
    District& u = districts[i];
    u.x *= scale;
    u.y *= scale;
    u.radii.first *= scale;
    u.radii.second *= scale;
  }

  /* Determine maximum radius. */
  mdreal rmax = rlnan;
  for(mdsize i = 0; i < ndistricts; i++) {
    mdreal r = districts[i].radii.second;
    if(rmax == rlnan) rmax = r;
    if(rmax < r) rmax = r;
  }

  /* Update object. */
  p->maxradius = rmax;
  p->levels = zpos;
  p->coord = districts;
  this->buffer = p;
}

/*
 *
 */
Topology::Topology(const vector<mdreal>& zpos,
		   const vector<District>& districts) {
  TopologyBuffer* p = new TopologyBuffer();
  mdreal rlnan = medusa::rnan();

  /* Check level positions. */
  if(zpos.size() < 1) return;
  if(zpos[0] == rlnan) return;
  for(mdsize k = 1; k < zpos.size(); k++) {
    if(zpos[k] <= zpos[k-1]) return;
    if(zpos[k] == rlnan) return;
  }

  /* Determine maximum radius. */
  mdreal rmax = rlnan;
  for(mdsize i = 0; i < districts.size(); i++) {
    mdreal r = districts[i].radii.second;
    if(rmax == rlnan) rmax = r;
    if(rmax < r) rmax = r;
  }

  /* Update object. */
  p->maxradius = rmax;
  p->levels = zpos;
  p->coord = districts;
  this->buffer = p;
}

/*
 *
 */
Topology::Topology(const Topology& t) {
  this->buffer = new TopologyBuffer(t.buffer);
}

/*
 *
 */
void
Topology::operator=(const Topology& t) {
  if(this == &t) return;
  TopologyBuffer* p = (TopologyBuffer*)buffer; delete p;
  this->buffer = new TopologyBuffer(t.buffer);
}

/*
 *
 */
Topology::~Topology() {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  delete p;
}
