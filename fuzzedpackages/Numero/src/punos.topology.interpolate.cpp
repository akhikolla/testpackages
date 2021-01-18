/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

static vector<mdsize> find_pivots(const vector<vector<mdreal> >&,
				  const mdsize);
static vector<mdreal> interp_plane(const vector<mdreal>&,
				   const vector<mdsize>&,
				   const vector<District>&,
				   const vector<mdsize>&);
static mdsize find_hermit(const vector<vector<mdreal> >&,
			  const vector<vector<mdreal> >&);
static mdreal calc_diff(const vector<vector<mdreal> >&,
			const vector<mdreal>&);

/*
 *
 */
vector<vector<mdreal> >
Topology::interpolate(const vector<vector<mdreal> >& seeds) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  mdsize nseeds = seeds.size();
  mdsize ndistricts = (p->coord).size();
  if(nseeds < 3) panic("Too few seeds.\n", __FILE__, __LINE__);
  if(seeds[0].size() < 1) panic("No data.\n", __FILE__, __LINE__);
  if(ndistricts < nseeds) panic("Too few districts.\n", __FILE__, __LINE__);

  /* Check if anything to do. */
  if(p->maxradius <= 0.0) return vector<vector<mdreal> >();
  
  /* Collect district coordinates. */
  vector<vector<mdreal> > loci(ndistricts);
  for(mdsize i = 0; i < ndistricts; i++) {
    loci[i].push_back(p->coord[i].x);
    loci[i].push_back(p->coord[i].y);
  }

  /* Spread points. */
  vector<mdsize> districtpivots = find_pivots(loci, nseeds);
  vector<mdsize> seedpivots = find_pivots(seeds, nseeds);

  /* Interpolate map according to point locations. */
  mdsize ndim = seeds[0].size();
  vector<vector<mdreal> > proto(ndistricts);
  for(mdsize j = 0; j < ndim; j++) {

    /* Collect data column. */
    vector<mdreal> x(nseeds);
    for(mdsize i = 0; i < nseeds; i++)
      x[i] = seeds[i][j];

    /* Update component plane. */
    vector<mdreal> y;
    y = interp_plane(x, seedpivots, p->coord, districtpivots);
    for(mdsize i = 0; i < ndistricts; i++)
      proto[i].push_back(y[i]);
  }
  return proto;
}

/*
 *
 */
vector<mdsize>
find_pivots(const vector<vector<mdreal> >& vectors0,
	    const mdsize nlimit) {
  mdsize sznan = medusa::snan();
  if(nlimit < 2) panic("Unusable input.", __FILE__, __LINE__);

  /* Local data structures. */
  vector<vector<mdreal> > vectors = vectors0;
  vector<vector<mdreal> > accepted;
  vector<mdsize> res;

  /* Find the most isolated vector. */
  mdsize pos = find_hermit(vectors, vectors);
  if(pos >= sznan) panic("Inconsistent state.", __FILE__, __LINE__);
  accepted.push_back(vectors[pos]);
  vectors[pos].clear();
  res.push_back(pos);

  /* Find the next most isolated vectors. */
  while(accepted.size() < nlimit) {
    mdsize pos = find_hermit(vectors, accepted);
    if(pos == sznan)
      panic("Inconsistent state.", __FILE__, __LINE__);
    accepted.push_back(vectors[pos]);
    vectors[pos].clear();
    res.push_back(pos);
  }
  return res;
}

/*
 *
 */
vector<mdreal>
interp_plane(const vector<mdreal>& data, const vector<mdsize>& dpivots,
	     const vector<District>& districts, const vector<mdsize>& upivots) {
  mdreal rlnan = medusa::rnan();

  /* Check inputs. */
  mdsize ndistricts = districts.size();
  mdsize npivots = dpivots.size();
  if(upivots.size() != dpivots.size())
    panic("Inconsistent state.", __FILE__, __LINE__);

  /* Distance-based weights for map values. */
  vector<mdreal> plane(ndistricts, rlnan);
  for(mdsize i = 0; i < ndistricts; i++) {
    double wsum = 0.0;
    for(mdsize k = 0; k < npivots; k++) {
      mdsize dpos = dpivots[k];
      mdsize upos = upivots[k];
      mdreal dx = (districts[i].x - districts[upos].x);
      mdreal dy = (districts[i].y - districts[upos].y);
      mdreal w = 1.0/(dx*dx + dy*dy + 0.1);
      if(plane[i] == rlnan) plane[i] = 0.0;
      plane[i] += w*(data[dpos]);
      wsum += w;
    }
    if(wsum <= 0.0) continue;
    plane[i] /= wsum;
  }
  return plane;
}

/*
 *
 */
mdsize
find_hermit(const vector<vector<mdreal> >& candidates,
	    const vector<vector<mdreal> >& selected) {
  mdsize pos = medusa::snan();
  mdreal dmax = 0.0;
  for(mdsize i = 0; i < candidates.size(); i++) {
    mdreal d = calc_diff(selected, candidates[i]);
    if(d < dmax) continue;
    dmax = d;
    pos = i;
  }
  return pos;
}

/*
 *
 */
mdreal
calc_diff(const vector<vector<mdreal> >& vectors,
	  const vector<mdreal>& x) {
  mdreal rlnan = medusa::rnan();
  mdsize ndim = x.size();
  if(ndim < 1) return -1.0;

  /* Calculate sum of squared differences. */
  mdsize nsum = 0;
  long double dsum = 0.0;
  for(mdsize i = 0; i < vectors.size(); i++) {
    const vector<mdreal>& y = vectors[i];
    if(y.size() < 1) continue;
    if(y.size() != ndim) panic("Inconsistent state.", __FILE__, __LINE__);
    for(mdsize j = 0; j < ndim; j++) {
      if(x[j] == rlnan) continue;
      if(y[j] == rlnan) continue;
      double d = (y[j] - x[j]);
      dsum += d*d;
      nsum++;
    }
  }

  /* Estimate average distance. */
  if(nsum < 1) return -1.0;
  return sqrt(dsum/nsum);
}
