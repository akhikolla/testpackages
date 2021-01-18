/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
vector<mdsize>
find_nearest(const vector<vector<mdreal> >& book,
	     const vector<mdreal>& x, const mdsize limit,
	     mt19937& twister) {
  mdreal rlnan = medusa::rnan();
  
  /* Check if missing values. */
  mdsize nval = 0;
  mdsize ndim = x.size();
  for(mdsize j = 0; j < ndim; j++) {
    if(x[j] == rlnan) continue;
    nval++;
  }
  if(nval == ndim) return vector<mdsize>();
  if(nval == 0) return vector<mdsize>();

  /* Allocate space. */
  mdsize nsampl = book.size();
  mdsize nelem = nsampl;
  if(nelem > limit) nelem = limit;
  vector<mdreal> delta(nelem, rlnan);
  vector<mdsize> neighbors(nelem);

  /* Set increment for random sampling. */
  mdsize incr = nsampl/nelem;
  if(nsampl > nelem) incr++;

  /* Calculate neighbor distances. */
  mdsize n = 0;
  for(mdsize i = 0; i < nsampl; i++) {
    i += twister()%incr;
    if(i >= nsampl) break;
    if(n >= nelem) {
      nelem += (limit + 1);
      delta.resize(nelem, rlnan);
      neighbors.resize(nelem);
    }

    /* Euclidean distance. */
    mdsize counter = 0;
    long double dsum = 0.0;
    const vector<mdreal>& page = book[i];
    for(mdsize j = 0; j < ndim; j++) {
      if(x[j] == rlnan) continue;
      if(page[j] == rlnan) continue;
      double d = (x[j] - page[j]);
      dsum += d*d;
      counter++;
    }
    if(counter < 1) continue;
    delta[n] = sqrt(dsum)/counter;
    neighbors[n] = i; 
    n++;
  }

  /* Trim arrays. */
  delta.resize(n);
  neighbors.resize(n);

  /* Sort neighbors. */
  vector<mdsize> order = medusa::sortreal(delta, 1);
  vector<mdsize> indices(n);
  for(mdsize i = 0; i < n; i++)
    indices[i] = neighbors[order[i]];
  return indices;
}

/*
 *
 */
bool
merge_values(vector<mdreal>& x, const vector<mdreal>& t) {
  mdsize n = 0;
  mdreal rlnan = medusa::rnan();
  for(mdsize j = 0; j < t.size(); j++) {
    if(x[j] != rlnan) continue;
    x[j] = t[j];
    n++;
  }
  return (n == 0);
}

/*
 *
 */
RcppExport SEXP
nro_impute(SEXP xdata_R, SEXP nsub_R, SEXP lag_R) {
  mdsize nsub = as<mdsize>(nsub_R);
  mdreal lag = as<mdreal>(lag_R);
  time_t stamp = time(NULL);

  /* Check data. */
  vector<vector<mdreal> > xvect = nro::matrix2reals(xdata_R, 0.0);
  if(xvect.size() < 1) return CharacterVector("Empty input.");

  /* Set up random numbers. */
  string seedval = (long2string(nsub) + real2string(xvect[0][0]));
  seed_seq seed(seedval.begin(), seedval.end());
  mt19937 twister(seed);

  /* Fill in missing values. */
  time_t reset = stamp;
  vector<vector<mdreal> > data = xvect;
  for(mdsize i = 0; i < xvect.size(); i++) {
    vector<mdreal>& x = xvect[i];
    vector<mdsize> neigh = find_nearest(data, x, nsub, twister);
    for(mdsize k = 0; k < neigh.size(); k++) {
      mdsize ind = neigh[k];
      if(merge_values(x, data[ind])) break;
    }

    /* Progress message. */
    if(lag < 0.0) continue;
    if(difftime(time(NULL), reset) < lag) continue;
    string dt = medusa::time2text(difftime(time(NULL), stamp));
    Rprintf("%.1f%% in %s\n", 100*(i + 1.0)/(xvect.size()), dt.c_str());
    reset = time(NULL);
  }

  /* Reduce memory footprint. */
  data.clear();
  
  /* Return results. */
  NumericMatrix output = reals2matrix(xvect);
  return output;
}
