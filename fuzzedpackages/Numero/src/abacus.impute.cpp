/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
static vector<mdsize>
find_nearest(const vector<mdreal>& x, const vector<mdreal*>& book,
	     const mdsize nsubsampl, mt19937& twister) {
  mdreal rlnan = medusa::rnan();
  vector<mdsize> output;
  
  /* Check if anything to do. */
  size_t nval = 0;
  size_t ndim = x.size();
  for(size_t j = 0; j < ndim; j++)
    nval += (x[j] != rlnan);
  if(nval == ndim) return output;
  if(nval == 0) return output;

  /* Number of neihgbors to collect. */
  size_t nbook = book.size();
  size_t nneigh = nbook;
  if(nneigh > nsubsampl) nneigh = nsubsampl;

  /* Set increment for random sampling. */
  size_t incr = nbook/nneigh;
  if(nbook > nneigh) incr++;
  
  /* Calculate neighbor distances. */
  vector<mdreal> delta;
  vector<mdsize> neighbors;
  for(size_t i = 0; neighbors.size() < nneigh; i++) {
        
    /* Euclidean distance. */
    mdreal dsum = 0.0;
    mdreal wsum = 0.0;
    const mdreal* page = book[i%nbook];
    for(size_t j = 0; j < ndim; j++) {
      mdreal bit = (x[j] != rlnan)*(page[j] != rlnan);
      mdreal d = (x[j] - page[j]);
      dsum += bit*d*d;
      wsum += bit;
    }

    /* Add new neighbor. */
    if(wsum > 0.0) {
      delta.push_back(sqrt(dsum)/wsum);
      neighbors.push_back(i%nbook);
    }

    /* Make sure loop will end. */
    i += twister()%incr;
    if(i > 2*nbook) break;
  }

  /* Sort neighbors. */
  output = medusa::sortreal(delta, 1);
  for(size_t k = 0; k < output.size(); k++)
    output[k] = neighbors[output[k]];
  return output;
}

/*
 *
 */
static bool
merge_values(vector<mdreal>& x, const mdreal* page) {
  size_t n = 0;
  mdreal rlnan = medusa::rnan();
  for(size_t j = 0; j < x.size(); j++) {
    if(x[j] == rlnan) { 
      x[j] = page[j];
      n++;
    }
  }
  return (n == 0);
}

/*
 *
 */
void
abacus::impute(vector<vector<mdreal> >& output, const mdsize nsub) {
  mdreal rlnan = medusa::rnan();
  mdsize nvect = output.size();
  if(nvect < 1) return;
  if(nsub < 1) return;

  /* Set up random numbers. */
  string seedval = (long2string(nsub) + long2string(nvect));
  seed_seq seed(seedval.begin(), seedval.end());
  mt19937 twister(seed);

  /* Find the longest vector. */
  mdsize nelem = 0;
  for(mdsize i = 0; i < output.size(); i++) {
    vector<mdreal>& values = output[i];
    for(size_t j = 0; j < values.size(); j++) {
      if((values[j] != rlnan) && (j >= nelem)) nelem = (j + 1);
    }
  }
  
  /* Allocate data buffer.*/
  vector<mdreal*> book(output.size(), NULL);
  mdreal* ptr = (mdreal*)malloc(nvect*nelem*sizeof(mdreal));
  for(mdsize i = 0; i < output.size(); i++) {
    book[i] = ptr;
    ptr += nelem;
  }
  
  /* Copy data values. */
  for(mdsize i = 0; i < output.size(); i++) {
     mdreal* page = book[i]; 
     vector<mdreal>& values = output[i];
     for(size_t k = 0; k < values.size(); k++)
       page[k] = values[k];
     for(size_t k = values.size(); k < nelem; k++)
       page[k] = rlnan;
  }
  
  /* Fill in missing values. */
  for(mdsize i = 0; i < output.size(); i++) {
    vector<mdreal>& values = output[i];
    vector<mdsize> neigh = find_nearest(values, book, nsub, twister);
    values.resize(nelem, rlnan);
    for(vector<mdsize>::iterator it = neigh.begin();
	it != neigh.end(); it++) {
      if(merge_values(values, book[*it])) break;
    }
  }

  /* Release data buffer. */
  free(book[0]);
}
