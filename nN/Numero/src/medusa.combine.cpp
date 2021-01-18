/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 * Set union, intersection or difference between two vectors.
 */
vector<mdsize>
medusa::combine(const vector<mdsize>& x, const vector<mdsize>& y,
		const int flag) {
  mdsize sznan = medusa::snan();
  mdsize n = (x.size() + y.size())/2;
  vector<mdsize> z; z.reserve(n);

  /* Determine data overlap. */
  unordered_map<mdsize, char> visits;
  for(mdsize i = 0; i < x.size(); i++)
    visits[x[i]] = 'x';
  for(mdsize i = 0; i < y.size(); i++) {
    mdsize ind = y[i];
    if(visits.count(ind) < 1) visits[ind] = 'y';
    if(visits[ind] != 'y') visits[ind] = 's';
  }

  /* Collect elements. */
  unordered_map<mdsize, char>::const_iterator pos;
  for(pos = visits.begin(); pos != visits.end(); pos++) {
    if(pos->first == sznan) continue;
    bool shared = (pos->second == 's');
    if((flag < 0) && shared) continue; /* difference */
    if((flag > 0) && !shared) continue; /* intersection */
    z.push_back(pos->first);
  }
  return z;
}
