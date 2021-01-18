/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
static mdsize match_exec(vector<mdsize>& x, vector<mdsize>& y,
			 const vector<string>& a, const vector<string>& b) {
  x.clear();
  y.clear();

  /* Construct lookup tables. */
  unordered_set<string> duplicates;
  map<string, mdsize> a2rank;
  map<string, mdsize> b2rank;
  for(mdsize i = 0; i < a.size(); i++) {
    const string& key = a[i];
    if(key.size() < 1) continue;
    if(a2rank.count(key) < 1) a2rank[key] = i;
    else duplicates.insert(key);
  }
  for(mdsize j = 0; j < b.size(); j++) {
    const string& key = b[j];
    if(key.size() < 1) continue;
    if(a2rank.count(key) < 1) continue;
    if(b2rank.count(key) > 0) return 0;
    if(duplicates.count(key) > 0) return 0;
    b2rank[key] = j;
  }

  /* Collect indices. */
  mdsize nmatched = b2rank.size();
  for(map<string, mdsize>::iterator it = b2rank.begin();
      it != b2rank.end(); it++) {
    x.push_back(a2rank[it->first]);
    y.push_back(it->second);
  }
  return nmatched;
}

/* Find the elements that have matching values between Vector A (3rd
   input) and Vector B (4th). The resulting arrays X (1st) and Y (2nd)
   contain indices such that A[X] = B[Y]. Empty strings are ignored.
   Duplicate matched strings are not allowed. The size of X is returned. */
mdsize
medusa::match(vector<mdsize>& x, vector<mdsize>& y,
	      const vector<string>& a, const vector<string>& b) {
  if(a.size() > b.size())
    return match_exec(y, x, b, a);
  return match_exec(x, y, a, b);
}
