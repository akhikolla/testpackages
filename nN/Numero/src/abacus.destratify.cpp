/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

#define N_DISTINCT_min 10

/*
 *
 */
class Group {
public:
  vector<mdsize> keys;
  vector<mdreal> values;
public:
  void add(const mdsize key, const mdreal value) {
    keys.push_back(key);
    values.push_back(value);
  };
};

/*
 *
 */
static vector<mdreal>
prep_bg(unordered_map<mdsize, Group>& groups,
	const vector<mdreal>& data) {

  /* Original range. */
  mdreal xmin = statistic(data, "min");
  mdreal xmax = statistic(data, "max");

  /* Calculate centers. */
  vector<mdsize> keys;
  vector<mdreal> weights;
  vector<mdreal> centers;
  for(unordered_map<mdsize, Group>::iterator it = groups.begin();
      it != groups.end(); it++) {    
    Group& g = it->second;
    mdreal mu = statistic(g.values, "center");
    keys.push_back(it->first);
    weights.push_back(g.values.size());
    centers.push_back(mu);
  }

  /* Collect shifted values. */
  vector<mdreal> bg;
  mdreal origin = statistic(centers, weights, "center");
  for(mdsize k = 0; k < keys.size(); k++) {
    Group& g = groups[keys[k]];
    vector<mdreal>& vals = g.values;
    for(mdsize i = 0; i < vals.size(); i++) {
      mdreal x = (vals[i] - centers[k] + origin);
      if(x < xmin) continue;
      if(x > xmax) continue;
      bg.push_back(x);
    }
  }
  
  /* Sort values. */
  sort(bg.begin(), bg.end());
  return bg;
}

/*
 * X must be sorted and must contain only valid values.
 */
static vector<mdreal>
calc_quantiles(const vector<mdreal>& bg, const vector<mdreal>& q) {
  mdsize n = bg.size();
  vector<mdreal> y = q;
  mdreal rlnan = medusa::rnan();
  if(n < 2) panic("Unusable background.", __FILE__, __LINE__);

  /* Check that there is enough diversity. */
  unordered_set<mdreal> distinct(q.begin(), q.end());
  distinct.erase(rlnan);
  if(distinct.size() < N_DISTINCT_min)
    return vector<mdreal>(q.size(), rlnan);

  /* Transform quantiles to values. */
  for(mdsize i = 0; i < q.size(); i++) {
    double qi = q[i];
    if(qi == rlnan) continue;
    if(qi < 0.0) panic("Invalid quantile.", __FILE__, __LINE__);
    if(qi > 1.0) panic("Invalid quantile.", __FILE__, __LINE__);
    
    /* Determine quantile indices. */
    mdsize a = (mdsize)(qi*(n - 1));
    mdsize b = (a + 1);
    if(b >= n) b = (n - 1);
    if(a == b) {
      y[i] = bg[a];
      continue;
    }
    
    /* Interpolate between indices. */
    double d = (qi*(n - 1) - a);
    y[i] = ((1.0 - d)*(bg[a]) + d*(bg[b]));
  }
  return y;
}

/*
 *
 */
vector<mdreal>
abacus::destratify(const vector<mdreal>& x, const vector<mdsize>& g) {
  mdreal rlnan = medusa::rnan();

  /* Check inputs. */
  if(x.size() != g.size())
    panic("Incompatible arguments.", __FILE__, __LINE__);

  /* Divide values into groups. */
  unordered_map<mdsize, Group> groups;
  for(mdsize i = 0; i < x.size(); i++) {
    if(x[i] == rlnan) continue;
    groups[g[i]].add(i, x[i]);
  }
  if(groups.size() < 2) return x;

  /* Collect background. */
  vector<mdreal> y = prep_bg(groups, x);
  if(y.size() < 2) return x;

  /* Convert group values to global values. */
  unordered_map<mdsize, Group>::iterator pos;
  for(pos = groups.begin(); pos != groups.end(); pos++) {
    vector<mdreal> tmp = transform((pos->second).values, "uniform");
    (pos->second).values = calc_quantiles(y, tmp);
  }

  /* Collect results. */
  y = x;
  for(pos = groups.begin(); pos != groups.end(); pos++) {
    vector<mdreal>& values  = (pos->second).values;
    vector<mdsize>& keys = (pos->second).keys;
    for(mdsize k = 0; k < keys.size(); k++)
      y[keys[k]] = values[k];
  }
  return y;
}

