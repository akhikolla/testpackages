/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
string
Engine::insert(const string& key, const mdsize unit,
	       const mdreal value) {
  vector<mdreal> array(1, value);
  return this->insert(key, unit, array);
}

/*
 *
 */
string
Engine::insert(const string& key, const mdsize unit,
	       const vector<mdreal>& values) {
  EngineBuffer* p = (EngineBuffer*)buffer;
  mdreal rlnan = medusa::rnan();
  
  /* Check inputs. */
  if(key.size() < 1) return "Empty identity.";
  if(unit >= (p->topology).size()) return "Unusable map unit.";

  /* Check if any unusable values. */
  mdsize nvalid = 0;
  for(mdsize j = 0; j < values.size(); j++)
    nvalid += (values[j] != rlnan);
  if(nvalid < 1) return "No usable data.";
 
  /* Check that dimensions match. */
  if(p->order < 1) p->order = values.size();
  if(values.size() != p->order)
    return "Incompatible input.";
			       
  /* Insert a new data point. */
  mdsize rank = (p->points).size();
  p->points[key] = Point(rank, values, unit);
  
  /* Reset engine state. */
  if(nvalid < values.size()) p->complete = false;
  (p->cache).clear();
  return "";
}

