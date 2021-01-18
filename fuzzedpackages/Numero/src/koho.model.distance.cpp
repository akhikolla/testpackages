/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
vector<mdreal>
Model::distance(const string& key) const { 
  ModelBuffer* p = (ModelBuffer*)buffer;

  /* Find the data point. */
  map<string, Point>::const_iterator pos;
  const map<string, Point>& points = p->points;
  if((pos = points.find(key)) == points.end())
    return vector<mdreal>();

  /* Check if a trainer exists. */
  if((p->trainer).size() < 1)
    p->trainer = Trainer(p->codebook, p->topology, 0, 0.0);
  
  /* Estimate distances. */
  return (p->trainer).distance(pos->second);
}
