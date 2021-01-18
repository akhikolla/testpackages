/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
vector<string>
Model::identities() const {
  ModelBuffer* p = (ModelBuffer*)buffer;
  vector<string> array;
  map<string, Point>& points = p->points;
  map<string, Point>::const_iterator it;
  for(it = points.begin(); it != points.end(); it++)
    array.push_back(it->first);
  return array;
}
