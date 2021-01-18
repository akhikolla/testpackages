/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"


/*
 * Returns an null pointer if capacity was not exceeded.
 * If capacity would be exceeded, removes the pointer
 * with the largest distance, and returns it.
 */
Point*
Subset::join(Point* pnt, const mdreal delta) {
  mdsize sznan = medusa::snan();
  if(occupancy > capacity)
    panic("Capacity exceeded.", __FILE__, __LINE__);
  
  /* Test capacity. */
  if(capacity == 0) return pnt;
  if(occupancy < capacity) {
    (this->contents[delta]).push_back(pnt);
    this->occupancy += 1;
    pnt->move(label);
    return NULL;
  }
  
  /* Test the largest distance, i.e. if inserting the item would
     decrease the maximum distance. */
  ContentMap::reverse_iterator it = contents.rbegin();
  mdreal dmax = it->first;
  if(delta >= dmax) {
    pnt->move(sznan);
    return pnt;
  }

  /* Remove an existing member with the lowest value. */
  vector<Point*>& batch = it->second;
  Point* prevpnt = batch[0];
  prevpnt->move(sznan);
  batch.resize(batch.size() - 1);
  if(batch.size() < 1) contents.erase(dmax);
  
  /* Add the new item. */
  (this->contents[delta]).push_back(pnt);
  pnt->move(label);
  return prevpnt;
}
