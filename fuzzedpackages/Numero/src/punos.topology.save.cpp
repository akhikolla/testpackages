/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
unsigned long
Topology::save(const string& fname) const {
  unsigned long n = 0;
  mdreal rlnan = medusa::rnan();
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  
  /* Open file. */
  File f; f.open(fname, "w");
  if(f.error().size() > 0) return 0;

  /* Save levels. */
  n += f.write("LEVEL\n");
  for(mdsize i = 0; i < (p->levels).size(); i++)
    n += f.write(real2string(p->levels[i]) + "\n");
  
  /* Save coordinate headings. */
  vector<string> array(7);
  array[0] = "\nDISTRICT";
  array[1] = "X";
  array[2] = "Y";
  array[3] = "RADIUSa";
  array[4] = "RADIUSb";
  array[5] = "ANGLEa";
  array[6] = "ANGLEb";
  n += f.write(array, '\t');

  /* Save coordinate data. */
  for(mdsize i = 0; i < (p->coord).size(); i++) {
    District district = p->coord[i];
    if(district.x == rlnan) panic("Unusable district.", __FILE__, __LINE__);
    array[0] = long2string(i);
    array[1] = real2string(district.x);
    array[2] = real2string(district.y);
    array[3] = real2string(district.radii.first);
    array[4] = real2string(district.radii.second);
    array[5] = real2string(district.angles.first);
    array[6] = real2string(district.angles.second);
    n += f.write(array, '\t');
  }

  /* Save neighborhood radius. */
  n += f.write("\nSIGMA\n" + real2string(p->sigma) + "\n");
  return n;
}
