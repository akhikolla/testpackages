/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
string
Model::configure(const mdsize unit, const vector<mdreal>& values) {
  ModelBuffer* p = (ModelBuffer*)buffer;
  mdreal rlnan = medusa::rnan();

  /* Check inputs. */
  mdsize nvalid = 0;
  for(mdsize j = 0; j < values.size(); j++)
    if(values[j] != rlnan) nvalid++;
  if(nvalid < 1) return "Unusable codebook."; 
  if(unit >= (p->topology).size()) return "Unusable unit index.";
  
  /* Replace codebook elements. */
  for(mdsize j = 0; j < values.size(); j++)
    (p->codebook).insert(unit, j, values[j]);

  /* Reset training state. */
  (p->state).clear();
  return "";
}
