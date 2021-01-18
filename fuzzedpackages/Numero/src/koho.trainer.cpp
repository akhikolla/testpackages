/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
Trainer::Trainer() {}

/*
 *
 */
Trainer::Trainer(const Matrix& codebook, const Topology& topo,
		 const mdsize ntrain, const mdreal eq) {
  
  /* Optimal subset capacities. */
  mdsize nunits = topo.size();
  vector<mdsize> subsizes(nunits, 0);
  for(mdsize i = 0; i < ntrain; i++)
    subsizes[nunits-(i%nunits)-1] += 1;

  /* Create subsets. */
  (this->subsets).resize(nunits);
  for(mdsize k = 0; k < nunits; k++) {
    mdsize cap = subsizes[k];
    mdreal rho = (exp(-5.0) - exp(-5.0*eq))/(exp(-5.0) - 1.0);
    cap += (mdsize)(rho*(ntrain - cap - nunits));
    (this->subsets[k]).configure(k, cap);
  }

  /* Check codebook size. */
  if(codebook.size() < 1) return;
  if(codebook.size() != nunits)
    panic("Incompatible inputs.", __FILE__, __LINE__);

  /* Set prototypes. */
  (this->prototypes).resize(nunits);
  for(mdsize i = 0; i < nunits; i++)
    this->prototypes[i] = codebook.row(i);
}

/*
 *
 */
Trainer::~Trainer() {}

/*
 *
 */
mdsize
Trainer::size() const {
  return prototypes.size();
}
