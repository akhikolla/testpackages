/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<mdreal>
Approximation::parameters() const {
  vector<mdreal> prm;

  /* Distribution mode. */
  prm.push_back(mode);

  /* Parameters for Gaussian components. */
  vector<mdreal> pgP = positive.parameters();
  vector<mdreal> pgN = negative.parameters();
  prm.insert(prm.end(), pgP.begin(), pgP.end());
  prm.insert(prm.end(), pgN.begin(), pgN.end());
  return prm;
}

/*
 *
 */
bool
Approximation::parameters(const vector<mdreal>& prm) {

  /* Clear old values. */
  this->positive = Gaussian();
  this->negative = Gaussian();
  
  /* Check input. */
  if(prm.size() < 15) return false;
  for(mdsize k = 0; k < 15; k++)
    if(prm[k] == medusa::rnan()) return false;
  
  /* Distribution mode. */
  this->mode = prm[0];
  
  /* Parameters for Gaussian components. */
  vector<mdreal> pgP;
  vector<mdreal> pgN;
  for(mdsize k = 1; k <= 7; k++)
    pgP.push_back(prm[k]);
  for(mdsize k = 8; k <= 14; k++)
    pgN.push_back(prm[k]);
 
  /* Update Gaussians. */
  this->positive = Gaussian(pgP);
  this->negative = Gaussian(pgN);
  return true;
}
