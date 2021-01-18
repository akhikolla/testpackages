/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 * 
 */
vector<mdreal>
BaseGaussian::parameters() const {
  vector<mdreal> output(7, 0.0);
  if(method == "exp") output[0] = 'E';
  if(method == "log") output[0] = 'L';
  if(method == "linear") output[0] = 'l';
  output[1] = center;
  output[2] = offset;
  output[3] = scale;
  output[4] = factor;
  output[5] = mu;
  output[6] = sigma;
  return output;
}

/*
 * 
 */
void
BaseGaussian::parameters(const vector<mdreal>& prm) {
  vector<mdreal> prmcopy = prm;
  prmcopy.resize(7, 0.0);
  (this->method).clear();
  if(prmcopy[0] == 'E') this->method = "exp";
  if(prmcopy[0] == 'L') this->method = "log";
  if(prmcopy[0] == 'l') this->method = "linear";
  this->center = prmcopy[1];
  this->offset = prmcopy[2];
  this->scale = prmcopy[3];
  this->factor = prmcopy[4];
  this->mu = prmcopy[5];
  this->sigma = prmcopy[6];
}
