/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
class GaussianMinimizer : public Minimizer {
public:
  mdreal mu;
  mdreal sigma;
  mdreal delta;
  vector<mdreal>* values;
  vector<mdreal>* weights;
  Gaussian* source;
public:
  GaussianMinimizer() : Minimizer() {
    this->mu = 0.0;
    this->sigma = 1.0;
    this->delta = -1.0;
  };
  ~GaussianMinimizer() {};
  mdreal value(const mdreal f) {
    const Gaussian& gauss = *source;
    vector<mdreal> x = *values; /* deep copy */

    /* Apply transform. */
    gauss.apply(x, f);
    
    /* Weighted statistics. */
    mdreal xmu = statistic(x, *weights, "mean");
    mdreal xsigma = statistic(x, *weights, "sd");

    /* Safeguard against zero variance. */
    if(xsigma < 1e-9) xsigma = 1e-9;
    
    /* Distance from Gaussian form. */
    mdreal d = gauss.distance(f, xmu, xsigma);

    /* Check if improved. */
    if((delta < 0.0) || (d < delta)) {
      this->mu = xmu;
      this->sigma = xsigma;
      this->delta = d;
    }
    return d;
  };
};

/*
 *
 */
mdreal
Gaussian::optimize(const string& mcode) {
  mdreal rlnan = medusa::rnan();
  vector<mdreal> x = values;
  vector<mdreal> w = weights;
  if(center == rlnan) return rlnan;

  /* Linear transformation. */
  if(mcode == "linear") {
    this->method = mcode;
    this->factor = 0.0;
    this->mu = statistic(x, w, "mean");
    this->sigma = statistic(x, w, "sd");
    return this->quality();
  }

  /* Non-linear transformation. */
  if((mcode == "exp") || (mcode == "log")) {
    this->method = mcode;
    
    /* Prepare optimization engine. */
    GaussianMinimizer engine;
    engine.algorithm(8, 1e-6);
    engine.space(0.0, 1.0);
    
    /* Link data resources. */
    engine.values = &values;
    engine.weights = &weights;
    engine.source = this;

    /* Find optimal parameter setting. */
    this->factor = Minimizer::optimize(engine);
    this->mu = engine.mu;
    this->sigma = engine.sigma;
    return this->quality();
  }
  panic("Unknown method.", __FILE__, __LINE__);
  return 0.0;
}
