/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
string
Model::train(vector<Resident>& layout, vector<mdreal>& trace,
	     const mdreal quota) {
  mdreal rlnan = medusa::rnan();
  ModelBuffer* p = (ModelBuffer*)buffer;
  map<string, Point>& points = p->points;
  mt19937& twister = p->twister;
  time_t stamp = time(NULL);
  
  /* Clear output containers. */
  layout.clear();
  trace.clear();

  /* Check resources. */
  Topology topocopy = p->topology;
  mdsize npoints = (p->points).size();
  if(topocopy.size() < 1) return "No map units.";
  if(npoints < 10) return "Too few points.";
  if(p->ntrain < 10) return "Too few training points.";

  /* Final neighborhood radius for trained map. */
  mdreal sigma = (p->topology).sigma();

  /* Set initial neighborhood radius. */
  vector<mdreal>& history = p->history;
  map<string, mdreal>& state = p->state;
  if(state.size() < 1) {
    state["rho"] = 0.5*(topocopy.radius());
    if(state["rho"] < sigma) state["rho"] = sigma;
    history.clear();
  }

  /* Create training engine. */
  Trainer trainer(p->codebook, topocopy, p->ntrain, p->equality);
  
  /* Make pointers to points. */
  vector<Point*> pointers;
  for(map<string, Point>::iterator it = points.begin();
      it != points.end(); it++)
    pointers.push_back(&(it->second));

  /* Prepare sampling mask. */
  vector<Point*> mask = pointers;
  if(p->ntrain < npoints) mask.resize(p->ntrain);
  
  /* Fit codebook to training data. */
  while(state["rho"] >= 0.0) {
    
    /* Set neighborhood radius. */
    topocopy.rewire(state["rho"]);

    /* Run training batch. */
    bool finished = false;
    while(!finished) {
      
      /* Shuffle pointers and sampling mask. */
      if(mask.size() < npoints) {
	for(mdsize i = 0; i < mask.size(); i++) {
	  mdsize rank = twister()%npoints;
	  Point* pnt = pointers[rank];
	  pointers[rank] = pointers[i];
	  pointers[i] = pnt;
	}	
	for(mdsize i = 0; i < mask.size(); i++)
	  mask[i] = pointers[i];
      }
 
      /* Perform a training cycle. */
      mdreal delta = trainer.cycle(mask, topocopy);

      /* Check if initial centroids were available. */
      if(delta == rlnan) {
	if(history.size() > 0) return "Training cycle failed.";
	if(trace.size() < 1) delta = trainer.cycle(mask, topocopy);
      }

      /* Store training error. */
      trace.push_back(delta);
      history.push_back(delta);
      finished = convergence(history, 0.01);

      /* Check time quota. */
      if(quota == rlnan) continue;
      if(difftime(time(NULL), stamp) >= quota) break; 
    }

    /* Check if batch was finished. */
    if(!finished) break;
    history.clear();
 
    /* Check if training is finished. */
    if(state["rho"] <= sigma) {
      state["rho"] = -1.0;
      break;
    }

    /* Update neighborhood radius. */
    state["rho"] *= 0.67;
    if(state["rho"] < sigma) state["rho"] = sigma;
  }

  /* Update codebook. */
  p->codebook = trainer.codebook();
  
  /* Return final layout. */
  for(map<string, Point>::iterator it = points.begin();
      it != points.end(); it++) {
    Resident res;
    res.identity = it->first;
    res.district = (it->second).location();
    res.residual = trainer.distance(it->second, res.district);
    layout.push_back(res);
  }
  return "";
}
