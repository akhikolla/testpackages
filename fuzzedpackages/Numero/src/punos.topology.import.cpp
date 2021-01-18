/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

static string load_levels(TopologyBuffer*, File&);
static string load_coord(TopologyBuffer*, File&);
static string load_sigma(TopologyBuffer*, File&);

/*
 *
 */
string
Topology::import(const string& fname) {

  /* Discard previous contents. */
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  delete p; p = new TopologyBuffer();
  this->buffer = p;

  /* Open file. */
  File f; f.open(fname, "r");
  if(f.error().size() > 0) return f.error();

  /* Import vertical levels. */
  string err = load_levels(p, f);
  if(err.size() > 0) return err;

  /* Import node coordinates. */
  err = load_coord(p, f);
  if(err.size() > 0) return err;

  /* Import neighborhood radius. */
  err = load_sigma(p, f);
  if(err.size() > 0) return err;

  /* Set neighborhood radius. */
  if(this->rewire(p->sigma)) return "";
  return "Unusable neighborhood radius.";
}

/*
 *
 */
string
load_levels(TopologyBuffer* p, File& f) {
  mdreal rlnan = medusa::rnan();

  /* Find first row. */
  while(f.error().size() < 1) {
    vector<string> fields = f.read('\t', 1);
    if(fields[0] == "LEVEL") break;
  }

  /* Import data rows. */
  vector<mdreal> levels;
  while(f.error().size() < 1) {
    vector<string> fields = f.read('\t', 0);
    if(fields.size() != 1) break;
    mdreal z = string2real(fields[0]);
    if(z == rlnan) return "Unusable level.";
    levels.push_back(z);
  }

  /* Check if enough data. */
  mdsize nlevels = levels.size();
  if(nlevels < 1) return "No levels.";
  if(nlevels != uniqreal(levels).size())
    return "Duplicate levels.";

  /* Make sure levels are sorted. */
  sort(levels.begin(), levels.end());
  p->levels = levels;
  return "";
}

/*
 *
 */
string
load_coord(TopologyBuffer* p, File& f) {
  mdreal rlnan = medusa::rnan();

  /* Find first row. */
  while(f.error().size() < 1) {
    vector<string> fields = f.read('\t', 1);
    if(fields[0] == "DISTRICT") break;
  }

  /* Import data rows. */
  vector<District> districts; District u;
  while(f.error().size() < 1) {
    vector<string> fields = f.read('\t', 0);
    if(fields.size() != 7) break;

    /* Convert to real values. */
    vector<mdreal> values(fields.size());
    for(mdsize j = 0; j < fields.size(); j++) {
      values[j] = string2real(fields[j]);
      if(values[j] == rlnan) return "Unusable value.";
    }

    /* Check district identity. */
    mdsize key = (mdsize)(values[0] + 0.5);
    if(key != districts.size()) return "Inconsistent district data.";
    
    /* Store district. */
    u.x = values[1];
    u.y = values[2];
    u.radii.first = values[3];
    u.radii.second = values[4];
    u.angles.first = values[5];
    u.angles.second = values[6];
    districts.push_back(u);
  }

  /* Check if any data. */
  if(districts.size() < 1) return "No districts.";

  /* Determine map radius. */
  mdreal rmax = rlnan;
  for(mdsize i = 0; i < districts.size(); i++) {
    if(rmax == rlnan) rmax = districts[i].radii.second;
    if(districts[i].radii.second < rmax) continue;
    rmax = districts[i].radii.second;
  }
  
  /* Update object. */
  p->maxradius = rmax;
  p->coord = districts;
  return "";
}

/*
 *
 */
string
load_sigma(TopologyBuffer* p, File& f) {
  mdreal rlnan = medusa::rnan();
  while(f.error().size() < 1) {
    vector<string> fields = f.read('\t', 1);
    if(fields[0] != "SIGMA") continue;
    fields = f.read('\t', 0);
    if(fields.size() != 1) return "Incompatible file.";
    mdreal sigma = string2real(fields[0]);
    if(sigma == rlnan) return "Unusable neighborhood radius.";
    if(sigma <= 0.0) return "Non-positive neighborhood radius.";
    p->sigma = sigma;
    return "";
  }
  return "No neighborhood radius.";
}
