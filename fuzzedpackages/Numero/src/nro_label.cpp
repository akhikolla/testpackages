/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_label(SEXP topo_R, SEXP data_R, SEXP binflags_R, SEXP sigma_R) {
  mdreal rlnan = medusa::rnan();
  mdreal sigma = as<mdreal>(sigma_R);

  /* Get map topology.*/
  vector<vector<mdreal> > topodata = nro::matrix2reals(topo_R, 0.0);
  punos::Topology topo = nro::reals2topology(topodata, 0.0);
  if(topo.size() < 1)
    return CharacterVector("Unusable topology.");

  /* Set neighborhood radius. */
  if(topo.rewire(sigma) == false)
    return CharacterVector("Unusable label gap.");

  /* Check data vectors. */
  vector<vector<mdreal> > vectors = nro::matrix2reals(data_R, 0.0);
  if(vectors.size() < 0) return CharacterVector("No data.");
  if(vectors.size() != topo.size())
    return CharacterVector("Incompatible inputs.");

  /* Transpose to column-major. */
  vector<vector<mdreal> > columns(vectors[0].size());
  for(mdsize i = 0; i < vectors.size(); i++) {
    for(mdsize j = 0; j < vectors[i].size(); j++)
      columns[j].push_back(vectors[i][j]);
    vectors[i].clear();
  }

  /* Indicators for binary variables. */
  IntegerVector binflags(binflags_R);

  /* Process columns. */
  List labels;
  List visible;
  for(mdsize j = 0; j < columns.size(); j++) {
    vector<mdreal>& values = columns[j];
    
    /* Calculate deviations from mean. */
    vector<mdreal> amps = abacus::transform(values, "z");
    for(mdsize i = 0; i < amps.size(); i++)
      if(amps[i] != rlnan) amps[i] = fabs(amps[i]);
    
    /* Sort map units according to amplitude. */
    vector<mdreal> tmp = amps;
    vector<mdsize> sorted = sortreal(tmp, -1);
    
    /* Select labeled units. */
    vector<bool> flags(amps.size(), true);
    for(mdsize k = 0; k < amps.size(); k++) {

      /* Skip units that have already been denied. */
      mdsize unit = sorted[k];
      if(!flags[unit]) continue;

      /* Deny label from too close neighbors. */
      vector<mdsize> neigh = topo.neighbors(unit);
      for(vector<mdsize>::iterator pos = neigh.begin();
	  pos != neigh.end(); pos++) {
	if(*pos == unit) continue;
	mdreal r = topo.distance(unit, *pos);
	if(r > sigma) continue;
	flags[*pos] = false;
      }
    }

    /* Check if only positive values. */
    mdreal xmin = abacus::statistic(values, "min");
    mdreal xmax = abacus::statistic(values, "max");

    /* Check if binary. */
    bool binary = (binflags[j] > 0);
    if(xmin < 0.0) binary = false;
    if(xmax > 1.0) binary = false;

    /* Parse label text. */
    vector<string> parsed(amps.size());
    for(mdsize k = 0; k < amps.size(); k++) {
      string& txt = parsed[k];
      if(!binary)
	txt = medusa::real2text(values[k]);
      else {
	mdreal percent = 100*(values[k]);
	if(percent < 0.1) {
	  if(percent < 0.05) txt = "0.0%";
	  else txt = "0.1%";
	}
	else {
	  txt = (medusa::real2text(percent) + "%");
	}
      }
      if((xmin >= 0) && (txt[0] == '+')) txt = txt.substr(1);
    }

    /* Save labels. */
    string name = ("Col" + long2string(j));
    labels.push_back(parsed, name);
    visible.push_back(flags, name);
  }

  /* Return results. */
  List output;
  output.push_back(labels, "labels");
  output.push_back(visible, "visible");
  return output;
}
