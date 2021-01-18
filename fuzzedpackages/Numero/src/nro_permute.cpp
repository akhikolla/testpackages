/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

#define nro_NMAX_FALSE 20

/*
 *
 */
static string
nro_permute_exec(vector<mdreal>& output, const vector<mdreal>& values,
		 const vector<mdsize>& bmus, const Topology& topo,
		 const mdsize ncycl) {
  output.clear();
  if(values.size() != bmus.size())
    medusa::panic("Incompatible inputs.", __FILE__, __LINE__);
  
  /* Create a simulation engine. R-style indexing switched to C++. */
  Engine eng(topo);
  mdsize ndata = 0;
  for(mdsize i = 0; i < values.size(); i++) {
    string e = eng.insert(long2string(i), (bmus[i] - 1), values[i]);
    ndata += (e.size() < 1); /* excluded if any missing value */
  }
  
  /* Observed component plane. */
  vector<vector<mdreal> > observed = eng.average();
  mdsize nvars = observed.size();
  if(nvars != 1) return "Observation failed.";

  /* Observed statistic. */
  mdreal score = abacus::statistic(observed[0], "sd");

  /* Simulate null distributions. */    
  mdsize nfalse = 0;
  mdsize ntotal = 0;
  abacus::Empirical fnull;  
  for(mdsize n = 0; n < ncycl; n++) {
    if(nfalse/nvars >= nro_NMAX_FALSE) break;
    
    /* Permute without replacement. */
    eng.shuffle(false);
    
    /* Permuted component planes. */
    vector<vector<mdreal> > permuted = eng.average();
    if(permuted.size() != observed.size())
      return "Permutation failed.";
    
    /* Permuted statistic. */
    mdreal x = abacus::statistic(permuted[0], "sd");
    
    /* Check if any false positives. */
    nfalse += (x >= score);
    
    /* Update null distributions. */
    fnull.add(x, 1.0);
    ntotal++;
  }

  /* Estimate statistics. */
  output.push_back(score);
  output.push_back(fnull.z(score));
  output.push_back(fnull.p(score, 1));
  output.push_back(ndata);
  output.push_back(fnull.size());
  return "";
}

/*
 *
 */
RcppExport SEXP
nro_permute(SEXP topo_R, SEXP sigma_R, SEXP bmus_R, SEXP data_R,
	    SEXP numcycl_R, SEXP lag_R) {
  vector<mdsize> numcycl = nro::vector2sizes(numcycl_R);
  mdreal sigma = as<mdreal>(sigma_R);
  mdreal lag = as<mdreal>(lag_R);
  time_t stamp = time(NULL);
  
  /* Check inputs. */
  vector<mdsize> bmus = nro::vector2sizes(bmus_R);
  vector<vector<mdreal> > vectors = nro::matrix2reals(data_R, 0.0);
  if(vectors.size() < 1) return CharacterVector("No data.");
  if(vectors.size() != bmus.size())
    return CharacterVector("Incompatible inputs.");
  if(vectors[0].size() != numcycl.size())
    return CharacterVector("Incompatible inputs.");
  
  /* Get map topology. */
  vector<vector<mdreal> > topodata = nro::matrix2reals(topo_R, 0.0);
  punos::Topology topo = reals2topology(topodata, sigma);
  if(topo.size() < 1) return CharacterVector("Unusable topology.");

  /* Asynchronous permutations. That is, separate synchronous
     permutations for each single dimension. */
  time_t reset = stamp;
  mdsize nvars = vectors[0].size();
  vector<vector<mdreal> > stats;
  for(mdsize j = 0; j < nvars; j++) {
    
    /* Extract data values. */
    vector<mdreal> column;
    for(mdsize i = 0; i < vectors.size(); i++)
      column.push_back(vectors[i][j]);
    
    /* Estimate statistics. */
    vector<mdreal> batch;
    string err = nro_permute_exec(batch, column, bmus, topo, numcycl[j]);
    if(err.size() > 0) return CharacterVector(err);
    
    /* Update results. */
    stats.resize(batch.size());
    for(mdsize i = 0; i < batch.size(); i++)
      stats[i].push_back(batch[i]);
    
    /* Progress message. */
    if(lag < 0.0) continue;
    if(difftime(time(NULL), reset) < lag) continue;
    string dt = medusa::time2text(difftime(time(NULL), stamp));
    Rprintf("%.1f%% in %s\n", 100*(j + 1.0)/nvars, dt.c_str());
    reset = time(NULL);
  }

  /* Return results. */
  List res;
  res.push_back(nro::reals2vector(stats[0]), "SCORE");
  res.push_back(nro::reals2vector(stats[1]), "Z");
  res.push_back(nro::reals2vector(stats[2]), "P.freq");
  res.push_back(nro::reals2vector(stats[3]), "N.data");
  res.push_back(nro::reals2vector(stats[4]), "N.cycles");
  return res;
}
