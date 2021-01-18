/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_webpage(SEXP fname_R, SEXP bytes_R) {
  string fname = as<string>(fname_R);
  vector<string> bytes = as<vector<string> >(bytes_R);
  
  /* Open output file. */
  medusa::File f; f.open(fname, "w");

  /* Save data. */
  for(mdsize i = 0; i < bytes.size(); i++) {
    if(f.write(bytes[i]) > 0) continue;
    return CharacterVector(f.error());
  }
  
  /* Return file size. */
  List output;
  output.push_back(long2string(f.size()), "nbytes");
  output.push_back(long2text(f.size()), "text");
  return output;
}
