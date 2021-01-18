/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_colorize(SEXP zvals_R, SEXP name_R) {
  scriptum::Color black = scriptum::colormap(0.0, "grey");
  scriptum::Color white = scriptum::colormap(1.0, "grey");
  string name = as<string>(name_R);
  mdreal rlnan = medusa::rnan();

  /* Import color scores. */
  vector<vector<mdreal> > zvals = nro::matrix2reals(zvals_R, 0.0);
  if(zvals.size() < 1) return CharacterVector("Empty input.");

  /* Set colors. */
  List flags;
  List colors;
  for(mdsize j = 0; j < zvals[0].size(); j++) {
    vector<bool> bits;
    vector<string> array;
    for(mdsize i = 0; i < zvals.size(); i++) {
      mdreal z = zvals[i][j];
      if(z == rlnan) z = 0.5;
      Color c = scriptum::colormap(z, name);
      mdreal cB = black.contrast(c);
      mdreal cW = white.contrast(c);
      bits.push_back(fabs(cB) > fabs(cW));
      array.push_back("#" + c.hex());
    }
    flags.push_back(bits);
    colors.push_back(array);
  }

  /* Return results. */
  List output;
  output.push_back(flags, "contrast");
  output.push_back(colors, "colors");
  return output;
}
