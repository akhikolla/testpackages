/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_circus_paint(SEXP offsets_R, SEXP topo_R, SEXP ccodes_R,
		 SEXP key_R, SEXP title_R) {
  string key = as<string>(key_R);
  string title = as<string>(title_R);
  
  /* Position offsets. */
  vector<mdreal> xy = nro::vector2reals(offsets_R);
  xy.resize(2, 0.0);
  
  /* Get map topology. */
  vector<vector<mdreal> > topodata = nro::matrix2reals(topo_R, 0.0);
  punos::Topology topo = reals2topology(topodata, 0.0);
  if(topo.size() < 1) return CharacterVector("Unusable topology.");

  /* Make sure all text elements are safe. */
  mdsize cap = (mdsize)(6*(topo.radius() + 1.0));
  if(key.size() >= cap) return CharacterVector("Unusable identifier.");
  if(key.size() > 8) return CharacterVector("Unusable identifier.");
  key = string2safe(key, (cap - 1));
  title = string2safe(title, (cap - key.size()));

  /* Convert color codes to color objects. */
  vector<Color> colors;
  vector<string> ccodes = as<vector<string> >(ccodes_R);
  for(mdsize i = 0; i < ccodes.size(); i++)
    colors.push_back(Color(ccodes[i]));

  /* Paint component plane. */
  scriptum::Style sty;
  sty.identity = key;
  sty.strokewidth = 0.0;
  scriptum::Frame frame = topo.paint(xy[0], xy[1], colors, sty);

  /* Get bounding box of the colored area. */
  vector<mdreal> inner(4);
  inner[0] = frame.horizontal().first;
  inner[1] = frame.vertical().first;
  inner[2] = frame.horizontal().second;
  inner[3] = frame.vertical().second;
  
  /* Full bounding box. */
  vector<mdreal> outer = inner;
  outer[0] = ceil(outer[0] - 1.0*(sty.fontsize) - 0.5);
  outer[1] = ceil(outer[1] - 1.9*(sty.fontsize) - 0.5);
  outer[2] = floor(outer[2] + 1.0*(sty.fontsize) + 0.5);
  outer[3] = floor(outer[3] + 1.3*(sty.fontsize) + 0.5);

  /* Determine title location. */
  mdreal x = (inner[0] + 1.4*(sty.fontsize));
  mdreal y = (inner[1] - 0.9*(sty.fontsize));
  for(mdsize k = 1; k < key.size(); k++)
    x += 0.67*(sty.fontsize);

  /* Set style attributes for title. */
  sty = Style(); sty.identity = (key + "_title");
  sty.fillcolor = scriptum::colormap(0.0, "grey");
  sty.anchor = "start";
  sty.strokewidth = 0.0;
  frame.stylize(sty);

  /* Write title. */
  if((title.size() > 0) && (frame.text(x, y, title) == false))
    return CharacterVector("Title failed.");

  /* Plot identifier location. */
  x = inner[0];

  /* Set style attributes for identifier. */
  sty.identity = (key + "_key");
  sty.fontweight = 900;
  sty.anchor = "start";
  frame.stylize(sty);
  
  /* Write identifier. */
  if((title.size() > 0) && (frame.text(x, y, key) == false))
    return CharacterVector("Identifier failed.");
  
  /* Return results. */
  List res;
  res.push_back(key, "key");
  res.push_back(frame.flush(), "code");
  res.push_back(inner, "bbox.inner");
  res.push_back(outer, "bbox");
  return res;
}
