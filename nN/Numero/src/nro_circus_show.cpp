/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_circus_show(SEXP offsets_R, SEXP topo_R, SEXP ccodes_R,
		SEXP labels_R, SEXP key_R) {
  string key = as<string>(key_R);
  key = string2safe(key, key.size());
  
  /* Position offsets. */
  vector<mdreal> xy = nro::vector2reals(offsets_R);
  xy.resize(2, 0.0);

  /* Check inputs. */
  vector<string> ccodes = as<vector<string> >(ccodes_R);
  vector<string> labtxt = as<vector<string> >(labels_R);
  if(ccodes.size() != labtxt.size())
    return CharacterVector("Incompatible inputs.");
  
  /* Get map topology. */
  vector<vector<mdreal> > topodata = nro::matrix2reals(topo_R, 0.0);
  punos::Topology topo = reals2topology(topodata, 0.0);
  if(topo.size() < 1) return CharacterVector("Unusable topology.");

  /* Convert color codes to color objects. */
  vector<Color> colors;
  for(mdsize i = 0; i < ccodes.size(); i++)
    colors.push_back(Color(ccodes[i]));
  
  /* Make sure all labels are single characters. */
  vector<char> labels;
  for(mdsize i = 0; i < labtxt.size(); i++) {
    string c = medusa::string2safe(labtxt[i], 1);
    if(c.size() < 1) labels.push_back('\0');
    else labels.push_back(c[0]);
  }
  
  /* Write subgroup labels. */
  scriptum::Style sty;
  sty.identity = (key + "_region");
  sty.fontsize *= 0.8;
  sty.fillcolor = scriptum::colormap(1.0, "grey");
  sty.strokecolor = scriptum::colormap(1.0, "grey");
  scriptum::Frame frame = topo.highlight(xy[0], xy[1], labels,
					 colors, sty);

  /* Get bounding box. */
  vector<mdreal> bbox(4);
  bbox[0] = frame.horizontal().first;
  bbox[1] = frame.vertical().first;
  bbox[2] = frame.horizontal().second;
  bbox[3] = frame.vertical().second;
  
  /* Return results. */
  List res;
  res.push_back(frame.flush(), "code");
  res.push_back(bbox, "bbox");
  return res;
}
