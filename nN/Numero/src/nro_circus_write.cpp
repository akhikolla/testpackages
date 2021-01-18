/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_circus_write(SEXP offsets_R, SEXP topo_R, SEXP labels_R, SEXP visible_R,
		 SEXP contrast_R, SEXP key_R, SEXP font_R) {
  scriptum::Color black = scriptum::colormap(0.0, "grey");
  scriptum::Color grey = scriptum::colormap(0.6, "grey");
  scriptum::Color white = scriptum::colormap(1.0, "grey");

  /* Check subplot key. */
  string key = as<string>(key_R);
  key = string2safe(key, key.size());

  /* Check font size. */
  mdreal fntsize = as<mdreal>(font_R);
  if(fntsize < 0.1) return CharacterVector("Unusable font.");
  if(fntsize > 100) return CharacterVector("Unusable font.");
  
  /* Position offsets. */
  vector<mdreal> xy = nro::vector2reals(offsets_R);
  xy.resize(2, 0.0);

  /* Check inputs. */
  vector<string> labels = as<vector<string> >(labels_R);
  vector<mdsize> visible = nro::vector2sizes(visible_R);
  vector<mdsize> contrast = nro::vector2sizes(contrast_R);
  if(visible.size() != labels.size())
    return CharacterVector("Incompatible inputs.");
  if(contrast.size() != labels.size())
    return CharacterVector("Incompatible inputs.");
  
  /* Get map topology. */
  vector<vector<mdreal> > topodata = nro::matrix2reals(topo_R, 0.0);
  punos::Topology topo = reals2topology(topodata, 0.0);
  if(topo.size() < 1) return CharacterVector("Unusable topology.");

  /* Make sure all text elements are safe. */
  for(mdsize i = 0; i < labels.size(); i++)
    labels[i] = medusa::string2safe(labels[i], 16);
  
  /* Set labels and colors according to visibility and contrast. */
  vector<scriptum::Color> colors(labels.size());
  vector<scriptum::Color> shadows(labels.size());  
  for(mdsize i = 0; i < labels.size(); i++) {
    if(visible[i] == 0) {
      labels[i].clear();
      continue;
    }
    if(contrast[i] != 0)
      colors[i] = black;
    else {
      colors[i] = white;
      shadows[i] = grey;
    }
  }

  /* Write label shadows. */
  scriptum::Style sty;
  sty.identity = (key + "_shadow");
  sty.fontsize *= fntsize;
  sty.strokewidth = 0.16*(sty.fontsize);
  scriptum::Frame frameB = topo.write(xy[0], xy[1], labels, shadows, sty);
 
  /* Write labels. */
  sty = Style();
  sty.identity = (key + "_label");
  sty.fontsize *= fntsize;
  sty.strokewidth = 0.0;
  scriptum::Frame frameT = topo.write(xy[0], xy[1], labels, colors, sty);

  /* Get bounding box. */
  vector<mdreal> bboxB(4);
  vector<mdreal> bboxT(4);
  bboxB[0] = frameB.horizontal().first;
  bboxB[1] = frameB.vertical().first;
  bboxB[2] = frameB.horizontal().second;
  bboxB[3] = frameB.vertical().second;
  bboxT[0] = frameT.horizontal().first;
  bboxT[1] = frameT.vertical().first;
  bboxT[2] = frameT.horizontal().second;
  bboxT[3] = frameT.vertical().second;
  
  /* Return results. */
  List res;
  res.push_back(frameB.flush(), "code.shadow");
  res.push_back(frameT.flush(), "code.label");
  res.push_back(bboxB, "bbox.shadow");
  res.push_back(bboxT, "bbox.label");
  return res;
}
