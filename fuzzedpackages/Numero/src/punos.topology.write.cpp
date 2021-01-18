/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
Frame
Topology::write(const mdreal xorig, const mdreal yorig,
		const vector<string>& labels,
		const vector<Color>& colors,
		const Style& base) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  mdreal rho = TopologyBuffer::scale();

  /* Check input. */
  const vector<District>& districts = p->coord;
  if((labels.size() != districts.size()) ||
     (colors.size() != districts.size())) {
    medusa::worry("Incompatible inputs.", __FILE__);
    return Frame();
  }

  /* Force labels to the center of district areas. */
  Style sty = base;
  sty.pointable = false;
  sty.anchor = "middle";

  /* Open group. */
  Frame fr;
  fr.group(base.identity);

  /* Write labels. */
  for(mdsize k = 0; k < labels.size(); k++) {
    const string& label = labels[k];
    const Color& c = colors[k];
    if(label.size() < 1) continue;
    if(c.opacity <= 0.0) continue;

    /* Scale positions to canvas coordinates. */
    const District& u = districts[k];
    mdreal x = (xorig + rho*(u.x));
    mdreal y = (yorig + rho*(u.y));

    /* Fine tune label position. */
    if(label[0] == '+') x -= 0.33*(sty.fontsize);
    if(label[0] == '-') x -= 0.33*(sty.fontsize);

    /* Set label color and identity. */
    sty.fillcolor = c;
    sty.strokecolor = c;
    if(base.identity.size() > 0)
      sty.identity = (base.identity + "_" + long2string(k));
    fr.stylize(sty);

    /* Write label. */
    if(!fr.text(x, y, label)) return Frame();
  }

  /* Close group. */
  fr.group();
  return fr;
}
