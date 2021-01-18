/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
Frame
Topology::highlight(const mdreal xorig, const mdreal yorig,
		    const vector<char>& labels,
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
  
  /* Set style parameters. */
  Style fgsty = base; fgsty.identity.clear();
  Style bgsty = base; bgsty.identity.clear();
  fgsty.anchor = "middle";
  bgsty.strokewidth = 0.1*(base.fontsize);
  fgsty.pointable = false;
  bgsty.pointable = false;

  /* Open group. */
  Frame fr;
  fr.group(base.identity);
  
  /* Write highlighted labels. */
  for(mdsize k = 0; k < labels.size(); k++) {
    char label = labels[k];
    if(label == '\0') continue;

    /* Open highlight group. */
    fr.group(base.identity + "_" + long2string(k));

    /* Scale positions to canvas coordinates. */
    const District& u = districts[k];
    mdreal x = (xorig + rho*(u.x));
    mdreal y = (yorig + rho*(u.y));
    
    /* Set background color. */
    bgsty.fillcolor = colors[k];
    fr.stylize(bgsty);
    
    /* Draw background circle. */
    if(!fr.shape(x, y, 0.8*(fgsty.fontsize), "circle")) return Frame();
    
    /* Fine-tune label position. */
    y += 0.05*(fgsty.fontsize);
    
    /* Set label style. */
    fr.stylize(fgsty);
    
    /* Write label. */
    string txt(1, label);
    if(!fr.text(x, y, txt)) return Frame();

    /* Close highlight group. */
    fr.group();
  }

  /* Close group. */
  fr.group();
  return fr;
}
