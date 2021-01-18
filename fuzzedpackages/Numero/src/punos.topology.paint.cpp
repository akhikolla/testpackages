/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
Frame
Topology::paint(const mdreal xorig, const mdreal yorig,
		const vector<Color>& colors, const Style& base) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  char textbuf[256];
  
  /* Check input. */
  const vector<District>& districts = p->coord;
  if(colors.size() != districts.size()) {
    medusa::worry("Incompatible input.", __FILE__);
    return Frame();
  }
  
  /* Set base attributes. */
  Style sty = base;
  sty.pointable = true;
  sty.strokewidth = 0.0;
  
  /* Open group. */
  Frame fr;
  fr.group(base.identity + "_map");
  
  /* Unit square. */
  mdreal rho = TopologyBuffer::scale();
  vector<mdreal> xcalibr(5, xorig);
  vector<mdreal> ycalibr(5, yorig);
  xcalibr[1] += rho;
  xcalibr[2] += rho;
  ycalibr[2] += rho;
  ycalibr[3] += rho;

  /* Create the calibration element. */
  sty.values.resize(2);
  sty.values[0] = base.identity;
  sty.values[1] = long2string(districts.size());
  sty.identity = (base.identity + "_calibration");
  sty.fillcolor = Color("#ffffff");
  sty.pointable = false;
  fr.stylize(sty);
  if(!fr.curve(xcalibr, ycalibr)) return Frame();

  /* Check if anything else to do. */
  if(p->maxradius <= 0.0) {
    fr.group();
    return fr;
  }

  /* Open subgroup. */
  fr.group(base.identity + "_paint");

  /* Draw slices. */
  mdreal rmax = 0.0;
  for(mdsize i = 0; i < districts.size(); i++) {
    if(colors[i].opacity <= 0.0) continue;

    /* Scale positions to canvas coordinates. */
    const District& u = districts[i];
    mdreal r1 = rho*(u.radii.first);
    mdreal r2 = rho*(u.radii.second);
    mdreal a1 = u.angles.first;
    mdreal a2 = u.angles.second;
    
    /* Apply a small inflation to fill gaps. */
    if(r1 > 1e-9) {
      r1 -= 0.07;
      r2 += 0.07;
      a1 -= 0.07;
      a2 += 0.07;
    }
    
    /* Update circle radius. */
    if(r2 > rmax) rmax = r2;

    /* Unique identity. */
    string key = long2string(i);
    sty.identity = (base.identity + "_paint_" + key);

    /* Additional positional information. */
    sty.values.resize(4);
    sty.values[0] = base.identity;
    sty.values[1] = key;
    sprintf(textbuf, "%.4f", rho*(u.x));
    sty.values[2] = string(textbuf);
    sprintf(textbuf, "%.4f", rho*(u.y));
    sty.values[3] = string(textbuf);
    
    /* Set slice attributes. */
    sty.pointable = true;
    sty.fillcolor = colors[i];
    sty.strokewidth = 0.0;
    fr.stylize(sty);

    /* Draw slice. */
    if(!fr.slice(xorig, yorig, r1, r2, a1, a2)) return Frame();
  }
  
  /* Close subgroup. */
  fr.group();

  /* Set line style. */
  sty.pointable = false;
  sty.strokewidth = 0.5;
  sty.strokecolor = scriptum::colormap(0.6, "grey");
  sty.fillcolor.opacity = 0.0;
  if(base.identity.size() > 0)
    sty.identity = (base.identity + "_perimeter");
  sty.values.clear();
  fr.stylize(sty);

  /* Enclose in a circle. */
  if(!fr.shape(xorig, yorig, rmax, "circle")) return Frame();

  /* Close group. */
  fr.group();
  return fr;
}
