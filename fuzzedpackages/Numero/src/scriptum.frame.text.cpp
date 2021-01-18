/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
bool
Frame::text(const mdreal x, const mdreal y, const string& s) {
  FrameBuffer* p = (FrameBuffer*)buffer;
  Style& sty = p->style;
  
  /* Check if unusable coordinates. */
  if(x == medusa::rnan()) return false;
  if(y == medusa::rnan()) return false;
  if(s.size() < 1) return false;

  /* Create element. */
  mdreal fs = sty.fontsize;
  sprintf(p->f(), "\n<text x=\"%.3f\" ", x);
  sprintf(p->f(), "y=\"%.3f\"\n", (y + 0.34*fs));
  p->append(p->textstycode);
  p->append(">\n"); p->append(s); 
  p->append("\n</text>\n");
 
  /* Estimate horizontal bounds. */
  mdsize len = s.size();
  double width = 0.58*len*fs;
  vector<mdreal> xlims(2, x);
  if(sty.anchor == "middle") {
    xlims[0] -= 0.5*width;
    xlims[1] += 0.5*width;
  }
  if(sty.anchor == "end") xlims[0] -= width;
  if(xlims[0] == xlims[1]) xlims[1] += width;

  /* Estimate element height. */
  vector<mdreal> ylims(2, y);
  ylims[0] -= 0.5*fs;
  ylims[1] += 0.6*fs;

  /* Take rotation angle into account for bounding box. */
  if(sty.angle != 0.0) {
    vector<mdreal> origin = sty.origin;
    origin.resize(2, 0.0);

    /* Convert to polar coordinates. */
    pair<mdreal, mdreal> pmin;
    pair<mdreal, mdreal> pmax;
    pmin = polarize(origin[0], origin[1], xlims[0], ylims[0]);
    pmax = polarize(origin[0], origin[1], xlims[1], ylims[1]);

    /* Apply rotation. */
    pmin.second += M_PI*(sty.angle)/180.0;
    pmax.second += M_PI*(sty.angle)/180.0;
    
    /* Update bounding box. */
    xlims[0] = (x + (pmin.first)*cos(pmin.second));
    ylims[0] = (y + (pmin.first)*sin(pmin.second));
    xlims[1] = (x + (pmax.first)*cos(pmax.second));
    ylims[1] = (y + (pmax.first)*sin(pmax.second));
  }

  /* Update limits and filesize. */
  (p->limits).first.update(xlims, sty);
  (p->limits).second.update(ylims, sty);
  return true;
}
