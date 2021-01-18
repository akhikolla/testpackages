/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
bool
Frame::curve(const vector<mdreal>& vx, const vector<mdreal>& vy) {
  mdreal rlnan = medusa::rnan();
  FrameBuffer* p = (FrameBuffer*)buffer;

  /* Check size. */
  mdsize nv = vx.size();
  if(nv < 2) return false;
  if(vy.size() != nv) return false;

  /* Check if closed path. */
  bool closeflag = ((vx[0] == vx[nv-1]) && (vy[0] == vy[nv-1]));
  if(closeflag && (nv < 4)) return false;
  if(closeflag) nv--;

  /* Check if unusable coordinates. */
  for(mdsize i = 0; i < nv; i++) {
    if(vx[i] == rlnan) return false;
    if(vy[i] == rlnan) return false;
  }

  /* Create polyline. */
  sprintf(p->f(), "\n<path d=\"\n");
  sprintf(p->f(), "M\t%.2f\t%.2f", vx[0], vy[0]);
  for(mdsize i = 1; i < nv; i++)
    sprintf(p->f(), "\nL\t%.2f\t%.2f", vx[i], vy[i]);
  
  /* Set closing flag. */
  if(closeflag) p->append("\nZ");
  p->append("\"\n");
  
  /* Apply style. */
  p->append(p->linestycode);
  p->append("/>\n");

  /* Update limits and filesize. */
  (p->limits).first.update(vx, p->style);
  (p->limits).second.update(vy, p->style);
  return true;
}

/*
 *
 */
bool
Frame::curve(const medusa::mdreal& xA, const medusa::mdreal& yA,
	     const medusa::mdreal& x0, const medusa::mdreal& y0,
	     const medusa::mdreal& xB, const medusa::mdreal& yB) {
  mdreal rlnan = medusa::rnan();
  FrameBuffer* p = (FrameBuffer*)buffer;
  
  /* Check inputs. */
  if(xA == rlnan) return false;
  if(yA == rlnan) return false;
  if(x0 == rlnan) return false;
  if(y0 == rlnan) return false;
  if(xB == rlnan) return false;
  if(yB == rlnan) return false;

  /* Create polyline. */
  sprintf(p->f(), "\n<path d=\"\n");
  sprintf(p->f(), "M\t%.2f\t%.2f", xA, yA);
  sprintf(p->f(), "\nQ\t%.2f\t%.2f", x0, y0);
  sprintf(p->f(), "\n\t%.2f\t%.2f\"\n", xB, yB);
  
  /* Apply style. */
  p->append(p->linestycode);
  p->append("/>\n");

  /* Update limits and filesize. */
  pair<Limes, Limes>& lims = p->limits;
  lims.first.update(xA, p->style);
  lims.first.update(x0, p->style);
  lims.first.update(xB, p->style);
  lims.second.update(yA, p->style);
  lims.second.update(y0, p->style);
  lims.second.update(yB, p->style);
  return true;
}

