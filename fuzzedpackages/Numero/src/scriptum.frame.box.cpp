/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
bool
Frame::box(const mdreal x, const mdreal y) {
  FrameBuffer* p = (FrameBuffer*)buffer;
  pair<Limes, Limes>& lims = p->limits;
  Style sty; sty.padding = 0.0;
  if(!lims.first.update(x, sty)) return false;
  if(!lims.second.update(y, sty)) return false;
  return true;
}
