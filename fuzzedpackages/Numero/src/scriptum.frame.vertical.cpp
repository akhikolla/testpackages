/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
pair<mdreal, mdreal>
Frame::vertical() const {
  FrameBuffer* p = (FrameBuffer*)buffer;
  mdreal rlnan = medusa::rnan();
  pair<mdreal, mdreal> res;
  res.first = (p->limits).second.alpha;
  res.second = (p->limits).second.omega;
  if((res.first == rlnan) || (res.second == rlnan)) {
    res.first = 0.0;
    res.second = 0.0;
  }
  return res;
}
