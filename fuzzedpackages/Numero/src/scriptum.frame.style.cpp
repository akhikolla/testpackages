/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
Style
Frame::style() const {
  FrameBuffer* p = (FrameBuffer*)buffer;
  return p->style;
}
