/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
mdsize
Frame::group() {
  FrameBuffer* p = (FrameBuffer*)buffer;
  if(p->ngroups > 0) { 
    p->append("</g>\n");
    p->ngroups -= 1;
  }
  return p->ngroups;
}

/*
 *
 */
mdsize
Frame::group(const string& key) {
  FrameBuffer* p = (FrameBuffer*)buffer;  
  if(key.size() < 1) p->append("\n<g>");
  else p->append("\n<g id=\"" + key + "\">");
  p->ngroups += 1;
  return p->ngroups;
}
