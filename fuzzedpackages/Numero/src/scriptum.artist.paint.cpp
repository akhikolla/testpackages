/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
bool
Artist::paint(Frame& f) {
  ArtistBuffer* p = (ArtistBuffer*)buffer;
  if(p->output == NULL) return false;

  /* Update limits. */
  pair<mdreal, mdreal> xlim = f.horizontal();
  pair<mdreal, mdreal> ylim = f.vertical();
  (p->limits).first.update(xlim.first);
  (p->limits).first.update(xlim.second);
  (p->limits).second.update(ylim.first);
  (p->limits).second.update(ylim.second);

  /* Save data. */
  string data = f.flush();
  p->filesize += fprintf(p->output, "%s", data.c_str());
  return (data.size() > 0);
}
