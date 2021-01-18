/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
mdsize
Artist::group() {
  ArtistBuffer* p = (ArtistBuffer*)buffer;
  if(p->ngroups > 0) { 
    p->filesize += fprintf(p->output, "\n</g>\n");
    p->ngroups -= 1;
  }
  return p->ngroups;
}

/*
 *
 */
mdsize
Artist::group(const string& key) {
  ArtistBuffer* p = (ArtistBuffer*)buffer;
  if(key.size() < 1) p->filesize += fprintf(p->output, "\n<g>\n");
  else p->filesize += fprintf(p->output, "\n<g id=\"%s\">\n", key.c_str());
  p->ngroups += 1;
  return p->ngroups;
}
