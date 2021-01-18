/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
void
Artist::background(const Color& c) {
  ArtistBuffer* p = (ArtistBuffer*)buffer;
  p->bgcolor = c;
}
