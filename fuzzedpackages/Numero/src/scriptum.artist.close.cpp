/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
unsigned long
Artist::close(const string& script) {
  ArtistBuffer* p = (ArtistBuffer*)buffer;
  if(p->output == NULL) return p->filesize;

  /* Close all open groups. */
  while(this->group() > 0) {};

  /* Write custom script segment. */
  if(script.size() > 0)
    p->filesize += fwrite(script.c_str(), sizeof(char),
			  script.size(), p->output);

  /* Close an SVG file. */
  string endsvg = "\n</g>\n</svg>\n";
  p->filesize += fwrite(endsvg.c_str(), sizeof(char),
			endsvg.size(), p->output);
  
  /* Update prolog. */
  string protext = p->prolog();
  if(protext.size() != p->prosize)
    medusa::panic("Prolog failed.", __FILE__, __LINE__);

  /* Overwrite prolog. */
  rewind(p->output);
  fwrite(protext.c_str(), sizeof(char), protext.size(), p->output);

  /* Close file. */
  closefile(p->output);
  p->output = NULL;
  return p->filesize;
}
