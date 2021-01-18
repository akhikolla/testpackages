/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
Artist::Artist() {
  this->buffer = new ArtistBuffer();
}

/*
 *
 */
Artist::Artist(const string& fname) {
  ArtistBuffer* p = new ArtistBuffer();
  this->buffer = p;

  /* Open output file. */
  p->output = openfile(fname, "w");
  if(p->output == NULL) {
    worry(("Cannot open '" + fname + "'."), "");
    return;
  }
  
  /* Print prolog. */
  string protext = p->prolog();
  p->filesize += fprintf(p->output, "%s", protext.c_str());
  p->prosize = protext.size();
}

/*
 *
 */
Artist::Artist(const Artist& t) {
  this->buffer = new ArtistBuffer(t.buffer);
}

/*
 *
 */
void
Artist::operator=(const Artist& t) {
  if(this == &t) return;
  ArtistBuffer* p = (ArtistBuffer*)buffer; delete p;
  this->buffer = new ArtistBuffer(t.buffer);
}

/*
 *
 */
Artist::~Artist() {
  ArtistBuffer* p = (ArtistBuffer*)buffer;
  delete p;
}

