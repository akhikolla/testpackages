/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
FrameBuffer::FrameBuffer() {
  this->ngroups = 0;
  this->bytes[0] = '\0';
}

/*
 *
 */
FrameBuffer::FrameBuffer(const void* ptr) {
  const FrameBuffer* p = (FrameBuffer*)ptr;
  this->ngroups = p->ngroups;
  this->data = (p->data + p->bytes);
  this->bytes[0] = '\0';
  this->limits = p->limits;
  this->style = p->style;
  this->linestycode = p->linestycode;
  this->textstycode = p->textstycode;
}

/*
 *
 */
FrameBuffer::~FrameBuffer() {}

/*
 *
 */
void
FrameBuffer::append(const string& s) {
  (this->data).append(bytes);
  (this->data).append(s);
  this->bytes[0] = '\0';
}

/*
 * Pointer for formatted printing.
 */
char*
FrameBuffer::f() {
  (this->data).append(bytes);
  this->bytes[0] = '\0';
  return bytes;
}

/*
 *
 */
string
FrameBuffer::flush() {
  string out = (data + bytes);
  this->ngroups = 0;
  (this->limits).first = Limes();
  (this->limits).second = Limes();
  this->style = Style();
  (this->linestycode).clear();
  (this->textstycode).clear();
  (this->data).clear();
  this->bytes[0] = '\0';
  return out;
}
