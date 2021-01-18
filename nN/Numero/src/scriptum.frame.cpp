/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "scriptum.local.h"

/*
 *
 */
Frame::Frame() {
  this->buffer = new FrameBuffer();
}

/*
 *
 */
Frame::Frame(const Frame& t) {
  this->buffer = new FrameBuffer(t.buffer);
}

/*
 *
 */
void
Frame::operator=(const Frame& t) {
  if(this == &t) return;
  FrameBuffer* p = (FrameBuffer*)buffer; delete p;
  this->buffer = new FrameBuffer(t.buffer);
}

/*
 *
 */
Frame::~Frame() {
  FrameBuffer* p = (FrameBuffer*)buffer;
  delete p;
}

