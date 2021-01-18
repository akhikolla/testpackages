/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
Engine::Engine() {
  this->buffer = new EngineBuffer();
}

/*
 *
 */
Engine::Engine(const Topology& tp) {
  EngineBuffer* p = new EngineBuffer();
  p->topology = tp;
  this->buffer = p;
}

/*
 *
 */
Engine::Engine(const Engine& t) {
  this->buffer = new EngineBuffer(t.buffer);
}

/*
 *
 */
void
Engine::operator=(const Engine& t) {
  if(this == &t) return;
  EngineBuffer* p = (EngineBuffer*)buffer; delete p;
  this->buffer = new EngineBuffer(t.buffer);
}

/*
 *
 */
Engine::~Engine() {
  EngineBuffer* p = (EngineBuffer*)buffer;
  delete p;
}
