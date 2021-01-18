/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
Point::Point() {
  this->key = medusa::snan();
  this->home = medusa::snan();
}

/*
 *
 */
Point::Point(const mdsize rank, const vector<mdreal>& array,
	     const mdsize unit) {
  this->key = rank;
  this->home = unit;
  this->contents = array;
}

/*
 *
 */
Point::~Point() {}

/*
 *
 */
vector<mdreal>
Point::data() const {
  return contents;
}

/*
 *
 */
mdsize
Point::location() const {
  return home;
}

/*
 *
 */
void
Point::move(const mdsize unit) {
  this->home = unit;
}

/*
 *
 */
mdsize
Point::rank() const {
  return key;
}
