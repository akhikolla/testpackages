/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
TwowayMap::erase(const mdsize r) {
  if(rank2name.count(r) < 1) return;
  string s = rank2name[r];
  (this->rank2name).erase(r);
  (this->name2rank).erase(s);
}

/*
 *
 */
void
TwowayMap::erase(const string& s) {
  if(name2rank.count(s) < 1) return;
  mdsize r = name2rank[s];
  (this->rank2name).erase(r);
  (this->name2rank).erase(s);
}

/*
 *
 */
void
TwowayMap::insert(const mdsize rnew, const string& snew) {
  if(rank2name.count(rnew) > 0) {
    string s = rank2name[rnew];
    mdsize r = name2rank[s];
    (this->rank2name).erase(r);
    (this->name2rank).erase(s);
  }
  if(name2rank.count(snew) > 0) {
    mdsize r = name2rank[snew];
    string s = rank2name[r];
    (this->rank2name).erase(r);
    (this->name2rank).erase(s);
  }
  this->rank2name[rnew] = snew;
  this->name2rank[snew] = rnew;
}

/*
 *
 */
string
TwowayMap::name(const mdsize r) {
  if(rank2name.count(r) < 1) return "";
  return rank2name[r];
}

/*
 *
 */
mdsize
TwowayMap::rank(const string& s) {
  if(name2rank.count(s) < 1) return medusa::snan();
  return name2rank[s];
}
