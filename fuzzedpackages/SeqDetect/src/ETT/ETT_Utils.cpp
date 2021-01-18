#include <unordered_map>
#include <sstream>
#include <string>
#include <algorithm>
#include <Rmath.h>
#include <R.h>
#include <random>
#include "ETT_Utils.hpp"
using namespace std;

string generate_hex(const unsigned int len) {
  stringstream ss;
  GetRNGstate();
  for(unsigned i = 0; i < len; i++) {
    const unsigned int rc = (unsigned)(unif_rand() * 255);
    //const unsigned int rc = (unsigned)(rand() % 255);
    stringstream hexstream;
    hexstream << hex << rc;
    auto hex = hexstream.str();
    ss << (hex.length() < 2 ? '0' + hex : hex);
  }
  PutRNGstate();
  return ss.str();
}

set<string> *ett_set_intersect(set<string> *s1,set<string> *s2,bool remove1,bool remove2) {
  set<string> *res=new set<string>();
  if(s1!=NULL && s2!=NULL) {
      for(set<string>::const_iterator p = s1->begin( );p != s1->end( ); ++p)
          if(s2->find(*p)!=s2->end()) res->insert(*p);
      if(remove1) delete s1;
      if(remove2) delete s2;
  }
  return res;
}

set<string> *ett_set_union(set<string> *s1,set<string> *s2,bool remove1,bool remove2) {
  set<string> *res=new set<string>();
  if(s1!=NULL) {
    res->insert(s1->begin(),s1->end());
    if(remove1) delete s1;
  }
  if(s2!=NULL) {
    res->insert(s2->begin(),s2->end());
    if(remove2) delete s2;
  }
  return res;
}

set<string> *ett_set_diff(set<string> *s1,set<string> *s2,bool remove1,bool remove2) {
  set<string> *res=new set<string>();
  if(s1!=NULL && s2!=NULL) {
      for(set<string>::const_iterator p = s1->begin( );p != s1->end( ); ++p)
          if(s2->find(*p)==s2->end()) res->insert(*p);
      if(remove1) delete s1;
      if(remove2) delete s2;
  }
  return res;
}

set<string> *ett_set_clone(set<string> *s1) {
  set<string> *res=new set<string>();
  *res=*s1;
  return res;
}

string formatVector(vector<string> *ss,bool remove) {
    if(ss==NULL) return "[NULL]";
    stringstream strs;
    strs << "[";
    unsigned i=0;
    for(string s:*ss) {
        if(i<(ss->size()-1)) strs << s << ",";
        else strs << s;
        i++;
    }
    strs << "]";
    if(remove) delete ss;
    return strs.str();
}

string formatSet(set<string> *ss,bool remove) {
    vector<string> *tv=new vector<string>(ss->begin(),ss->end());
    string res=formatVector(tv);
    delete tv;
    if(remove) delete ss;
    return res;
}
