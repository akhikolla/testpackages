/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef scriptum_local_INCLUDED
#define scriptum_local_INCLUDED

#include <cstdlib>
#include <cstdio>
#include <cfloat>
#include <cctype>
#include <ctime>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include "medusa.h"
#include "abacus.h"
#include "scriptum.h"

#define SAFETY_scriptum 65535
#define MINCOORD_scriptum -49999.0
#define MAXCOORD_scriptum 49999.0

using namespace std;
using namespace medusa;
using namespace abacus;
using namespace scriptum;

/* Encapsulate with redundant namespace in case in a collection
   of modules another module has the same class name(s) in use. */
namespace scriptum_local {
  
  /*
   *
   */
  class Limes {
  public:
    mdreal alpha;
    mdreal omega;
  public:
    Limes();
    ~Limes();
    bool update(const mdreal);
    bool update(const mdreal, const Style&);
    bool update(const vector<mdreal>&, const Style&);
  };

  /*
   *
   */
  class FrameBuffer {
  private:
    char bytes[SAFETY_scriptum];
    string data;
  public:
    mdsize ngroups;
    Style style;
    string linestycode;
    string textstycode;
    pair<Limes, Limes> limits;
  public:
    FrameBuffer();
    FrameBuffer(const void*);
    ~FrameBuffer();
    void append(const string&);
    char* f(); /* pointer for formatted printing */
    string flush();
  };

  /*
   *
   */
  class ArtistBuffer {
  public:
    mdsize ngroups;
    Color bgcolor;
    unsigned long counter;
    unsigned long filesize;
    unsigned long prosize;
    pair<Limes, Limes> limits;
    FILE* output;
  public:
    ArtistBuffer() {
      this->ngroups = 0;
      this->counter = 0;
      this->filesize = 0;
      this->prosize = 0;
      this->bgcolor = Color("#ffffff");
      this->output = NULL;
    };
    ArtistBuffer(const void* ptr) {
      ArtistBuffer* p = (ArtistBuffer*)ptr;
      if(p->output != NULL)
	panic("Cannot copy active object.\n", __FILE__, __LINE__);
      this->ngroups = p->ngroups;
      this->counter = p->counter;
      this->filesize = p->filesize;
      this->prosize = p->prosize;
      this->limits = p->limits;
      this->output = p->output;
    };
    ~ArtistBuffer() {
      if(output != NULL) closefile(output);
    };
    string prolog() const;
  };

  /* Utility functions. */
  extern void style2code(string&, string&, const Style&);
}

using namespace scriptum_local;

#endif /* scriptum_local_INCLUDED */
