/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef medusa_local_INCLUDED
#define medusa_local_INCLUDED

#define _FILE_OFFSET_BITS 64 /* large file support */
#define IOBUFCAP_medusa 262144

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <climits>
#include <ctime>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "medusa.h"

using namespace std;
using namespace medusa;

/* Encapsulate with redundant namespace in case in a collection
   of modules another module has the same class name(s) in use. */
namespace medusa_local {
  
  /*
   *
   */
  class FileBuffer {
  public:
    FILE* fid;
    string name;
    string errtext;
    unsigned long nread;
    unsigned long nwritten;
    char iobuf[IOBUFCAP_medusa];
    char bytes[IOBUFCAP_medusa];
    static unordered_map<FILE*, char*> handles;
  public:
    FileBuffer() {
      this->fid = NULL;
      this->nread = 0;
      this->nwritten = 0;
    };
    ~FileBuffer() {this->clear();};
    void clear() {
      if(fid != NULL) {fclose(fid); this->fid = NULL;}
      this->nread = 0;
      this->nwritten = 0;
      (this->errtext).clear();
    };
    void abort(const string& s) {
      if(fid != NULL) {fclose(fid); this->fid = NULL;}
      if(errtext.size() < 1) this->errtext = s;
    };
  };    

  /*
   *
   */
  class TableBuffer {
  private:
    mdsize counter;
    unordered_map<string, pair<mdsize, mdsize> > word2rank;
  public:
    unordered_map<mdsize, string> words;
    unordered_map<mdsize, unordered_map<mdsize, mdsize> > data;
  public:
    TableBuffer() {};
    TableBuffer(void* ptr) {
      TableBuffer* p = (TableBuffer*)ptr;
      this->counter = p->counter;
      this->word2rank = p->word2rank;
      this->words = p->words;
      this->data = p->data;
    };
    ~TableBuffer() {};
    mdsize grow(const string&);
    mdsize shrink(const mdsize);
  };

  /* Global list of file pointers. */
}

using namespace medusa_local;

#endif /* medusa_local_INCLUDED */
