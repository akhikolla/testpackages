/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
vector<string>
medusa::splitstr(const string& s, const char delim) {
  vector<string> fields;
  
  /* Replace carriage returns and delimiters. */
  mdsize nbytes = 0;
  char* data = new char[s.size()];
  for(mdsize i = 0; i < s.size(); i++) {
    char c = s[i];
    if(c == '\r') continue;
    if(c == delim) c = '\0';
    data[nbytes] = c;
    nbytes++;
  }

  /* Ensure correct termination. */
  if(data[nbytes-1] == '\n') nbytes--;
  if(nbytes < 1) {
    delete [] data;
    return fields;
  }
  data[nbytes] = '\0';

  /* Collect string segments. */
  char* ptr = data;
  for(mdsize i = 0; i <= nbytes; i++) {
    if(data[i] != '\0') continue;
    fields.push_back(string(ptr));
    ptr = (char*)(data + i + 1);
  }
  delete [] data;
  return fields;
}
