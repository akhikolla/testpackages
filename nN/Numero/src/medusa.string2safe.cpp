/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 * Replaces most non-alphanumeric characters.
 * Truncates the string to the given maximum length. 
 */
string
medusa::string2safe(const std::string& orig, const mdsize maxlen) {

  /* Adjust length. */
  string s = orig;
  if(s.size() > maxlen) {
    s = s.substr(0, maxlen);
    if(maxlen > 1) {
      s[maxlen-1] = '.';
      s[maxlen-2] = '.';
    }
  }

  /* Replace characters. */
  for(mdsize i = 0; i < s.size(); i++) {
    if(isalnum(s[i])) continue;
    switch(s[i]) {
    case '-': continue;
    case '+': continue;
    case '=': continue;
    case '%': continue;
    case '.': continue;
    case ':': continue;
    case ';': continue;
    case ',': continue;
    case '/': continue;
    case '\\': continue;
    case '@': continue;
    case '(': continue;
    case ')': continue;
    case ' ': continue;
    default: s[i] = '_';
    }
  }  
  return s;
}
