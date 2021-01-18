/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

static string convert_date(const string&);

/*
 *
 */
mdreal
medusa::string2real(const string& x) {
  mdsize i; 
  mdreal rl_nan = medusa::rnan();

  /* Check if date (YYYY-MM-DD). */
  string s = convert_date(x);
  mdsize n = s.size();
  if(n < 1) return rl_nan;

  /* Replace commas. */
  for(i = 0; i < n; i++)
    if(s[i] == ',') s[i] = '.';
  
  /* Remove preceding space. */
  const char* buf = s.c_str(); 
  for(i = 0; i < n; i++)
    if(isspace(buf[i]) == false) break;
  switch(buf[i]) {
  case '-': {i++; break;}
  case '+': {i++; break;}
  case '.': {i++; break;}
  }

  /* Check if valid representation. */
  if(isdigit(buf[i]) == false) return rl_nan;

  /* Check if within range. */
  double value = atof(buf);
  if(fabs(value) >= rl_nan) return rl_nan;
  return value;
}

/*
 *
 */
string
convert_date(const string& x) {
  if(x.size() != 10) return x;
  if(x[4] != '-') return x;
  if(x[7] != '-') return x;
  double year = atof(x.c_str());
  double month = atof((x.substr(5, 2)).c_str());
  double day = atof((x.substr(8, 2)).c_str());
  if((month < 1) || (month > 12)) return "";
  if((day < 1) || (day > 31)) return "";
  char buffer[32];
  double stamp = (year + (month - 1.0)/12.0 + (day - 1.0)/365.25);
  sprintf(buffer, "%.10e", stamp);
  return string(buffer);
}
