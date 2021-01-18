/* $Id: global.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */
#include "misc.h"
#include <R.h>

#ifndef GLOBAL_H
#define GLOBAL_H

//using namespace std; //df

// bit flags
// 0 - no debug
// 1 - all
// 2 - time
extern int debug;
extern bool trace;
extern bool verbose;
extern bool quiet;
extern bool adaptive;

#ifndef _MSC_VER
#include <Rcpp.h> 
#define cout Rcpp::Rcout
#endif


#ifdef TO_FILE
#define LOG(a) {}
#define FLOG(a) {if(debug&1) fout << a << endl;}
#else
#define LOG(a) {if(debug&1) cout << a << endl;}
#define FLOG(a) {}
#endif

#define LOGTIME(a) {if (debug&2) cout << a << endl;}

#define TRACE(a) {}//if(trace)cerr << a << endl;}
#define LOGIF(a,b) {if(debug && (a)) cout << b << endl;}

//#define LOG(a) {}
//#define TRACE(a) {}
//#define LOGIF(a,b) {}

bool parse_arguments(int argc,char** argv);

const char* program_name();

int options();
const char* option(int idx);

int parameters();
const char* parameter(int idx);
	 
bool has_option(const char* opt);
int int_option(const char* opt,int def=0);
double double_option(const char* opt,double def=0);
const char* string_option(const char* opt,const char* def="");

#endif
