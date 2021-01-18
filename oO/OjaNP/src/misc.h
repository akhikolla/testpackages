/* $Id: misc.h,v 1.1 2008/01/25 11:47:50 ruthe Exp $ */

#ifndef MISC_H
#define MISC_H

//#define TO_FILE
//#define DEEPDEBUG

#include <list>
#include <vector>
#include <iostream>

using namespace std;

#define check_bit(var,pos) ((var) & (1<<(pos)))

int fact(int k);
unsigned long choices(int n, int k);
double round(double x,double prec);
bool is_file(const char* s);

void get_partitions(list<vector<int> >& P,int n);
#endif
