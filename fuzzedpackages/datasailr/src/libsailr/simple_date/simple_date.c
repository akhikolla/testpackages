#include <R_ext/Print.h>
#include "simple_date.h"
#include "cpp_date.hpp"
#include <string.h>
#include <stdio.h>
#include <ctype.h>

int
simple_date_ymd( int y , int m, int d )
{
  return cpp_date_ymd(y, m, d);
}

int
simple_date_ym_weekday_nth( int int_y, unsigned int int_m, const char* c_wd , unsigned int int_nth )
{
  unsigned int int_wd;

  // Upper case
  int idx = 0;
  char c;
  char str[10];
  while(c_wd[idx]) {
    // printf("%c\n", c_wd[idx]);
    c = ( (char) toupper(c_wd[idx]) );
    // printf("%c\n", c);
    str[idx] = c;
    // printf("%c\n", str[idx]);
    idx++;
  }
  str[idx] = '\0';

  // Assign integer to int_wd;
  if( strcmp( str, "SUN" ) == 0 ){
    int_wd = 0;
  }else if( strcmp( str, "MON" ) == 0 ){
    int_wd = 1;
  }else if( strcmp( str, "TUE" ) == 0 ){
    int_wd = 2;
  }else if( strcmp( str, "WED" ) == 0 ){
    int_wd = 3;
  }else if( strcmp( str, "THU" ) == 0 ){
    int_wd = 4;
  }else if( strcmp( str, "FRI" ) == 0 ){
    int_wd = 5;
  }else if( strcmp( str, "SAT" ) == 0 ){
    int_wd = 6;
  }else{
    int_wd = 0; // This branch should never be executed.
    Rprintf("ERROR: Specified symbol is not valid for weekday. %s\n", str );
  }
  return cpp_date_ym_weekday_nth( int_y, int_m, int_wd, int_nth); 
}

int
simple_date_add_n_years( int unix_date , int years)
{
  return cpp_date_add_n_years( unix_date, years);
}

int
simple_date_add_n_months( int unix_date , int months)
{
  return cpp_date_add_n_months( unix_date , months);
}

int
simple_date_add_n_days( int unix_date , int days )
{
  return cpp_date_add_n_days( unix_date , days );
}

char*
simple_date_new_cstr_format ( int unix_date, const char* fmt  )
{
  return cpp_date_new_cstr_format ( unix_date, fmt );
}

