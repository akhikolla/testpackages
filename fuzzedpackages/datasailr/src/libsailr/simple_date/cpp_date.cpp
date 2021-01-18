#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include "date.h"
#include "cpp_date.hpp"
#include <chrono>

// I am using date.h that is a great C++ header only library based on the standard library of chrono.
// It seems that this library will be 
// 
// Roughly speaking time has two aspects, 1. time point and 2. period (time interval) .
// 
// 1. For time points the following types can be used.
//
// * (field-based) date::year_month_day
// * (serial-based) date::sys_days
// * (serial-based) date::year, date::month, date::day  (The other field information is not held.)
// 
// date::year_month_day holds inforamtion year, month, day seperately.
// date::sys_days holds one information, how many days have passed from specific date (UNIX time).
// 
// About calculation,
// date::sys_days is used for day level calculation. 
// date::year_month_day is used for year or month level calculation.
// 
// date::sys_days is appropriate for date calculation, because it's serial based. 
// date::year, date::month and date::day can be used for calculation.
// From serial based types, you can obtain values using .count() method.
// 
// 
// 2. For time intervals 
//
// * date::days
// * date::years, date::months
// 
// Calculating intervals from date::sys_days, you can get date::days.
// You can also get date::days from date::day calculation.
// date::years from date::year calculation, date::months from date::month calculation.
// 
// You cannot get date::days from date::year_month_day directly. Before calculation you need conversions, such as 
// a. date::sys_days{ date_ymd_obj } => date::sys_days
// b. date_ymd_obj.year() => date::year
// c. date_ymd_obj.month() => date::month
// d. date_ymd_obj.day() => date::day
// 
//
// 3. Constructors of these types
// 
// 3.1 For date::year_month_day
// 
// date::year{2019}/1/1 
// date::year_month_day{ date::year{ 2019 }, date::month{ 1 } , date::day{ 1 } }
// 
// 3.2 For date::sys_days
// 
// date::sys_days{ date_ymd_obj };
// 
// 


// Private functions 
date::sys_days
obtain_unix_epoch_sys_days()
{
  return date::sys_days{ date::year{1970}/1/1 };
}

date::days
convert_ymd_to_unix_date(date::year_month_day ymd)
{
  date::sys_days specified = date::sys_days{ ymd } ;
  date::sys_days unix_base = obtain_unix_epoch_sys_days();
  return (specified - unix_base);
}

date::days
convert_ymdi_to_unix_date(date::year_month_weekday ymdi)
{
  date::sys_days specified = date::sys_days{ ymdi } ;
  date::sys_days unix_base = obtain_unix_epoch_sys_days();
  return (specified - unix_base);
}

date::days
convert_sys_days_to_unix_date( date::sys_days sd )
{
  date::sys_days unix_base = obtain_unix_epoch_sys_days();
  date::days result = sd - unix_base; 
  return result;
}


date::sys_days
convert_unix_date_to_sys_days(int unix_date)
{
  date::sys_days unix_base = obtain_unix_epoch_sys_days();
  date::sys_days sys_day = unix_base + date::days{ unix_date };
  return sys_day ;
}

// Public functions

char*
cpp_date_new_cstr_format ( int unix_date, const char* fmt  )
{
  date::sys_days base_day = obtain_unix_epoch_sys_days();
  date::sys_days new_day = base_day + date::days{unix_date};
  std::stringstream ss;
  ss << date::format( fmt, new_day ) ;
  const char* const_str = ss.str().c_str();
  char* new_str = (char*) malloc( (strlen(const_str) + 1) * sizeof(char) );
  strcpy(new_str, const_str);
  return new_str;
}

int
cpp_date_ymd( int y , int m, int d )
{
  date::year_month_day ymd_obj = date::year{y}/m/d;;
  date::days unix_date = convert_ymd_to_unix_date(ymd_obj);
  return unix_date.count() ;
}

int
cpp_date_ym_weekday_nth( int int_y, unsigned int int_m, unsigned int int_wd , unsigned int int_nth)
{
// For wd, normal range is [0, 6], for Sunday through Saturday.
  date::year y = date::year{int_y};
  date::month m = date::month{int_m};
  date::weekday wd = date::weekday{int_wd};
  date::weekday_indexed wdi = wd[int_nth]; /* index n in the range [1, 5]. It represents the first, second, third, fourth, or fifth weekday of some month.  */
  date::year_month_weekday ymdi = y/m/wdi;
  date::days unix_date = convert_ymdi_to_unix_date(ymdi);
  return unix_date.count() ;
}

int
cpp_date_add_n_years( int unix_date , int years)
{
  date::sys_days sd = convert_unix_date_to_sys_days(unix_date);
  date::year_month_day ymd = date::year_month_day{ sd };
  date::years y = date::years{years};
  date::year_month_day new_ymd = ymd + y;
  date::days new_days = convert_ymd_to_unix_date( new_ymd );
  return new_days.count() ;
}

int
cpp_date_add_n_months( int unix_date , int months)
{
  date::sys_days sd = convert_unix_date_to_sys_days(unix_date);
  date::year_month_day ymd = date::year_month_day{ sd };
  date::months m = date::months{months};
  date::year_month_day new_ymd = ymd + m;
  date::days new_days = convert_ymd_to_unix_date( new_ymd );
  return new_days.count() ;
}

int
cpp_date_add_n_days( int unix_date , int days )
{
  date::sys_days sd = convert_unix_date_to_sys_days(unix_date);
  date::sys_days new_sd = sd + date::days{days};
  date::days unix_day = convert_sys_days_to_unix_date( new_sd ) ;
  return ( unix_day.count() );

}





