#ifndef DATE_HPP
#define DATE_HPP

#ifdef __cplusplus
extern "C" {
#endif

int cpp_date_ymd( int y , int m, int d ); 
int cpp_date_ym_weekday_nth( int int_y, unsigned int int_m, unsigned int int_nth , unsigned int int_wd );
int cpp_date_add_n_years( int unix_date , int years);
int cpp_date_add_n_months( int unix_date , int months);
int cpp_date_add_n_days( int unix_date , int days );
char* cpp_date_new_cstr_format( int unix_date, const char* fmt  );

#ifdef __cplusplus
}
#endif

#endif // DATE_HPP






