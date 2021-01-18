#ifndef SIMPLE_DATE_H
#define SIMPLE_DATE_H

int simple_date_ymd( int y , int m, int d ); 
int simple_date_ym_weekday_nth( int int_y, unsigned int int_m, const char* c_wd , unsigned int int_nth );
int simple_date_add_n_years( int unix_date , int years);
int simple_date_add_n_months( int unix_date , int months);
int simple_date_add_n_days( int unix_date , int days );
char* simple_date_new_cstr_format ( int unix_date, const char* fmt  );

#endif // SIMPLE_DATE_H
