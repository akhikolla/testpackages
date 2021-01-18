#include <R_ext/Print.h>
#include <stdio.h>
#include "cpp_string.hpp"


int main(int argc, char** argv)
{
	cpp_object* obj = cpp_string_new("Hello ");
	cpp_object* obj2 = cpp_string_new("World!");
	cpp_object* obj3 = cpp_string_concat(obj, obj2);
	Rprintf( "%s\n" , cpp_string_read(obj3));
	Rprintf( "%p\n" , obj);
	Rprintf( "%p\n" , obj2);
	Rprintf( "%p\n" , obj3);
	cpp_string_free(obj);
	cpp_string_free(obj2);

	cpp_object* obj4 = cpp_string_new("KONNICHIWA");
	cpp_string_move_ptr ( &obj3, &obj4);
	Rprintf( "%s\n" , cpp_string_read (obj3));
	Rprintf( "%s\n" , cpp_string_read (obj4)); // It's deleted

	cpp_object* obj6 = cpp_string_subset(obj3, 3, 5);
	Rprintf( "%s\n" , cpp_string_read (obj6));
	Rprintf( "%p\n" , obj6);

	cpp_object* obj7 = cpp_string_repeat(obj6, 3);
	Rprintf( "%s\n" , cpp_string_read (obj7));
	Rprintf( "%p\n" , obj7);

	cpp_object* obj8 = cpp_string_int2str( 123 );
	cpp_object* obj9 = cpp_string_double2str( 123.456 );
	Rprintf( "%s\n" , cpp_string_read (obj8));
	Rprintf( "%s\n" , cpp_string_read (obj9));

	cpp_object* obj10 = cpp_string_new("  Hello World   \n");
	cpp_object* obj11 = cpp_string_strip(obj10);
	cpp_object* obj12 = cpp_string_lstrip(obj10);
	cpp_object* obj13 = cpp_string_rstrip(obj10);
	Rprintf("XXX");
	Rprintf( "%s" , cpp_string_read (obj10));
	Rprintf("XXX");
	Rprintf( "%s" , cpp_string_read (obj11));
	Rprintf("XXX");
	Rprintf( "%s" , cpp_string_read (obj12));
	Rprintf("XXX");
	Rprintf( "%s" , cpp_string_read (obj13));
	Rprintf("XXX");

	cpp_string_free(obj3); // obj5 is also deleted.
	cpp_string_free(obj6);
	cpp_string_free(obj7);
	cpp_string_free(obj8);
	cpp_string_free(obj9);
	cpp_string_free(obj10);
	cpp_string_free(obj11);
	cpp_string_free(obj12);
	cpp_string_free(obj13);
}
