#ifndef STRUCT_STRING_H
#define STRUCT_STRING_H

/* Struct & Functions for String */

typedef struct struct_string {
  int len;
  char buf[0];
} struct_string, StructString, *pStructString ;

pStructString new_struct_string( char* , int  );
StructString* cat_struct_2strings( struct_string*, struct_string*);
struct_string* cat_struct_strings( int, ... );
void free_struct_string(struct_string*);

#endif /* STRUCT_STRING_H */
