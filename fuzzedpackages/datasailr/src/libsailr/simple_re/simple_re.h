#ifndef SIMPLE_RE_H
#define SIMPLE_RE_H

#include "onigmo.h"

#define OUT 

// This regular expression program is intended to match once, not multiple times.
// That's because this is called "simple" re (regular expression).

typedef struct simple_re_ {
	OnigRegexType* regexp;
	char* pattern;
	const char* encoding;
	UChar* str;
	OnigRegion* matched;
} simple_re ;

// Function Prototypes
simple_re* simple_re_compile( const char* pattern, const char* enc );
int simple_re_match ( simple_re* re, const char* str, simple_re** ptr_for_last_rexp);
int simple_re_reset( simple_re* re);
int simple_re_free( simple_re* re);
int simple_re_matched_group_num( simple_re* re);
char* simple_re_matched_str( simple_re* re , int idx );

const char* simple_re_read_pattern( simple_re* re);

// Utility function
OnigEncoding simple_re_obtain_onig_encoding( const char* );

#endif  /* SIMPLE_RE_H */

