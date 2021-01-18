#include <R_ext/Print.h>
#include "simple_re.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "helper.h"


simple_re*
simple_re_new()
{
	simple_re* new_re = (simple_re *) malloc( 1 * sizeof(simple_re));
	return new_re;
}

simple_re*
simple_re_compile( const char* pattern, const char* enc )
{
	OnigRegexType* regexp = malloc(sizeof(OnigRegexType));
	UChar* pattern_start = (UChar*) pattern; 
	int len = strlen(pattern);
	UChar* pattern_end = pattern_start + len ;
	OnigEncoding onig_enc = simple_re_obtain_onig_encoding( enc );
	OnigSyntaxType* onig_syntax = (OnigSyntaxType *) ONIG_SYNTAX_RUBY ;
	OnigErrorInfo einfo;

	int return_value;
	return_value = onig_new_without_alloc( regexp, pattern_start , pattern_end , ONIG_OPTION_SINGLELINE , onig_enc, onig_syntax , &einfo);

	if(return_value != ONIG_NORMAL){
		Rprintf("ERROR: Invalied regular epxression: %s \n", pattern);
	}

	simple_re* new_re = simple_re_new();
	new_re->regexp = regexp;
	char* copy_pattern = (char*) malloc(strlen(pattern) * sizeof(char) + 1 );
	memcpy( copy_pattern, pattern, strlen(pattern));
	copy_pattern[strlen(pattern)] = '\0';
	new_re->pattern = copy_pattern;
	new_re->encoding = enc;
	new_re->str = NULL;
	new_re->matched = NULL;
	return new_re;
}

int
simple_re_set_str( simple_re* re, const char* str)
{
	int len = strlen(str);
	char* new_str = (char*) malloc( (len + 1) * sizeof(char));
	strcpy(new_str, str);
	new_str[len] = '\0';
	re->str = (UChar*) new_str;
	return 1;
}

int
simple_re_set_matched_region( simple_re* re, OnigRegion* region)
{
	re->matched = region;
	return 1;
}

int
simple_re_match ( simple_re* re , const char* str , simple_re** pptr_for_last_rexp)
{
	DEBUG_PRINT ("onig_search is going to be executed. REGEXP: %s, STRING: %s\n", simple_re_read_pattern(re), str );

	simple_re_reset(re);
	simple_re_set_str( re, str );
	simple_re_set_matched_region( re, onig_region_new() );

	OnigRegexType* regexp = re->regexp;
	UChar* text =  re->str;
	unsigned int len , idx; // Here, what this code does is same as strlen(text)
	idx = 0;
	while( text[idx] != '\0'){
		idx = idx + 1;
	}
    len = idx; 
	const UChar* end_ptr = re->str + len ;
	const UChar* start_ptr = re->str ;
	OnigRegion* region = re->matched ; 

	OnigPosition return_value;
	return_value = onig_search(regexp, text, end_ptr, start_ptr, end_ptr, region , ONIG_OPTION_NONE);
	// printf("onig_search is executed: %d \n", (int)return_value);

	// This is important for back reference.
	*pptr_for_last_rexp = re ;

	if( return_value == ONIG_MISMATCH){
		DEBUG_PRINT("Regular expression is unmatched\n");
		return -1;
	}else {
		DEBUG_PRINT("Regular expression is matched\n");
		return 1;
	}
}

int
simple_re_matched_group_num( simple_re* re)
{
	if(re->matched != NULL){
		return re->matched->num_regs;
	} else {
		return -1 ;
	}
}

char*
simple_re_matched_str( simple_re* re , int idx )
{
	OnigRegion* matched;
	int num_groups;
	int length;
	char* matched_str;
	char* new_str;

	if(re->matched == NULL){
		Rprintf("WARNING: No available matched information. \n ");
		matched_str = NULL ;
		return matched_str ;
	}

	matched = re->matched;
	num_groups = simple_re_matched_group_num( re );
	if((idx < 0) | (idx > num_groups)){
		Rprintf("ERROR: Index is not within matched groups. \n");
		matched_str = NULL;
		return matched_str ;
	}else{
		length = matched->end[idx] - matched->beg[idx]  ; 
//		printf("Length: %d \n", length);
//		printf("from:%ld to:%ld length:%d \n", matched->beg[idx], matched->end[idx], length );
		new_str = (char*) malloc( (length + 1) * sizeof(char)) ;
		memcpy( new_str , re->str + matched->beg[idx] , length );
		new_str[length] = '\0';
		matched_str = new_str;
		// printf("Matched String: %s \n", new_str);
		return matched_str;
	}
}

int
simple_re_reset( simple_re* re )
{
	/* Reset only matching information */
	if(re->str != NULL){
		free(re->str );
		re->str = NULL;
	}
	if(re->matched != NULL){
		onig_region_free( re->matched , 1);
		re->matched = NULL;
	}
	return 1;
}


int
simple_re_free( simple_re* re )
{
	if(re->regexp != NULL){
		onig_free( re->regexp );
		re->regexp = NULL;
	}
	if(re->pattern != NULL){
		free(re->pattern );
		re->pattern = NULL;
	}
	if(re->str != NULL){
		free(re->str );
		re->str = NULL;
	}
	if(re->matched != NULL){
		onig_region_free( re->matched , 1);
		re->matched = NULL;
	}
	free(re);
	return 1;
}

const char*
simple_re_read_pattern( simple_re* re )
{
	return re->pattern;
}


/* --------------------------------------------------------------- */
// Not strictly part of the API, but useful for case-insensitive string comparison
extern int onigenc_with_ascii_strnicmp (OnigEncoding enc, const UChar *p, const UChar *end, const UChar *sascii, int n);

int ore_strnicmp (const char *str1, const char *str2, size_t num)
{
    return onigenc_with_ascii_strnicmp(ONIG_ENCODING_ASCII, (const UChar *) str1, (const UChar *) str1 + num, (const UChar *) str2, num);
}

// This function code is quoted from "ore" R package by jonclayden on Github.
// The original function signature is "OnigEncoding ore_name_to_onig_enc (const char *enc)"
OnigEncoding ore_name_to_onig_enc (const char *enc)
{
    if (ore_strnicmp(enc,"ASCII",5) == 0 || ore_strnicmp(enc,"US-ASCII",8) == 0)
        return ONIG_ENCODING_ASCII;
    else if (ore_strnicmp(enc,"UTF-8",5) == 0 || ore_strnicmp(enc,"UTF8",4) == 0)
        return ONIG_ENCODING_UTF8;
    else if (ore_strnicmp(enc,"ISO_8859-1",10) == 0 || ore_strnicmp(enc,"ISO-8859-1",10) == 0 || ore_strnicmp(enc,"ISO8859-1",9) == 0 || ore_strnicmp(enc,"LATIN1",6) == 0)
        return ONIG_ENCODING_ISO_8859_1;
    else if (ore_strnicmp(enc,"ISO_8859-2",10) == 0 || ore_strnicmp(enc,"ISO-8859-2",10) == 0 || ore_strnicmp(enc,"ISO8859-2",9) == 0 || ore_strnicmp(enc,"LATIN2",6) == 0)
        return ONIG_ENCODING_ISO_8859_2;
    else if (ore_strnicmp(enc,"ISO_8859-3",10) == 0 || ore_strnicmp(enc,"ISO-8859-3",10) == 0 || ore_strnicmp(enc,"ISO8859-3",9) == 0 || ore_strnicmp(enc,"LATIN3",6) == 0)
        return ONIG_ENCODING_ISO_8859_3;
    else if (ore_strnicmp(enc,"ISO_8859-4",10) == 0 || ore_strnicmp(enc,"ISO-8859-4",10) == 0 || ore_strnicmp(enc,"ISO8859-4",9) == 0 || ore_strnicmp(enc,"LATIN4",6) == 0)
        return ONIG_ENCODING_ISO_8859_4;
    else if (ore_strnicmp(enc,"ISO_8859-5",10) == 0 || ore_strnicmp(enc,"ISO-8859-5",10) == 0 || ore_strnicmp(enc,"ISO8859-5",9) == 0 || ore_strnicmp(enc,"LATIN5",6) == 0)
        return ONIG_ENCODING_ISO_8859_5;
    else if (ore_strnicmp(enc,"ISO_8859-6",10) == 0 || ore_strnicmp(enc,"ISO-8859-6",10) == 0 || ore_strnicmp(enc,"ISO8859-6",9) == 0 || ore_strnicmp(enc,"LATIN6",6) == 0)
        return ONIG_ENCODING_ISO_8859_6;
    else if (ore_strnicmp(enc,"ISO_8859-7",10) == 0 || ore_strnicmp(enc,"ISO-8859-7",10) == 0 || ore_strnicmp(enc,"ISO8859-7",9) == 0 || ore_strnicmp(enc,"LATIN7",6) == 0)
        return ONIG_ENCODING_ISO_8859_7;
    else if (ore_strnicmp(enc,"ISO_8859-8",10) == 0 || ore_strnicmp(enc,"ISO-8859-8",10) == 0 || ore_strnicmp(enc,"ISO8859-8",9) == 0 || ore_strnicmp(enc,"LATIN8",6) == 0)
        return ONIG_ENCODING_ISO_8859_8;
    else if (ore_strnicmp(enc,"ISO_8859-9",10) == 0 || ore_strnicmp(enc,"ISO-8859-9",10) == 0 || ore_strnicmp(enc,"ISO8859-9",9) == 0 || ore_strnicmp(enc,"LATIN9",6) == 0)
        return ONIG_ENCODING_ISO_8859_9;
    else if (ore_strnicmp(enc,"ISO_8859-10",11) == 0 || ore_strnicmp(enc,"ISO-8859-10",11) == 0 || ore_strnicmp(enc,"ISO8859-10",10) == 0 || ore_strnicmp(enc,"LATIN10",7) == 0)
        return ONIG_ENCODING_ISO_8859_10;
    else if (ore_strnicmp(enc,"ISO_8859-11",11) == 0 || ore_strnicmp(enc,"ISO-8859-11",11) == 0 || ore_strnicmp(enc,"ISO8859-11",10) == 0 || ore_strnicmp(enc,"LATIN11",7) == 0)
        return ONIG_ENCODING_ISO_8859_11;
    else if (ore_strnicmp(enc,"ISO_8859-13",11) == 0 || ore_strnicmp(enc,"ISO-8859-13",11) == 0 || ore_strnicmp(enc,"ISO8859-13",10) == 0 || ore_strnicmp(enc,"LATIN13",7) == 0)
        return ONIG_ENCODING_ISO_8859_13;
    else if (ore_strnicmp(enc,"ISO_8859-14",11) == 0 || ore_strnicmp(enc,"ISO-8859-14",11) == 0 || ore_strnicmp(enc,"ISO8859-14",10) == 0 || ore_strnicmp(enc,"LATIN14",7) == 0)
        return ONIG_ENCODING_ISO_8859_14;
    else if (ore_strnicmp(enc,"ISO_8859-15",11) == 0 || ore_strnicmp(enc,"ISO-8859-15",11) == 0 || ore_strnicmp(enc,"ISO8859-15",10) == 0 || ore_strnicmp(enc,"LATIN15",7) == 0)
        return ONIG_ENCODING_ISO_8859_15;
    else if (ore_strnicmp(enc,"ISO_8859-16",11) == 0 || ore_strnicmp(enc,"ISO-8859-16",11) == 0 || ore_strnicmp(enc,"ISO8859-16",10) == 0 || ore_strnicmp(enc,"LATIN16",7) == 0)
        return ONIG_ENCODING_ISO_8859_16;
    else if (ore_strnicmp(enc,"UTF-16BE",8) == 0)
        return ONIG_ENCODING_UTF16_BE;
    else if (ore_strnicmp(enc,"UTF-16LE",8) == 0)
        return ONIG_ENCODING_UTF16_LE;
    else if (ore_strnicmp(enc,"UTF-32BE",8) == 0)
        return ONIG_ENCODING_UTF32_BE;
    else if (ore_strnicmp(enc,"UTF-32LE",8) == 0)
        return ONIG_ENCODING_UTF32_LE;
    else if (ore_strnicmp(enc,"BIG5",4) == 0 || ore_strnicmp(enc,"BIG-5",5) == 0 || ore_strnicmp(enc,"BIGFIVE",7) == 0 || ore_strnicmp(enc,"BIG-FIVE",8) == 0)
        return ONIG_ENCODING_BIG5;
    else if (ore_strnicmp(enc,"CP932",5) == 0)
        return ONIG_ENCODING_CP932;
    else if (ore_strnicmp(enc,"CP1250",6) == 0 || ore_strnicmp(enc,"WINDOWS-1250",12) == 0)
        return ONIG_ENCODING_WINDOWS_1250;
    else if (ore_strnicmp(enc,"CP1251",6) == 0 || ore_strnicmp(enc,"WINDOWS-1251",12) == 0)
        return ONIG_ENCODING_WINDOWS_1251;
    else if (ore_strnicmp(enc,"CP1252",6) == 0 || ore_strnicmp(enc,"WINDOWS-1252",12) == 0)
        return ONIG_ENCODING_WINDOWS_1252;
    else if (ore_strnicmp(enc,"CP1253",6) == 0 || ore_strnicmp(enc,"WINDOWS-1253",12) == 0)
        return ONIG_ENCODING_WINDOWS_1253;
    else if (ore_strnicmp(enc,"CP1254",6) == 0 || ore_strnicmp(enc,"WINDOWS-1254",12) == 0)
        return ONIG_ENCODING_WINDOWS_1254;
    else if (ore_strnicmp(enc,"CP1257",6) == 0 || ore_strnicmp(enc,"WINDOWS-1257",12) == 0)
        return ONIG_ENCODING_WINDOWS_1257;
    else if (ore_strnicmp(enc,"EUC-JP",6) == 0 || ore_strnicmp(enc,"EUCJP",5) == 0)
        return ONIG_ENCODING_EUC_JP;
    else if (ore_strnicmp(enc,"EUC-KR",6) == 0 || ore_strnicmp(enc,"EUCKR",5) == 0)
        return ONIG_ENCODING_EUC_KR;
    else if (ore_strnicmp(enc,"EUC-TW",6) == 0 || ore_strnicmp(enc,"EUCTW",5) == 0)
        return ONIG_ENCODING_EUC_TW;
    else if (ore_strnicmp(enc,"GB18030",7) == 0)
        return ONIG_ENCODING_GB18030;
    else if (ore_strnicmp(enc,"KOI8-R",6) == 0)
        return ONIG_ENCODING_KOI8_R;
    else if (ore_strnicmp(enc,"KOI8-U",4) == 0)
        return ONIG_ENCODING_KOI8_U;
    else if (ore_strnicmp(enc,"SHIFT_JIS",9) == 0 || ore_strnicmp(enc,"SHIFT-JIS",9) == 0 || ore_strnicmp(enc,"SJIS",4) == 0)
        return ONIG_ENCODING_SJIS;
    else
    {
//        warning("Encoding \"%s\" is not supported by Oniguruma - using ASCII", enc);
        Rprintf("Encoding \"%s\" is not supported by Oniguruma - using ASCII", enc);        
		return ONIG_ENCODING_ASCII;
    }
}

// Wrapper function for ore_name_to_onig_enc(const char *enc)

OnigEncoding
simple_re_obtain_onig_encoding( const char* enc)
{
	return ore_name_to_onig_enc(enc);
}
