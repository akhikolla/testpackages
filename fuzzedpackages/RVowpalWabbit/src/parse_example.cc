/*
Copyright (c) 2009 Yahoo! Inc.  All rights reserved.  The copyrights
embodied in the content of this file are licensed under the BSD
(revised) open source license
 */

#include <cstring>		// bzero()
#include <math.h>
#include <ctype.h>
#include "parse_example.h"
#include "hash.h"
#include "cache.h"
#include "unique_sort.h"

#include <Rcpp.h>
#define VWCOUT Rcpp::Rcout

using namespace std;

size_t hashstring (substring s, unsigned long h)
{
  size_t ret = 0;
  //trim leading whitespace
  for(; *(s.start) <= 0x20 && s.start < s.end; s.start++);
  //trim trailing white space
  for(; *(s.end-1) <= 0x20 && s.end > s.start; s.end--);

  char *p = s.start;
  while (p != s.end)
    if (isdigit(*p))
      ret = 10*ret + *(p++) - '0';
    else
      return uniform_hash((unsigned char *)s.start, s.end - s.start, h);

  return ret + h;
}

size_t hashall (substring s, unsigned long h)
{
  return uniform_hash((unsigned char *)s.start, s.end - s.start, h);
}

hash_func_t getHasher(const string& s){
  if (s=="strings")
    return hashstring;
  else if(s=="all")
    return hashall;
  else{
    Rf_error("Unknown hash function: %s. Exiting", s.c_str());
    return hashstring;// not reached, just for -Wall
  }
  return hashstring;// not reached, just for -Wall
}

void feature_value(substring &s, v_array<substring>& name, float &v)
{
  tokenize(':', s, name);
  
  switch (name.index()) {
  case 0:
  case 1:
    v = 1.;
    break;
  case 2:
    v = float_of_substring(name[1]);
    if ( ISNAN(v))
      {
	VWCOUT << "error NaN value for feature: ";
	//cout.write(name[0].start, name[0].end - name[0].start);
	VWCOUT << " terminating." << endl;
	Rf_error("");
      }
    break;
  default:
    VWCOUT << "example with a wierd name.  What is ";
    //cerr.write(s.start, s.end - s.start);
    VWCOUT << "\n";
  }
}

// new should be empty. --Alex
template<class T> void copy_v_array(v_array<T>& old, v_array<T>& new_va)
{
  // allocate new memory for new.
  new_va.begin = (T *) malloc(sizeof(T)*old.index());

  // copy the old to the new.
  memcpy(new_va.begin,old.begin,sizeof(T)*old.index());
  
  // update pointers
  new_va.end = new_va.begin + old.index();
  new_va.end_array = new_va.end;
}

char* c_string_of_substring(substring s)
{
  size_t len = s.end - s.start+1;
  char* ret = (char *)calloc(len,sizeof(char));
  memcpy(ret,s.start,len-1);
  return ret;
}

char* copy(char* base)
{
  size_t len = 0;
  while (base[len++] != '\0');
  char* ret = (char *)calloc(len,sizeof(char));
  memcpy(ret,base,len);
  return ret;
}

int read_features(parser* p, void* ex)
{
  example* ae = (example*)ex;
  char *line=NULL;
  int num_chars = readto(*(p->input), line, '\n');
  if (num_chars == 0)
    return num_chars;
  
  substring example = {line, line + num_chars};

  tokenize('|', example, p->channels);
  p->lp->default_label(ae->ld);
  substring* feature_start = &(p->channels[1]);
  if (*line == '|')
    feature_start = &(p->channels[0]);
  else 
    {
      substring label_space = p->channels[0];
      char* tab_location = safe_index(label_space.start, '\t', label_space.end);
      if (tab_location != label_space.end)
	label_space.start = tab_location+1;
      
      tokenize(' ',label_space,p->words);
      if (p->words.index() > 0 && p->words.last().end == label_space.end) //The last field is a tag, so record and strip it off
	{
	  substring tag = p->words.pop();
	  push_many(ae->tag, tag.start, tag.end - tag.start);
	}
      
      p->lp->parse_label(ae->ld, p->words);
    }
  
  size_t mask = global.mask;
  bool audit = global.audit;
  for (substring* i = feature_start; i != p->channels.end; i++) {
    substring channel = *i;
    
    tokenize(' ',channel, p->words);
    if (p->words.begin == p->words.end)
      continue;
    
    float channel_v = 1.;
    size_t channel_hash;
    char* base=NULL;
    size_t index = 0;
    bool new_index = false;
    size_t feature_offset = 0;
    if (channel.start[0] != ' ')//we have a nonanonymous namespace
      {
	feature_offset++;
	feature_value(p->words[0], p->name, channel_v);

	if (p->name.index() > 0) {
	  index = (unsigned char)(*p->name[0].start);
	  if (ae->atomics[index].begin == ae->atomics[index].end)
	    new_index = true;
	}
	if (audit)
	  base = c_string_of_substring(p->name[0]);
	channel_hash = p->hasher(p->name[0], hash_base);
      }
    else
      {
	index = (unsigned char)' ';
	if (ae->atomics[index].begin == ae->atomics[index].end)
	  new_index = true;
	if (audit)
	  {
	    base = (char *)calloc(2,sizeof(char));
	    base[0]=' ';
	    base[1]='\0';
	  }
	channel_hash = 0;
      }
 
    for (substring* i = p->words.begin+feature_offset; i != p->words.end; i++) {
      float v = 0.;
      feature_value(*i, p->name, v);
      v *= channel_v;

      size_t word_hash = (p->hasher(p->name[0], channel_hash)) & mask;
      feature f = {v,(uint32_t)word_hash};
      ae->sum_feat_sq[index] += v*v;
      push(ae->atomics[index], f);
    }

    if (new_index && ae->atomics[index].begin != ae->atomics[index].end)
      push(ae->indices, index);
    
    if (audit)
      {
	for (substring* i = p->words.begin+feature_offset; i != p->words.end; i++) {
	  float v = 0.;
	  feature_value(*i, p->name, v);
	  v *= channel_v;

	  size_t word_hash = (p->hasher(p->name[0], channel_hash)) & mask;
      
	  char* feature = c_string_of_substring(p->name[0]);
	  audit_data ad = {copy(base), feature, word_hash, v, true};
	  push(ae->audit_features[index], ad);
	}
	free(base);
      }
  }

  return num_chars;
}

