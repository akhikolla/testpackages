// This is a function which does nothing with examples.  Used when VW is used as a compressor.

#include <cstring>		// bzero()
#include "example.h"
#include "parser.h"
#include "gd.h"

void start_noop()
{
  example* ec = NULL;
  
  while ( !parser_done()){
    ec = get_example(0);
    if (ec != NULL)
      finish_example(ec);
  }
}

void end_noop()
{
}
