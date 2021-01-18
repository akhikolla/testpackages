#include <R_ext/Print.h>
#ifndef SAILR_HELPER_H
#define SAILR_HELPER_H

#include "stdio.h"
#include "stdbool.h"

// DEBUG_PRINT macro

#ifdef DEBUG
#define DEBUG_PRINT(...) do{ fprintf( stderr, __VA_ARGS__); } while(false)
#else
#define DEBUG_PRINT(...) do{ } while( false )
#endif

#ifdef DEBUG
#define DEBUG_FLUSH() do{ fflush(stdout); } while(false)
#else
#define DEBUG_FLUSH() do{ } while( false )
#endif

// Default Encoding
// "UTF8" or "LATIN1" is supported.

#define SAILR_DEFAULT_ENCODING  "UTF8"
#define SAILR_DEFAULT_REXP_ENCODING  SAILR_DEFAULT_ENCODING
#define SAILR_DEFAULT_VM_CHARACTER_ENCODING  SAILR_DEFAULT_ENCODING

#endif /* SAILR_HELPER_H */
