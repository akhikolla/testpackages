#include <R_ext/Print.h>
#include "script_loc.h"

#include <stdio.h>

// Init loc variable
struct script_loc
loc_init()
{
	struct script_loc loc;
	loc.first_line = 0;
	loc.first_column = 0;
	loc.last_line = 0;
	loc.last_column = 0;
	return loc;
}

// Show location
void
loc_show( struct script_loc loc)
{
	Rprintf("approximate script position: from line %d col %d to line %d col %d \n", loc.first_line, loc.first_column, loc.last_line, loc.last_column );
}

