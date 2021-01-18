#ifndef SCRIPT_LOC_H
#define SCRIPT_LOC_H

// Used to hold script location
struct script_loc {
	int first_line;
	int first_column;
	int last_line;
	int last_column;
};

// Init loc variable
struct script_loc loc_init();

// Show location
void loc_show( struct script_loc );

#endif // SCRIPT_LOC_H
