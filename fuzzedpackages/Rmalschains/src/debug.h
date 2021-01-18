/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "stdarg.h"
#include "real.h"

/**
 * Enable the info messages
 */
void enable_print_info(void);
/**
 * Disable the info messages
 */
void disable_print_info(void);
/**
 * Print info in printf format if is was enable
 * @param output_str string in printf format
 */
void print_info(const char *output_str, ...);

/**
 * Print error msg in error output with printf format
 * @param error_str string in printf format
 */
void print_error(const char *error_str, ...);


/**
 * Enable the debug messages
 */
void enable_print_debug(void);
/**
 * Disable the debug messages
 */
void disable_print_debug(void);
/**
 * Print debug in printf format if is was enable
 * @param output_str string in printf format
 */
void print_debug(const char *output_str, ...);

void init_output_convergence(const char *fname);

void open_output_convergence(int num);
void print_convergence(const tGen *sol, int dim, tFitness fitness); 
void close_output_convergence(void);

