
#include <R_ext/Print.h>

/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This fconvergence is part of software Realea
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
#include "debug.h"
#include <cstdio>
#include <stdarg.h>
#include <cstdlib>
#include <string>

void print_error(const char *error_str, ...) {
   va_list ap;
   va_start (ap, error_str);
   REvprintf( error_str, ap);
   va_end (ap);
}


static bool info=false;
static bool debug=false;

void enable_print_info(void) {
    info = true;
}

void disable_print_info(void) {
    info = false;
}

void enable_print_debug(void) {
    debug = true;
}

void disable_print_debug(void) {
    debug = false;
}

void print_info(const char *output_str, ...) {
   va_list ap;

   if (info) {
    va_start (ap, output_str);
    Rvprintf(output_str, ap);
    va_end (ap);
   }
}

void print_output(const char *output_str, ...) {
   va_list ap;

   va_start (ap, output_str);
   Rvprintf(output_str, ap);
   va_end (ap);
}




void print_debug(const char *output_str, ...) {
   va_list ap;

   if (debug) {
    va_start (ap, output_str);
    Rvprintf (output_str, ap);
    va_end (ap);
   }
}
static int countFitness=0;
static string templatefname;
static double best;
static FILE *fconvergence;

void init_output_convergence(const char *fname) {
    countFitness = 0;
    templatefname = fname;
}

static bool foutput=false;

void open_output_convergence(int num) {
   char fname[300];

   if (templatefname == "output") {
      foutput = true;
      //fconvergence = stdout;
      fconvergence = NULL;
   }
   else if (templatefname != "") {
	sprintf(fname, "%s_%d.dat", templatefname.c_str(), num);
	fconvergence = fopen(fname, "w");
   }

   if (fconvergence == NULL) {
      print_error("Error escribiendo en el fichero %s\n", fname);
//      exit(1);
   }
}

void close_output_convergence(void) {
    if (templatefname != "output") {
        fclose(fconvergence);
    }
    countFitness = 0;
}

void print_convergence(const tGen *sol, int dim, tFitness fitness) {
    countFitness++; 
	
    if ( (countFitness == 1) || (fitness < best && (countFitness < 10000*dim)))  {
       best = fitness;

       if (foutput) {
          print_output("%d %.30Lf\n", countFitness, fitness);
       }
       else {
          fprintf(fconvergence,"%d %.30f\n", countFitness, fitness);
          fflush(fconvergence);
       }
    }   
}
