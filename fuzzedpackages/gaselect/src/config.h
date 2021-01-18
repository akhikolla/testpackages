//
//  config.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 08.05.2013.
//
//

#ifndef GenAlgPLS_config_h
#define GenAlgPLS_config_h

#include "autoconfig.h"

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#define MIN_MUTATION_PROBABILITY 0.000001 // 1e-6
#define TAB_DELIMITER "    "
#define PRECISION 8
#define WIDTH PRECISION + 5
#define DELIMITER_POSITION 4
#define BITS_PER_BYTE 8

#ifdef HAVE_UNSIGNED_LONG_LONG
	typedef unsigned long long IntChromosome;
#else
	typedef unsigned long IntChromosome;
#endif

#ifndef HAVE_UINT8_16_MAX
	typedef uint8_t unsigned short;
	typedef uint16_t unsigned short;
#endif

#if (defined HAVE_CLIMITS || defined HAVE_LIMITS_H)
	#ifdef HAVE_UNSIGNED_LONG_LONG
		#define INT_CHROMOSOME_MAX_VAL ULLONG_MAX
	#else
		#define INT_CHROMOSOME_MAX_VAL ULONG_MAX
	#endif
#endif

// Mutation algorithm
// If the ratio of set to unset bits is greater than this number, a random position
// is likely to result in a bit with the desired state
#define RATIO_RANDOM_SEARCH 0.2

#ifndef _REENTRANT
#undef HAVE_PTHREAD_H
#endif

#ifndef ENABLE_DEBUG_VERBOSITY
#define ARMA_DONT_PRINT_RUNTIME_ERRORS 1
#define ARMA_DONT_PRINT_ERRORS 1
#undef ARMA_PRINT_ERRORS
#define ARMA_NO_DEBUG 1
#endif

#define ARMA_DONT_USE_CXX11 1
#define ARMA_DONT_USE_WRAPPER 1

#endif
