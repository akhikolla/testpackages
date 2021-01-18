AC_INIT(SplitReg, 1.0.1)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

AC_PROG_CXX

AC_LANG(C++)

AC_CHECK_HEADERS(stdint.h inttypes.h limits.h climits)

## ----------------------------------------------------------------
# Check if unsigned long long is supported (it's not standard yet)
## ----------------------------------------------------------------

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  unsigned long long x = 100LL; return 0; }]])
	],
	[
		AC_MSG_RESULT([Type unsigned long long is supported])
		AC_DEFINE([HAVE_UNSIGNED_LONG_LONG], [1], [Define to 1 if unsigned long long is supported])
	],
	[AC_MSG_RESULT([Type unsigned long long is not supported])]
)

## -----------------------------------------------------------
# Check if support for uint16_t and uint8_t is complete
## -----------------------------------------------------------

if test "$ac_cv_header_stdint_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_PROGRAM([[#define __STDC_LIMIT_MACROS 1
		#include <stdint.h>]],
					[[uint16_t x = UINT16_MAX; uint8_t y = UINT8_MAX;]])
	],
	[
		AC_MSG_RESULT([uint16_t and unit8_t are available in stdint.h])
		AC_DEFINE([HAVE_UINT8_16_MAX], [1], [Define to 1 if types uint16_t and uint8_t are in stdint.h])
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint16_t and uint8_t are not available in stdint.h])]
	)
fi
if test "$ac_cv_header_inttypes_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_PROGRAM([[#define __STDC_LIMIT_MACROS 1
		#include <inttypes.h>]],
					[[uint16_t x = UINT16_MAX; uint8_t y = UINT8_MAX;]])
	],
	[
		AC_MSG_RESULT([uint16_t and unit8_t are available in inttypes.h])
		AC_DEFINE([HAVE_UINT8_16_MAX], [1], [Define to 1 if types uint16_t and uint8_t are in inttypes.h])
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint16_t and uint8_t are not available in inttypes.h])]
	)
fi

# if test "$suitable_integer" != "yes"; then
# 	AC_MSG_FAILURE([no suitable integral type available])
# fi;

## -----------------------------------------------------------
# Check for builtin function to "count trailing zeros" (ctz)
## -----------------------------------------------------------

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctzl(2UL); return 0; }]])
	],
	[
		AC_MSG_RESULT([Found builtin CTZ method (for unsigned long)])
		AC_DEFINE([HAVE_GCC_CTZL], [1], [Define to 1 if you have function __builtin_ctzl(unsigned long)])
	],
	[AC_MSG_RESULT([No builtin CTZ method (for unsigned long)])]
)

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctzll(2ULL); return 0; }]])
	],
	[
		AC_MSG_RESULT([Found builtin CTZ method (for unsigned long long)])
		AC_DEFINE([HAVE_GCC_CTZLL], [1], [Define to 1 if you have function __builtin_ctzll(unsigned long long)])
	],
	[AC_MSG_RESULT([No builtin CTZ method (for unsigned long long)])]
)

AC_CHECK_FUNCS(ffsl)
AC_CHECK_FUNCS(ffsll)

## -----------------------------------------------------------
# Check for builtin function to "count set bits" (popcountl)
## -----------------------------------------------------------

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountl(10UL); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method (for unsigned long)])
		AC_DEFINE([HAVE_BUILTIN_POPCOUNTL], [1], [Define to 1 if you have function __builtin_popcountl(unsigned long)])
	],
	[AC_MSG_RESULT([No builtin popcount method (for unsigned long)])]
)

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountll(10ULL); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method (for unsigned long long)])
		AC_DEFINE([HAVE_BUILTIN_POPCOUNTLL], [1], [Define to 1 if you have function __builtin_popcountl(unsigned long long)])
	],
	[AC_MSG_RESULT([No builtin popcount method (for unsigned long long)])]
)

## --------------------------
# Check for pthreads support
## --------------------------

oldCFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -pthread"
AC_CHECK_HEADERS(pthread.h)

CFLAGS="$oldCFLAGS"

AC_SUBST(CXX11FLAGS)

AC_CONFIG_HEADERS([src/autoconfig.h])

## --------------------------
## Check for OpenMP support
## --------------------------

dnl this the meat of R's m4/openmp.m4
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [AC_LINK_IFELSE([_AC_LANG_OPENMP],
   	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
	  for ac_option in -fopenmp -xopenmp -qopenmp \
                           -openmp -mp -omp -qsmp=omp -homp \
			   -fopenmp=libomp \
                           -Popenmp --openmp; do
	    ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
	    _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
	    AC_LINK_IFELSE([_AC_LANG_OPENMP],
	      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
	    _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
	    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
	      break
	    fi
	  done])])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
      "none needed" | unsupported)
	;; #(
      *)
	OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
    esac
  fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT