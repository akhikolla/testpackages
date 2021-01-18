dnl To recreate this configure.ac file and produce a configure file:
dnl  call autoscan
dnl  use the code from ARTP2's configure.ac file to create a system-dependent Makevars file
dnl Additional code from https://stackoverflow.com/a/42281495.
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([AnaCoDa], m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]), [cedric.landerer@gmail.com])

## Determine install location of R
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

RBIN="${R_HOME}/bin/R"
CXX11=`"${RBIN}" CMD config CXX11`
CXX11STD=`"${RBIN}" CMD config CXX11STD`
CXX="${CXX11} ${CXX11STD}"
CXXFLAGS=`"${RBIN}" CMD config CXX11FLAGS`

## Package requires C++
AC_LANG(C++)
AC_REQUIRE_CPP

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

dnl this the meat of R's m4/openmp.m4
dnl copied from ARTP2's configure.ac file
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

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
if test "${ac_cv_header_stdlib_h}" = no; then
    AC_MSG_ERROR("stdlib.h not found")
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
