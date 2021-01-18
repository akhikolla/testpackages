### configure.ac					-*- Autoconf -*-

dnl Licensed under GPL-3 with the Autoconf Configure
dnl Script Exception, version 3.0.

AC_PREREQ(2.62)
AC_INIT([icr],[0.5.6],[],[icr],[])

dnl Ensure C/C++ compiler and compiler flags are the same in configure test
dnl and when compiling the package (see section 1.2 of 'Writing R Extensions')
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

dnl Test for OpenMP
AC_LANG(C++)
AC_OPENMP

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
