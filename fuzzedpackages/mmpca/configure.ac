AC_INIT(mmpca, 1.0)
AC_LANG(C++)
AC_PROG_CXX

AX_CXX_COMPILE_STDCXX_11()

dnl Following lines kindly supplied by Dirk Eddelbuettel
dnl Use gsl-config to find arguments for compiler and linker flags

dnl Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
dnl If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   dnl Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

dnl If GSL is in a non-standard location we have to append it to CXXFLAGS
dnl before compiling the snippet below
CXXFLAGS_backup="${CXXFLAGS}"
CXXFLAGS="${GSL_CFLAGS} ${CXXFLAGS}"
AC_MSG_CHECKING([if GSL version >= 1.12])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit (major == 1 && minor < 12);
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.12])
else
  AC_MSG_RESULT([yes])
fi
dnl Reset CXXFLAGS
CXXFLAGS="${CXXFLAGS_backup}"

dnl Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_OUTPUT(src/Makevars)
