## R interface to CSDP semidefinite programming library
##
## By default the included CSDP code is built and linked, otherwise the user can supply
## the location of the csdp library and includes, or the location of the csdp installation directory
##
## 1) If the user sets PKG_CFLAGS, PKG_CPPFLAGS and PKG_LIBS (all three) then these are used
##
## 2) If user provides location of the csdp installation we use the lib and include directories there
##     --with-csdp=<location of csdp DIR>
## 3) If none of these are given the included CSDP code is built and linked
##
## Compiler flags for csdp can be given using --with-csdp-cflags

## require at least autoconf 2.61
AC_PREREQ(2.61)

AC_INIT([Bayesrel], 0.7.0)

FOUND_ALL=""

## Check if global variables are already set, if so we are done
if test [-n "$PKG_CFLAGS" -a -n "$PKG_CPPFLAGS" -a -n "$PKG_LIBS"] ; then
    FOUND_ALL="yes"
else
    if test [-n "$PKG_CFLAGS" -o -n "$PKG_CPPFLAGS" -o -n "$PKG_LIBS"] ; then
	AC_MSG_ERROR([All three of PKG_CFLAGS PKG_CPPFLAGS and PKG_LIBS must be given together])
    fi
fi

## get arguments
AC_ARG_WITH([csdp-cflags],
    [AS_HELP_STRING([--with-csdp-cflags],[compiler flags for csdp])],
    CSDP_CFLAGS=$withval,
    CSDP_CFLAGS="")

AC_ARG_WITH([csdp-dir],
    [AS_HELP_STRING([--with-csdp-dir],[specify location of csdp installation])],
    CSDP_DIR=$withval,
    CSDP_DIR="")

if test [-z "$FOUND_ALL" -a -n "$CSDP_DIR"] ; then
    if test [! -d "$CSDP_DIR"] ; then
	AC_MSG_ERROR([CSDP installation directory $CSDP_DIR does not exist])
    fi

    PKG_CPPFLAGS="-I$CSDP_DIR/include"
    PKG_LIBS="-L$CSDP_DIR/lib"
    FOUND_ALL="yes"
fi

## No location given so build and link provided csdp
if test [-z "$FOUND_ALL"] ; then
    echo "Using bundled CSDP"
    PKG_CPPFLAGS="-ICsdp/include"
    PKG_LIBS="-LCsdp/lib"
    CSDP_TS="Csdp.ts"
else
    CSDP_TS=""
fi

## Find the R home directory.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "Could not determine R_HOME."
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`

PKG_LIBS="${PKG_LIBS} -lsdp ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} -lm"
PKG_CFLAGS="-I${R_HOME}/R/library/RcppArmadillo/include -DNOSHORTS ${CSDP_CFLAGS}"

## add some stuff for RcppArmadillo
## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_PROG_CC
AC_PROG_CPP
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

## if given by user, check header and library
if test [-n "$FOUND_ALL"] ; then
    LIBS="${PKG_LIBS}"
    CFLAGS="${CFLAGS} ${PKG_CFLAGS}"
    CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"

    AC_CHECK_HEADER(declarations.h)
    if test "${ac_cv_header_declarations_h}" == no; then
	AC_MSG_ERROR([CSDP header declarations.h not found in provided location])
    fi

    AC_SEARCH_LIBS(easy_sdp,sdp,,AC_MSG_ERROR([Could not link CSDP library using $PKG_LIBS]))
fi

## substitute in src/Makevars.in
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CXXFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(CSDP_TS)
AC_OUTPUT(src/Makevars)

exit 0
