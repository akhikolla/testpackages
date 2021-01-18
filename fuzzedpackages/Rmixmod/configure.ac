#https://nerdland.net/2009/07/detecting-c-libraries-with-autotools/
AC_INIT(Rmixmod, 2.1.4)
AC_PATH_PROG([PKG_CONFIG], [pkg-config])
if test "${RMIXMOD_WITH_XML}" = "1"  ; then
	if test -z "${PKG_CONFIG}"; then
		AC_MSG_WARN([cannot find pkg-config, Rmixmod will be installed without optional XML features])
	else
		PKG_CHECK_MODULES([LIBXMLXX], [libxml++-2.6],[],[AC_MSG_WARN([No package 'libxml++-2.6' found, Rmixmod will be installed without optional XML features])])
	fi
fi
echo "LIBXMLXX_CFLAGS: $LIBXMLXX_CFLAGS"
echo "LIBXMLXX_LIBS: $LIBXMLXX_LIBS"
RMIXMOD_XML=""
if test ! -z "$LIBXMLXX_CFLAGS"  ; then
	RMIXMOD_XML="-DRMIXMOD_XML=1"
fi
OMP_FLAGS=""
if test "${RMIXMOD_WITH_OMP}" = "1"  ; then
	OMP_FLAGS="-fopenmp"
fi
AC_CHECK_HEADER([omp.h], [AC_SUBST(OMP_STATUS,"${OMP_FLAGS}")], [AC_SUBST(OMP_STATUS,"-DWITHOUT_OMP")], [])
AC_SUBST(LIBXMLXX_CFLAGS)
AC_SUBST(LIBXMLXX_LIBS)
AC_SUBST(RMIXMOD_XML)
AC_OUTPUT(src/Makevars)
