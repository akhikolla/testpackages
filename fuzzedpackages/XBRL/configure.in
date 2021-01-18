AC_INIT([XBRL], 0.99.17OC)

## Use Rscript to query Rcpp for compiler and linker flags
## link flag providing library as well as path to library,
## and optionally rpath
RCPP_LDFLAGS="`${R_HOME}/bin/Rscript -e 'Rcpp:::LdFlags()'`"


## Use xml2-config to find arguments for compiler + linker flags
##
## Check for non-standard programs: xml2-config(1)
AC_PATH_PROG([XML2_CONFIG], [xml2-config])
## If xml2-config was found, let's use it
if test "${XML2_CONFIG}" != ""; then
   # Use xml2-config for header and linker arguments
   #AC_MSG_NOTICE([xml2 version `${XML2_CONFIG} --version` found!])
   echo "xml2 version `${XML2_CONFIG} --version` found!"
   XML2_CPPFLAGS="`${XML2_CONFIG} --cflags`"
   XML2_LIBS="`${XML2_CONFIG} --libs`"
else
   AC_MSG_ERROR([xml2-config not found, is libxml2 installed?])
fi


# Now substitute these variables
# in src/Makevars.in to create src/Makevars
AC_SUBST(RCPP_LDFLAGS)
AC_SUBST(XML2_CPPFLAGS)
AC_SUBST(XML2_LIBS)

AC_OUTPUT(src/Makevars)
