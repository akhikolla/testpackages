AC_INIT([tlrmvnmvt], [1.0], [jian.cao@kaust.edu.sa])


F77=`$R_HOME/bin/R CMD config F77`
if test "$F77" = "" ; then
	AC_MSG_ERROR("F77 compiler not found")
fi


AC_LANG(C++)


m4_include([tools/ax_compiler_vendor.m4])
m4_include([tools/ax_check_compile_flag.m4])
m4_include([tools/ax_cc_maxopt.m4])
CXX=`$R_HOME/bin/R CMD config CXX14`
if test "$CXX" = "" ; then
	AC_MSG_ERROR("CXX14 compiler not found")
fi
AX_CC_MAXOPT()


AC_SUBST(OPTFLAG)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
