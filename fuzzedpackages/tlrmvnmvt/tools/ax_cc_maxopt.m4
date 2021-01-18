# ===========================================================================
#       https://www.gnu.org/software/autoconf-archive/ax_cc_maxopt.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CC_MAXOPT
#
# DESCRIPTION
#
#   Try to turn on "good" C optimization flags for various compilers and
#   architectures, for some definition of "good". (In our case, good for
#   FFTW and hopefully for other scientific codes. Modify as needed.)
#
#   The user can override the flags by setting the OPTFLAG environment
#   variable. The user can also specify --enable-portable-binary in order to
#   disable any optimization flags that might result in a binary that only
#   runs on the host architecture.
#
#   Note also that the flags assume that ANSI C aliasing rules are followed
#   by the code (e.g. for gcc's -fstrict-aliasing), and that floating-point
#   computations can be re-ordered as needed.
#
#   Requires macros: AX_CHECK_COMPILE_FLAG, AX_COMPILER_VENDOR,
#   AX_GCC_ARCHFLAG, AX_GCC_X86_CPUID.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 18

AC_DEFUN([AX_CC_MAXOPT],
[
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AX_COMPILER_VENDOR])

AC_ARG_ENABLE(portable-binary, [AS_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])],
	acx_maxopt_portable=$enableval, acx_maxopt_portable=no)

# Try to determine "good" native compiler flags if none specified via OPTFLAG
  case $ax_cv_cxx_compiler_vendor in
    dec) OPTFLAG="$OPTFLAG -newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
	 if test "x$acx_maxopt_portable" = xno; then
           OPTFLAG="$OPTFLAG -arch host"
         fi;;

    sun) OPTFLAG="$OPTFLAG -native -fast -xO5 -dalign"
	 if test "x$acx_maxopt_portable" = xyes; then
	   OPTFLAG="$OPTFLAG -xarch=generic"
         fi;;

    hp)  OPTFLAG="$OPTFLAG +Oall +Optrs_ansi +DSnative"
	 if test "x$acx_maxopt_portable" = xyes; then
	   OPTFLAG="$OPTFLAG +DAportable"
	 fi;;

    intel) OPTFLAG="$OPTFLAG -O3 -ansi_alias"
	;;

    gnu)
     # default optimization flags for gcc on all systems
     OPTFLAG="$OPTFLAG -O3 -fomit-frame-pointer"
     ;;

    microsoft)
     # default optimization flags for MSVC opt builds
     OPTFLAG="$OPTFLAG -O2"
     ;;
  esac

  if test -z "$OPTFLAG"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best OPTFLAG for this system  *"
        echo "* Use ./configure OPTFLAG=... to specify your own flags *"
	echo "* (otherwise, a default of OPTFLAG=-O3 will be used)    *"
	echo "********************************************************"
	echo ""
        OPTFLAG="$OPTFLAG -O3"
  fi

  AX_CHECK_COMPILE_FLAG($OPTFLAG, [], [
	echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed OPTFLAG don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure OPTFLAG=... to specify your own flags *"
        echo "********************************************************"
        echo ""
  ])

])
