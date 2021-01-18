#!/bin/sh
# Anticonf (tm) script by Jeroen Ooms & Murat Tasan (2017)
# This script will prefer cflags (specifically includefile dirs) and lib dirs
# in the following order of precedence:
#   (1) INCLUDE_DIR or LIB_DIR entered explicitly on the command line, e.g.
#       R CMD INSTALL --configure-vars='INCLUDE_DIR=/.../include LIB_DIR=/.../lib'
#   (2) Values found via 'pkg-config' for the libpq package.
#   (3) Values found via 'pg_config' given a PostgreSQL installation.

# Library settings
PKG_CONFIG_NAME="libpq"
PKG_DEB_NAME="libpq-dev"
PKG_RPM_NAME="postgresql-devel"
PKG_AMZ_RPM_NAMES="postgreql8-devel, psstgresql92-devel, postgresql93-devel, or postgresql94-devel"
PKG_CSW_NAME="postgresql_dev"
PKG_BREW_NAME="libpq"
PKG_TEST_HEADER="<libpq-fe.h>"
PKG_LIBS="-lpq"

# pkg-config values (if available)
if [ $(command -v pkg-config) ]; then
  PKGCONFIG_CFLAGS=$(pkg-config --cflags --silence-errors ${PKG_CONFIG_NAME})
  PKGCONFIG_LIBS=$(pkg-config --libs --silence-errors ${PKG_CONFIG_NAME})
  PKGCONFIG_MODVERSION=$(pkg-config --modversion --silence-errors ${PKG_CONFIG_NAME})

  # MacOS seems to ship a broken libpq.pc file
  if [[ "$OSTYPE" == "darwin"* ]]; then
    if [[ "$PKGCONFIG_CFLAGS" == *"Internal.sdk"* ]] || [[ -z "$PKGCONFIG_CFLAGS" ]]; then
      unset PKGCONFIG_CFLAGS
      unset PKGCONFIG_LIBS
    fi
  fi
fi

# pg_config values (if available)
if [ $(command -v pg_config) ]; then
  PG_INC_DIR=$(pg_config --includedir)
  PG_LIB_DIR=$(pg_config --libdir)
  PG_VERSION=$(pg_config --version)
fi

# Note that cflags may be empty in case of success
if [ "$INCLUDE_DIR" ] || [ "$LIB_DIR" ]; then
  echo "Found INCLUDE_DIR and/or LIB_DIR!"
  PKG_CFLAGS="-I$INCLUDE_DIR $PKG_CFLAGS"
  PKG_LIBS="-L$LIB_DIR $PKG_LIBS"
elif [ "$PKGCONFIG_CFLAGS" ] || [ "$PKGCONFIG_LIBS" ]; then
  echo "Found pkg-config cflags and libs ($PKG_CONFIG_NAME $PKGCONFIG_MODVERSION)!"
  PKG_CFLAGS=${PKGCONFIG_CFLAGS}
  PKG_LIBS=${PKGCONFIG_LIBS}
elif [ "$PG_INC_DIR" ] || [ "$PG_LIB_DIR" ]; then
  echo "Found pg_config includedir and libdir ($PG_VERSION)!"
  if [[ "$PG_VERSION" == "PostgreSQL 8"* ]]; then
    echo "This version of libpq is too old! We need at least 9.0!"
    exit 1
  fi
  PKG_CFLAGS="-I${PG_INC_DIR}"
  PKG_LIBS="-L${PG_LIB_DIR} ${PKG_LIBS}"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if [ $(command -v brew) ]; then
    BREWDIR=$(brew --prefix)
  else
    curl -sfL "https://jeroen.github.io/autobrew/$PKG_BREW_NAME" > autobrew
    source autobrew
  fi
  PKG_CFLAGS="-I$BREWDIR/opt/$PKG_BREW_NAME/include"
  PKG_LIBS="-L$BREWDIR/opt/$PKG_BREW_NAME/lib $PKG_LIBS"
fi

# For debugging
echo "Using PKG_CFLAGS=$PKG_CFLAGS"
echo "Using PKG_LIBS=$PKG_LIBS"

# Find compiler
CC=$(${R_HOME}/bin/R CMD config CC)
CFLAGS=$(${R_HOME}/bin/R CMD config CFLAGS)
CPPFLAGS=$(${R_HOME}/bin/R CMD config CPPFLAGS)

# Test configuration
echo "#include $PKG_TEST_HEADER" | ${CC} ${CPPFLAGS} ${PKG_CFLAGS} ${CFLAGS} -E -xc - >/dev/null 2>&1 || R_CONFIG_ERROR=1;

# Customize the error
if [ $R_CONFIG_ERROR ]; then
  echo "------------------------- ANTICONF ERROR ---------------------------"
  echo "Configuration failed because $PKG_CONFIG_NAME was not found. Try installing:"
  echo " * deb: $PKG_DEB_NAME (Debian, Ubuntu, etc)"
  echo " * rpm: $PKG_RPM_NAME (Fedora, EPEL)"
  echo " * rpm: $PKG_AMZ_RPM_NAMES (Amazon Linux)"
  echo " * csw: $PKG_CSW_NAME (Solaris)"
  echo " * brew: $PKG_BREW_NAME (OSX)"
  echo "If $PKG_CONFIG_NAME is already installed, check that either:"
  echo "(i)  'pkg-config' is in your PATH AND PKG_CONFIG_PATH contains"
  echo "     a $PKG_CONFIG_NAME.pc file; or"
  echo "(ii) 'pg_config' is in your PATH."
  echo "If neither can detect $PGK_CONFIG_NAME, you can set INCLUDE_DIR"
  echo "and LIB_DIR manually via:"
  echo "R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'"
  echo "--------------------------------------------------------------------"
  exit 1;
fi

# Write to Makevars
sed -e "s|@cflags@|$PKG_CFLAGS|" -e "s|@libs@|$PKG_LIBS|" src/Makevars.in > src/Makevars.new
if [ ! -f src/Makevars ] || (which diff > /dev/null && ! diff -q src/Makevars src/Makevars.new); then
  cp -f src/Makevars.new src/Makevars
fi
rm -f src/Makevars.new

# Success
exit 0
