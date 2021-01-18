#!/bin/bash

R_HOME=${1-/CHBS/apps/R/3.2.3}

RPKG=OncoBayes2

# by default we deploy the latest stable release which is under trunk,
# to deploy a specific tag and not the latest stable, the user can
# specify this as second argument, for example, tags/REL-1.0-0
SVN_URL=${2-trunk}
SVN_BASE=http://chbslx0132.eu.novartis.net/svn/OncoBayes2/  # Adapt to your directory structure!
URL=$SVN_BASE/$SVN_URL

REPORT_DIR=/tmp/${RPKG}_install
R_CMD=$R_HOME/bin/R
WORK=$(mktemp -dt "$(basename $0).XXXXXXXXXX")

if [ -d $REPORT_DIR ]
then
    echo "Please remove $REPORT_DIR from an old installation of ${RPKG} first."
    exit 1
fi

mkdir $REPORT_DIR

echo Installing ${RPKG} from:           $URL
echo Using R_HOME:                    $R_HOME
echo Saving validation material in:   $REPORT_DIR
echo Running installation routine in: $WORK

pushd $WORK

# checkout sources
svn co $URL ${RPKG}_source 2>&1

# extract version
RPKG_VERSION=`grep Version ${RPKG}_source/DESCRIPTION | cut -d : -f 2 | tr -d [[:space:]]`

SVN_REVISION=`svnversion -n ${RPKG}_source`

echo Installing ${RPKG} version:        ${RPKG_VERSION}, SVN revision ${SVN_REVISION}

# build source package
$R_CMD CMD build ${RPKG}_source 2>&1

SOURCE_TAR=`ls ${RPKG}*.tar.gz`

echo Source package created         : ${SOURCE_TAR}

# source tar copied for convenience, not needed to store in PROTON
cp -v $SOURCE_TAR $REPORT_DIR

# tests are run later as OQ, run only IQ related matters here
$R_CMD CMD check --no-tests $SOURCE_TAR 2>&1

# install it with tests
$R_CMD CMD INSTALL --install-tests $SOURCE_TAR 2>&1

# IQ done

# now copy documents needed for PROTON into report directory

# vignette
echo Copying vignette to ${REPORT_DIR}
$R_CMD --silent --vanilla -e "file.copy(system.file(\"doc/introduction.html\", package=\"${RPKG}\"), \"${REPORT_DIR}\")"
# reference manual
echo Copying reference manual to ${REPORT_DIR}
$R_CMD --silent --vanilla -e "file.copy(system.file(\"doc/${RPKG}.pdf\", package=\"${RPKG}\"), \"${REPORT_DIR}\")"

popd

echo "R_HOME=${R_HOME}" > $REPORT_DIR/R_HOME

# delete temporary directory $WORK
rm -rf $WORK

echo "Please find all validation related documents for PROTON under $REPORT_DIR"
