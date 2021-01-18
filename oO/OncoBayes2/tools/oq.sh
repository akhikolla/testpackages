#!/bin/bash

RPKG=OncoBayes2

REPORT_DIR=/tmp/${RPKG}_install

if [ ! -d $REPORT_DIR ]
then
    echo "Please first run install.sh to install ${RPKG}."
    exit 1
fi

# import R_HOME defined during installation
. $REPORT_DIR/R_HOME

LOG=$REPORT_DIR/${RPKG}-OQ.log

R_CMD=$R_HOME/bin/R

echo Using R_HOME: $R_HOME
echo Running OQ tests for ${RPKG}...

$R_CMD --silent --vanilla -e "source(system.file(\"extra/run-oq.R\", package=\"${RPKG}\"))" > $LOG 2>&1

echo OQ log file created by R:

cat $LOG

# Any test which fails will stop R immediatley and hence prevent
# further execution. All tests were successful if "ALL TESTS
# SUCCESSFUL" is displayed at the end of the log file.

echo "OQ log file in $LOG"

exit 0
