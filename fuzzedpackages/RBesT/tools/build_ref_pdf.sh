#!/bin/bash -x

[ -z "$R_HOME" ] && R_HOME=`R RHOME`

## create documentation reference PDF
install -d inst/doc
"${R_HOME}/bin/R" CMD Rd2pdf --batch --no-preview --force --output=inst/doc/RBesT.pdf .

## make PDF small
"${R_HOME}/bin/R" --vanilla --slave -e 'library(tools); tools::compactPDF("inst/doc/RBesT.pdf")'

