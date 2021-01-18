#!/bin/bash -x

[ -z "$R_HOME" ] && R_HOME=`R RHOME`

## create documentation reference PDF
install -d inst/doc
"${R_HOME}/bin/R" CMD Rd2pdf --batch --no-preview --force --output=inst/doc/OncoBayes2.pdf .

## create SBC report
"${R_HOME}/bin/R" --slave -e "library(rmarkdown); setwd('inst/sbc/'); rmarkdown::render('sbc_report.R')"

## make PDF small
"${R_HOME}/bin/R" --vanilla --slave -e "library(tools); tools::compactPDF('inst/doc/OncoBayes2.pdf')"

