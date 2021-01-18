\name{indelmiss-package}
\alias{indelmiss-package}
\alias{indelmiss}
\docType{package}
\title{
Insertion deletion analysis while accounting for possible missing data.
}
\description{
Genome-wide gene insertion and deletion rates can be modelled in a maximum 
    likelihood framework with the additional flexibility of modelling potential missing 
    data using the models included within. These models simultaneously estimate insertion 
    and deletion (indel) rates of gene families and proportions of "missing" data for 
    (multiple) taxa of interest. The likelihood framework is utilized for parameter 
    estimation. A phylogenetic tree of the taxa and gene presence/absence patterns 
    (with data ordered by the tips of the tree) are required. For more details, see 
    Utkarsh J. Dang, Alison M. Devault, Tatum D. Mortimer, Caitlin S. Pepperell, 
    Hendrik N. Poinar, G. Brian Golding (2016). Gene insertion deletion analysis 
    while accounting for possible missing data. Genetics (accepted).
}

\details{
\tabular{ll}{
Package: \tab indelmiss\cr
Type: \tab Package\cr
Version: \tab 1.0.7\cr
Date: \tab 2016-08-21\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@binghamton.edu}
}
\references{  
Eddelbuettel, Dirk and Romain Francois (2011). Rcpp: Seamless R and C++ Integration. Journal of Statistical Software, 40(8), 1--18.

Felsenstein, Joseph. Inferring phylogenies. Vol. 2. Sunderland: Sinauer Associates, 2004.
    
Gilbert, Paul and Ravi Varadhan (2012). numDeriv: Accurate Numerical Derivatives. R package version 2012.9-1.

Hao, Weilong, and G. Brian Golding. "The fate of laterally transferred genes: life in the fast lane to adaptation or death." Genome Research 16.5 (2006): 636--643.

Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290. R package version 3.2.
  
Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4) 592-593. R package version 1.99.11.

Wickham, Hadley (2012). stringr: Make it easier to work with strings. R package version 0.6.2.
}