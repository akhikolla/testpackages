%system("R CMD Rdconv --type=html --output=/home/cipollini/local/R/library/TAQMNGR/html/TAQ-package.html /home/cipollini/Dropbox/francesco/R/MultiPlatform/TAQ/TAQ_R_Project/man/TAQ-package.Rd")

%%%% Package name
\name{TAQMNGR}

%%%% Aliases
\alias{TAQMNGR}


%\docType{package}

%%%% Title
\title{
  TAQ Manager
}

%%%% Description
\description{
  The package manages tick-by-tick transaction data, performing \emph{cleaning}, \emph{aggregation} and \emph{import}. 
}


\author{
  Francesco Calvori <francesco.calvori@gmail.com>, 
  
  Fabrizio Cipollini <cipollini.fabrizio@gmail.com>, 
  
  Giampiero M. Gallo <giampiero.gallo@gmail.com>.
  
  \strong{Maintainer}: <fabrizio.cipollini@gmail.com>
}

\details{

  The package manages tick-by-tick transaction data, performing \emph{cleaning}, \emph{aggregation} and \emph{import} in an efficient and fast way (the package engine is developed in \code{C++}).\cr
  Cleaning and Aggregation are performed according to Brownlees and Gallo (2006).

 \tabular{ll}{
 Package: \tab TAQMNGR\cr
 Type:    \tab Package\cr
 Version: \tab 2015.2-1\cr
 Date:    \tab 2015-02-21\cr
 License: \tab GPL (>= 2)\cr
 }
}

\note{
 Currently, the package processes raw data from WRDS (Wharton Research Data Service).\cr
 They have to satisfy the following requirements:\cr{}
 \itemize{%
   \item all fields have to be included (select the '\code{Check All}' button at the WRDS downloading page);\cr
   \item select the fixed-width text and '\code{G zip}' as output format and compression type, respectively, at the WRDS downloading page.
 }
An example with fake raw data can be downloaded at \url{http://local.disia.unifi.it/cipollini/webpage-new/data/data_sample.txt.gz}.\cr
The package uses the following libraries: 'Gzstream' (available at 'http://www.cs.unc.edu/Research/compgeom/gzstream/' under LGPL license), and 'zlib' (freely available at 'http://www.zlib.net/').
}

\references{
Brownlees, C. T., and Gallo, G. M. (2006). Financial Econometric Analysis at Ultra{--}High Frequency: Data Handling Concerns, \emph{Computational Statistics and Data Analysis} \bold{51}, 2232--2245.
}
