\name{NEWS} 
\title{News for Package 'hsar'}

\section{Changes in version 0.5.1 (2020-6-27)}{
 \itemize{
 \item Small changes to get rid of Warnings.
 }
}

\section{Changes in version 0.5 (2020-6-1)}{
 \itemize{
 \item The draws for the coefficients are returned.
 \item The summary has been improved. The SD and the quantiles are printed.
  \item spatialreg is used instead of depreciated functions of spdep.
  \item Correct a bug, making the hsar(lambda=0) routine to crash.
  \item A new vignette with real estate dataset from Athens.
  \item A new dataset with real estate properties in municipality of Athens for 2017.
  \item Minor bugs noticed.
 }
}

\section{Changes in version 0.4.2 (2017-10-16)}{
 \itemize{
  \item Acknowledgement of our financial support (in the description file).
  \item Minor but important compatibility changes. 
 }
}

\section{Changes in version 0.4.0 (2016-04-20)}{
  \itemize{
     \item Documentation has been added as vignette that explains, how you can use the hsar() with synthetic data created for that purpose
     \item The dataframe landprice has been updated and the SpatialPointsDataFrame landSPDF and the SpatialPolygonsDataFrame Beijingdistricts have been added.These are used in the function's examples.
     \item Formulae of the form y~x1+x2+x3 and basic data checks have been implemented
     \item The examples in documentation have been improved
     \item New starting parameters available to users
     \item A thin parameter has been inserted
     \item Bug fix regarding the burnin and NSim parameters 
     \item The parameter Unum in HSAR model and in the data has been removed
     \item The z definition in documentation is corrected
     \item The output has been improved and the print functions have been added
  }
}
