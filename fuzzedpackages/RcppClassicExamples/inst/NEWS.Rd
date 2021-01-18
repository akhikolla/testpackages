\name{NEWS}
\title{News for Package \pkg{RcppClassicExamples}}
\newcommand{\cpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.1.2 (2018-03-15)}{
  \itemize{
    \item Registered S3 print method [per CRAN request]
    \item Added \code{src/init.c} with registration and updated all
    \code{.Call} usages taking advantage of it
    \item Updated http references to https
    \item Updated DESCRIPTION conventions
  }
}

\section{Changes in version 0.1.1 (2012-12-30)}{
  \itemize{
    \item Added versioned Depends: in DESCRIPTION to not build under
    older versions of \cpkg{Rcpp} and \cpkg{RcppClassic}
  }
}

\section{Changes in version 0.1.0 (2012-12-27)}{
  \itemize{
    \item Initial release, split from package \cpkg{RcppExamples}
  }
}