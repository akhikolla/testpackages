\name{NEWS}
\title{News for Package \pkg{RcppGetconf}}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\ghpr}{\href{https://github.com/eddelbuettel/rcppgetconf/pull/#1}{##1}}
\newcommand{\ghit}{\href{https://github.com/eddelbuettel/rcppgetconf/issues/#1}{##1}}

\section{Changes in inline version 0.0.3 (2018-11-16)}{
  \itemize{
    \item Examples no longer run on Solaris where they appear to fail.
  }
}

\section{Changes in inline version 0.0.2 (2016-08-01)}{
  \itemize{
    \item A new function \code{getConfig} for single values was added.
    \item The \code{struct vars} is now defined more portable allowing
    compilation on OS X (PR \ghpr{1} by Qiang Kou).
  }
}

\section{Changes in version 0.0.1 (2016-07-26)}{
  \itemize{
    \item Initial CRAN upload with \code{getAll} function.
    \item Package is functional on Linux.
  }
}
