\name{NEWS}
\title{News for package \pkg{smurf}}
\encoding{UTF-8}


\section{Changes in version 1.0.7 (2020-11-24)}{

  \subsection{Changes in \code{R}-code:}{
      \itemize{
        \item \code{plot_lambda}: Properly check whether the object \code{x} of class \code{glmsmurf} contains components related to the selection of the sparsity parameter \code{lambda}.
        \item Use named families for \pkg{glmnet} whenever possible as this is faster.
        }
  }
  
}


\section{Changes in version 1.0.6 (2020-05-16)}{

  \subsection{Changes in \code{R}-code:}{
      \itemize{
        \item Require at least version \code{4.0} of \pkg{glmnet} as this allows for more families to be used.
        }
  }
  
  \subsection{Changes in documentation:}{
      \itemize{
        \item \code{plot.glmsmurf}, \code{plot_reest.glmsmurf} and \code{plot_lambda.glmsmurf}: Link to \code{plot.default} instead of \code{plot} (\pkg{graphics}) to avoid warnings on R-4.0.0.
      }
  }
}


\section{Changes in version 1.0.5 (2020-03-19)}{

  \subsection{Changes in \code{R}-code:}{
      \itemize{
        \item Use numerical tolerance that is also suitable on platforms without a long double.
        \item Use five-fold instead of two-fold cross-validation in tests.
      }
  }
  
    \subsection{Miscellaneous changes:}{
      \itemize{
        \item Use \code{rent} dataset from \pkg{catdata} again as this package will no longer be archived.
        \item Use default date format for vignette date.
      }
  }
  
}

\section{Changes in version 1.0.4 (2020-02-09)}{

  \subsection{Changes in \code{R}-code:}{
      \itemize{
        \item Adapt test for input for prior weights to catch multiple warnings.
      }
  }
  
  \subsection{Documentation changes:}{
    \itemize{
      \item Use only one core to select the optimal value of \code{lambda} in the example of \code{plot_lambda}.
    }
  }
  
}

\section{Changes in version 1.0.3 (2020-02-07)}{

 \subsection{Changes in \code{R}-code:}{
  \itemize{
    \item Add extra tests for proximal operators, number of cross-validation folds and prior weights. 
   }
  }
  
  \subsection{Miscellaneous changes:}{
      \itemize{
        \item Add \code{rent} dataset from \pkg{catdata} as this package is scheduled to be archived on CRAN on 14 February 2020.
        \item Correct README to display correct pipeline status on GitLab.
      }
  }
  
}


\section{Changes in version 1.0.2 (2019-11-25)}{

 \subsection{Changes in \code{R}-code:}{
  \itemize{
    \item Select first element of \code{class(obj)} when \code{obj} might be a \code{matrix} object as \code{matrix} objects will also inherit from class \code{array} in R 4.0.0 (and hence \code{class(obj)} will be of length > 1).
    \item Add an extra test for the output of a \code{glmsmurf} object where a Graph-Guided Fused Lasso is used.
   }
  }
  
  \subsection{Documentation changes:}{
    \itemize{
      \item Update to \pkg{roxygen2} version 7.0.0.
      \item Add empty first line to \code{Rent_example2.R} to avoid problems with \pkg{roxygen2} version 7.0.0.
      \item Move example files from \code{/inst} to \code{/inst/examples}.
    }
  }
} 


\section{Changes in version 1.0.1 (2019-05-11)}{

  \subsection{Changes in \code{R}-code:}{
    \itemize{
      \item \code{glmsmurf}: Improve handling of coefficient names.
    }
  }
  
  \subsection{Documentation changes:}{
    \itemize{
      \item \code{plot_lambda}: Replace \code{\\dontrun} by \code{\\donttest} as requested by CRAN.
    }
  }
} 
  
\section{Changes in version 1.0.0 (2018-12-03)}{

  \subsection{Changes in \code{R}-code:}{
      \itemize{
        \item \code{p}: change order of \code{group} and \code{refcat} arguments.
        \item \code{glmsmurf}: catch errors when computing the maximum value of lambda (\href{https://gitlab.com/TReynkens/smurf/-/issues/2}{#2}).
        \item Fix bug in standardization when a continuous predictor is penalized with a Lasso or Group Lasso penalty (\href{https://gitlab.com/TReynkens/smurf/-/issues/4}{#4}).
      }
  }
  
  \subsection{Documentation changes:}{
    \itemize{
      \item Move examples for S3 methods into example for \code{glmsmurf}.
    }
  }
  
  \subsection{Miscellaneous changes:}{
      \itemize{
        \item Remove maintainer field in DESCRIPTION as it is already set using Authors@R.
        \item Change GitLab URL in README.
        \item Update Devriendt et al. (2018) reference.
        \item Add reference to Devriendt et al. (2018) in DESCRIPTION.
        \item First release on CRAN.
      }
  }
}


\section{Changes in version 0.4.1 (2018-09-25)}{

  \itemize{
    \item First public release on GitLab.
  }
  
}


\section{Changes in version 0.4.0 (2018-09-24)}{

  \subsection{Changes in \code{R}-code:}{
      \itemize{
        \item \code{glmsmurf}: Use "\code{cv1se}" to indicate selection of lambda using cross-validation with the one standard error rule. E.g. \code{"cv.dev.1se"} is renamed to \code{"cv1se.dev"}.
      }
  }
  
  \subsection{Changes in documentation:}{
      \itemize{
        \item \code{glmsmurf}: Add note that selected value of lambda for out-of-sample selection and cross-validation is not (always) deterministic.
        \item General documentation updates.
      }
  }
  
  \subsection{Miscellaneous changes:}{
      \itemize{
        \item Add continuous integration (CI) on GitLab.
        \item Add tests for \code{plot}, \code{plot_lambda} and \code{summary} functions.
        \item Add tests for elements of \code{glmsmurf}-class related to selection of lambda.
        \item General vignette update.
        \item Add LICENSE file.
      }
  }
  
}

\section{Changes in version 0.3.0 (2018-09-17)}{

  \itemize{
    \item First release on GitLab.
  }
  
}

