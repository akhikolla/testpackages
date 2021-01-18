\name{NEWS}
\title{blackbox News}
\encoding{UTF-8}
\section{CHANGES UP TO VERSION 1.1.29}{
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Checks of invalid input (and resulting screen output) adapted to non-Migraine usage.
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item 'lbfgsb3' is replaced by 'lbfgsb3c' as a suggested package.
\item Requests a more recent (>= 3.1.0) version of 'spaMM'.
}
}
}
\section{CHANGES UP TO VERSION 1.1.27}{
\subsection{CORRECTED BUGS}{
\itemize{
\item Long-standing glitch that substantially slowed blackbox down (though with correct results).
}
}
}
\section{CHANGES UP TO VERSION 1.1.25}{
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item (Migraine-specific) New composite parameters and profile computations for migration parameters.
\item Dependency on 'ff' package removed, assuming a recent version of spaMM (here >= 2.4.8).
\item Dependencies changed for parallel computations. See help("blackbox.options") for details.
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item (Migraine-specific) writeFinalInfo() now correctly writes the results when sampleByResp() has not been run.
\item (Migraine-specific) CKrigcoefs() could fail when minSmoothness=maxSmoothness.
\item (Migraine-specific) user-given 'CovFnParam' were ignored.
\item (Migraine-specific) now computes Dgmu parameter only when Dg/2N is a parameter of the model (e.g. not a for a discrete change in pop size).
\item (Migraine-specific) resetCHull() handles input points in linear space of lower dimension than the number of fitted parameters, now avoiding a bug in the rare case where there are many profile points, all of them in such a space. 
}
}
} 
\section{CHANGES UP TO VERSION 1.1.7}{
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item No user-level changes.
}
}
\subsection{OTHER CHANGES}{
\itemize{
\item Changes for compatibility wth Rcpp::compileAttributes() in Rcpp >= 0.12.12
\item Other internal changes for compatibility with planned changes in other packages.
}
}
}
\section{CHANGES UP TO VERSION 1.1.2}{
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item (Migraine-specific) Parallel computation is implemented for profile point computations.
}
}
}
\section{CHANGES UP TO VERSION 1.1.0}{
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Less primitive convergence diagnostics are implemented in bboptim().
\item (Migraine-specific) New mesures of ''goodness of prediction'' implemented.
\item (Migraine-specific) add logscale for 1D profile plots
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item blackbox.options(<new values>) did not return old values.
\item (Migraine-specific) Return code modified for consistency with documentation (but this does not really affect the analyses).
\item (Migraine-specific) variables names could be lost in generateNewPoints().
\item (Migraine-specific) Invalid requests by users for CIs for some composite variables could lead to an execution error.
}
}
\subsection{OTHER CHANGES}{
\itemize{
\item Added native routine registration.
}
}
}
\section{CHANGES UP TO VERSION 1.0.18}{
\subsection{BUGS}{
\itemize{
\item More C++ portability issues fixed.
\item (Migraine-specific, without any impact on inferences) optimization convergence diagnostic (actually nloptr()'s return status) was missing from the output_...txt file, resulting in lines with different number of values.
\item (Migraine-specific) toCanonical() corrected to handle single-candidate input matrix.
\item (Migraine-specific) fixed parameter space issues in generateNewPoints() and canonize().
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Modified init values of some smoothing parameters may change bboptim() results for rough response surfaces.
}
}
\subsection{OTHER CHANGES}{
\itemize{
\item Some function definitions has been moved from the spaMM package to blackbox.
}
}
}
\section{CHANGES UP TO VERSION 1.0.8}{
\subsection{BUGS}{
\itemize{
\item (Migraine-specific) Wrong variable spaces in generateNewPoints() and toCanonical() corrected.
\item (Migraine-specific) Calls to Rf_error() removed from Rcpp code.
\item (Migraine-specific) An invalid cbind() that affected the generation of 'nextpoints' has been corrected.
\item (Migraine-specific) Within subHullWrapper(), the scdd() results with double precision are now checked and replaced by a computation with rational arithmetic if they are inaccurate.  
\item (Migraine-specific) "samplingSpace=,Nb," has been broken for a long time. This is corrected.
\item (Migraine-specific) sampleNextPars() could stop when called with non-default 'samplingSpace'. This is corrected.
}
}
}


