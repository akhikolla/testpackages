\name{NEWS}
\title{genepop News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
\section{CHANGES UP TO VERSION 1.1.7}{
\subsection{CORRECTED BUGS}{
\itemize{
\item Incorrect tables of multi-locus and multi 'pop' tests of Hardy-Weinberg [test_HW() which default which="Proba"] when some individual tests had p-value estimates equal to 0, similar to the problem fixed in version 1.1.4 for differentiation tests. Per-locus per-pair-of-'pop' tests were correctly reported. More display glitches fixed in this and other options.
}
}
}
\section{CHANGES UP TO VERSION 1.1.4}{
\subsection{CORRECTED BUGS}{
\itemize{
\item Incorrect final table of multi-locus tests of pairwise differentiation [test_diff(...,pairs = TRUE)] particularly in cases where all individual tests had p-value estimates equal to 0 (per-locus per-pair tests were correctly reported). Several minor display glitches in similar tables have also been fixed (as a result, scientific format for numbers may now appear in such outputs).
}
}
}
\section{CHANGES UP TO VERSION 1.1.3}{
\subsection{NEW FEATURES}{
\itemize{
\item New convenience function genepopExample() replacing potentially confusing code in the Examples
}
}
}
\section{CHANGES UP TO VERSION 1.1.2}{
\subsection{NEW FEATURES}{
\itemize{
\item New keyword `intra_all_types` for setting "popTypeSelection" allows one to perform a single spatial regression for all pairs of individuals or populations within several types of subpopulations.
\item The IBD bootstrap no longer writes many temporary files. This removes an exposure to operating system-dependent crashes. 
\item (R package only) Progress bar for some long computations.
\item New function clean_workdir() to remove easily files written by Genepop (including output files!).
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item Missing pairwise genetic information was no correctly detected in Mantel tests. 
\item Fixed documentation bug referring to PopTypeSelection="inter_all_pairs" instead of the correct PopTypeSelection="inter_all_types".
}
}
}
\section{CHANGES UP TO VERSION 1.0.5}{
\subsection{NEW FEATURES}{
\itemize{
\item The name extension and accessibility of the genepop input file is now checked within the R code.
\item Existence of the output file is checked.
\item More accurate messages when differentiation options are run on an input file containing only one 'pop'.
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item Linkage disequilibrium tests incorrectly indicated 'No data' for tests involving the last locus.
\item The ibd() last line of output in the <>.IBD file mis-reported an uninteresting computation. Both the computation and the line of output have been removed.
\item (R package only) Mantel test always used rank correlation statistic.
}
}
}

