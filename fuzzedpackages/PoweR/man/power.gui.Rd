\name{power.gui}
\alias{power.gui}
\title{PoweR GUI}
\description{
Graphical user interface (GUI) for the package.
}
\usage{
power.gui()
}
\details{
  This GUI is a 5-tabbed notebook whose goal is to make our package
  easier to use :

- Tab 1 \code{\link{gensample}} : generate random samples from a law added in the package;

- Tab 2 \code{\link{statcompute}} : perform the test for a given index
  value of test statistic;

- Tab 3 \code{\link{many.crit}} : computation of critical values for several test statistics;

- Tab 4 \code{\link{powcomp.fast}} : computation of power and level tables for hypothesis tests;

- Tab 5 \emph{Examples} : reproduce results from published articles.

Important note concerning 'Iwidgets': for the GUI to work, a third party software has to be installed.

Under Microsoft Windows:

First, install ActiveTcl following indications given here:
'http://www.sciviews.org/_rgui/tcltk/TabbedNotebook.html'

After the installation of ActiveTcl and the modification of the PATH
variable, launch from an MsDOS terminal (accessible through typing 'cmd'
in the Start Menu) the following command:
C:\\Tcl\\bin\\teacup.exe install Iwidgets

You can then check the existence of a directory called 'Iwidgets4.0.2'
in 'C:\\Tcl\\lib\\teapot\\package\\tcl\\lib'.

Under Linux:

Install 'iwidgets'.


}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
