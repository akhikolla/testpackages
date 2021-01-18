\name{scatterPlot}
\alias{scatterPlot}
\title{Draw Scatter Plots of Edited Dataset}
\description{
Compare scatter plots of the input dataset and the edited dataset. 
}
\usage{
scatterPlot(model.obj, data.obj=NULL, xvar=NULL, yvar=NULL)
}
\arguments{
	\item{model.obj}{\code{\link{Rcpp_bei}} model-object generated from \code{\link{createModel}}.}
	\item{data.obj}{\code{\link{EditIn.data}} object generated from \code{\link{readData}}. This is used to draw the lines to indicate ratio edits.}
	\item{xvar}{variable to draw on X-axis. If \code{NULL}, plots are drawn for all variables in turn.}
	\item{yvar}{variable to draw on Y-axis. If \code{NULL}, plots are drawn for all variables in turn.}
}
\details{
Draw the scatter plots of log transformed values of the input data (\code{model.obj$Y.input}, left panel) and the edit-imputed data (\code{model.obj$Y.edit}, right panel). The sky-blue dots on the background represent edit-passing records which are identical for both datasets. The blue dots on the left panel shows the edit-failing records in the original data and those on the right panel shows their edited values. The red dotted lines (if any) show the ratio edits of two log-transformed variables.
}