\name{plot.ksIRT}
\alias{plot.ksIRT}
\alias{plot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Method for ksIRT - kernel smoothing in Item Response Theory
}
\description{
The plot method for \code{ksIRT} objects includes a variety of exploratory plotting tools.
}
\usage{
\method{plot}{ksIRT}(x, plottype = c("OCC", "EIS", "density", "expected", "sd", 
"triangle", "tetrahedron", "RCC", "EISDIF", "OCCDIF", "PCA", "expectedDIF", 
"densityDIF"), items = "all", subjects, axistype = c("scores", "distribution"), 
alpha, main, xlab, ylab, xlim, ylim, cex, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
\item{x}{a \code{ksIRT} created with \code{ksIRT()}}

\item{plottype}{
the type of plot to be created (see section Details below). 
With the default value, \code{plottype="OCC"}, option characteristic curves are plotted.
}

\item{items}{a vector containing the items to be plotted. 
With the default value, \code{items="all"}, all items are plotted.}

\item{alpha}{either \code{alpha=FALSE} or a scalar indicating the confidence level to be used for creating confidence intervals. 
It is used with \code{plottype="EIS"}, and the default value is \code{alpha=.05}, or with \code{plottype="OCC"}, and the default value is \code{alpha=FALSE}.}

\item{subjects}{
a vector specifying the subjects to plot. 
This argument is only used when \code{plottype="RCC"}.
}

\item{axistype}{
a character string specifying the display variable to be used on the \emph{x}-axis. 
The default is \code{axistype="distribution"}, which uses the subjscoresummary of the distribution specified in thetadist. 
The alternative is \code{axistype="scores"} which displays the expected score.
}

\item{main, xlab, ylab, xlim, ylim, cex}{plotting parameters (see \code{plot()}) with useful defaults.}

\item{...}{further plotting parameters.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
Possible values for plottype are:
\describe{

\item{\code{plottype="density"}}{
produces a simple kernel density plot of the observed scores.
}

\item{\code{plottype="EIS"}}{
plot of the expected item scores for each of the item numbers in the \code{items} argument.
}

\item{\code{plottype="OCC"}}{
plot of the option characteristic curves for each of the item numbers in the \code{items} argument. 
}

\item{\code{plottype="expected"}}{
plot of the observed vs. expected scores.
}

\item{\code{plottype="sd"}}{
plot of the standard deviation of observed scores.
}

\item{\code{plottype="RCC"}}{
plots the RCC and actual score for each subject specified by the \code{subjects} argument.
}

\item{\code{plottype="triangle"}}{
produces a triangle simplex plot with the highest 3 probability options for each item specified by the \code{items} argument. 
}

\item{\code{plottype="tetrahedron"}}{
produces a tetrahedron simplex plot with the highest 4 probability options for each item specified by the \code{items} argument. The tetrahedron plot can be rotated by using the mouse. 
}

\item{\code{plottype="PCA"}}{
produces Principle Component Analysis plot of the test. 
}

Below are values for \code{plottype} used for Differential Item Functioning (DIF) plots. 
They are available only if the \code{groups} argument is specified when creating the \code{ksIRT} object:

\item{\code{plottype="densityDIF"}}{
plots density of observed scores for each of the different groups.
}

\item{\code{plottype="expectedDIF"}}{
plots pairwise expected value comparison plots for each of the different groups.
}

\item{\code{plottype="EISDIF"}}{
plots expected item scores for each of the different groups. 
Accepts the same arguments as \code{plottype="EIS"}, but by default does not show confidence intervals. 
This can be changed with the \code{alpha} argument.
}

\item{\code{plottype="OCCDIF"}}{
plots option characteristic curves for each of the different groups. 
Accepts the same arguments as \code{plottype="OCC"}
}

}


}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

No values are returned from the plot function.
}


\references{
 Mazza A, Punzo A, McGuire B. (2014). 
\pkg{KernSmoothIRT}: An {R} Package for Kernel Smoothing in Item Response Theory.
\emph{Journal of Statistical Software}, 58 6, 1-34. URL: \url{http://www.jstatsoft.org/v58/i06/}.

Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 
   
}




