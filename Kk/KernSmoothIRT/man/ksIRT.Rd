\name{ksIRT}
\alias{ksIRT}
\alias{print.ksIRT}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ksIRT - kernel smoothing in Item Response Theory
}
\description{
Fits nonparametric item and options characteristic curves using kernel smoothing techniques.
Within the \pkg{KernSmoothIRT} package, it provides the relevant data for the graphical analysis of multiple choice test and questionnaire data. 
}
\usage{
ksIRT(responses, key, format, kernel = c("gaussian","quadratic","uniform"), itemlabels,
weights,miss = c("option","omit","random.multinom","random.unif"), NAweight = 0, 
evalpoints, nevalpoints, bandwidth = c("Silverman","CV"), RankFun = "sum", SubRank, 
thetadist = list("norm",0,1), groups = FALSE, nsubj)

\method{print}{ksIRT}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{responses}{
input data matrix with options selected from each individual for each item. Rows represent individuals, columns represent items. Alternatively, a data.frame or list can be specified. Missing values are inserted as \code{NA}.
}
\item{key}{
%%     ~~Describe \code{key} here~~
a numeric vector or  a scalar. If \code{key} is a vector, its length must match the number of items; if it is a scalar, its value is used for all items. 

If the items are multiple choice, \code{key} should contain the option that corresponds to the correct response. 

If the data are rating-scale, \code{key} should contain the largest option value for each item. In this case, the weight assigned to each option is equal to its option number. 

More complicated weighting schemes, such as partial credit, can be specified in the \code{weights} argument. If \code{weights} is specified, \code{key} must be left blank. 

}
\item{format}{
a numeric scalar or vector specifying the type of items. 
If all of the items are multiple choice, then \code{format = 1}. 
If all of the items are rating-scale or partial credit, then \code{format = 2}. 
If all of the items are nominal items, then \code{format = 3}. 
If the test has a mixture of items of different formats, then format is a vector with length equal to the number of items with entries of 1 for each multiple choice item and 2 for each rating-scale item. For more complicated weighting schemes use the \code{weights} argument.

}
\item{kernel}{
a character string specifying the kernel function. \code{kernel} must be either \code{"gaussian"}, \code{"quadratic"} or \code{"uniform"}. The default is \code{"gaussian"}.
}
\item{itemlabels}{
optional list of labels for each item. 
If omitted, each item will be labelled according to its numerical order. 
These labels will be used in plotting.
}
\item{weights}{
%%     ~~Describe \code{weights} here~~
optional list that may be used in lieu of including \code{key}. Specifying \code{weights} allows for more complicated weighting schemes than the default. Its length must be equal to the number of items and each entry must be a matrix with option numbers in the first row and option weights in the second row. If weights is omitted and \code{format=1}, then weights are given according to \code{key}. If \code{weights} is omitted and \code{format=2}, then an option weight equals  the option number is given to each response. If \code{weights} is omitted and and \code{format=3}, then weights are set to zero.
}

\item{miss}{
%%     ~~Describe \code{miss} here~~
a character string specifying the method used to manage missing responses.

The default value, \code{miss="option"}, considers the missing responses as a further option, labeled as \code{NA}, with zero weight. 
Such \code{NA} option will be added to the plot of the Option Characteristic Curves. 
Alternatively, a different weight for the \code{NA} option may be specified through the \code{NAweight} argument.

\code{miss="random.unif"} substitutes \code{NA}s with options randomly chosen from the possible ones for the corresponding  item.

\code{miss="random.multinom"} does the same substitution as \code{miss="random.unif"} but each option has a probability of being selected proportional to its relative frequency.

\code{miss="omit"}  excludes from the analysis all the subjects with at least one omitted response.
}

\item{NAweight}{
a scalar value that specifies the weight given to missing responses when \code{miss="option"}. 
The default is zero.
}

\item{evalpoints}{
%%     ~~Describe \code{evalpoints} here~~
an optional numeric vector that specifies the quantiles at which to estimate the Option Characteristic Curves. 
If unspecified, the default is \code{nevalpoints} evenly spaced values with end points determined according to the number of subjects and the distribution specified with the \code{thetadist} argument.

}
	\item{nevalpoints}{
an optional scalar value that specifies the number of evenly spaced points at which curves are estimated. This value is used as an alternative to a user defined vector in the \code{evalpoints} argument. The default value is 51.
The end points are determined according to the number of subjects and to the distribution specified for the \code{thetadist} argument. 
If both \code{nevalpoints} and \code{evalpoints} are specified, then \code{evalpoints} takes precedence.
 
}
\item{bandwidth}{
%%     ~~Describe \code{bandwidth} here~~
either \code{"Silverman"}, \code{"CV"} or a numeric vector specifying, for each item, the bandwidth to use for kernel smoothing. The default value, \code{bandwidth="Silverman"}, is a numeric vector computed following the well-known Silverman's rule of thumb. If \code{bandwidth="CV"}, then the bandwidth is chosen for each item through cross-validation.

 }
 
\item{RankFun}{
a function that is used to rank subjects. The default value is \code{"sum"}. Another common choice is \code{"mean"}.
}

\item{SubRank}{
a numeric vector specifying the rank of each of the subjects. If unspecified and \code{format=1} or \code{format= 2}, subjects will be ranked according to the function passed through the argument \code{RankFun}. When \code{format=3} this argument must be provided.
}
 
\item{thetadist}{
a list specifying the distribution to be used to thetadist (see Ramsay, 1991, p. 615) the subjects. 
By default a standard normal distribution is used. 
A different distribution can be adopted by specifying the first element of the list as \code{"norm"}, \code{"beta"}, \code{"unif"}, \code{"gamma"}, etc. where the character string is the same as used in the subjscoresummary function  \code{qnorm()},\code{qbeta()}, \code{qunif()}, \code{qgamma()}. 
The other elements of the list should be the distribution parameters as required by the subjscoresummary function chosen.
}

\item{groups}{
an optional vector of length equal to the number of subjects containing the group designation of each subject. Adding this option allows for comparisons between groups using the Differential Item Functioning tools (see details section).

}

\item{nsubj}{
an optional numeric value with the number of subjects.
}

\item{x}{
a \code{ksIRT} object to be printed.
}

\item{...}{further parameters}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
When \code{bandwidth="Silverman"}, the rule of thumb of Silverman (1986, p. 45) is implemented with the formula: \code{1.06*sigma.hat*nsubj^(-.2)}, where \code{nsubj} is the number of subjects and \code{sigma.hat} is the standard deviation of the subjscoresummary associated to the subjects according to the distribution specified with \code{thetadist}. 
Note that when \code{thetadist=list("norm",mean,sd)}, \code{sigma.hat} is the value specified for \code{sd}.

Printing the \code{ksIRT} object shows the point polyserial correlation correlation between each item and the overall test score.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

Returned from this function is a \code{ksIRT} object which is a list with the following components:

\item{nitem}{
an integer indicating the number of items.
}

\item{nsubj}{
an integer indicating the number of subjects.
}

\item{nevalpoints}{
an integer indicating the number of points for curve estimation.
}

\item{binaryresp}{
a matrix of binary responses. Each row corresponds to a single option. The first three columns
specify the item, the option, and the corresponding weight. Each additional column is a binary indicator of whether a subject selected that option.
}
\item{OCC}{
a matrix with the first 3 columns the same as \code{binaryresp} and an additional column for each quantile at which the option characteristic curves have been estimated. 
The additional columns contain the kernel smoothed probabilities of selecting each option.
}
\item{stderrs}{
a matrix as \code{OCC} containing the standard errors of \code{OCC}. 
}

\item{subjscore}{
a vector containing the observed score of each subject.
}

\item{itemlabels}{
a list containing the label for each item.
}

\item{thetadist}{
a list indicating the distribution used to rank subjects (see \code{thetadist} in Arguments).
}

\item{subjtheta}{
a vector of quantile ranks for each subject on the distribution specified in \code{thetadist}.
}

\item{evalpoints}{
a vector with the subjscoresummary used in curve estimation.
}

\item{subjscoresummary}{
a vector of subjscoresummary, of probabilities \code{.05}, \code{.25}, \code{.50}, \code{.75}, \code{.95}, for the observed overall scores. 
}

\item{subjscoresummaryevalpoints}{
a vector as \code{subjscoresummary} but computed on \code{subjtheta}.
}

\item{SmthWgts}{
a matrix containing the kernel weights.
}

\item{scale}{
a vector indicating whether each item is multiple-choice, rating-scale or nominal; \code{1} indicates multiple-choice, \code{0} indicates rating-scale, \code{3} indicates nominal.
}
\item{format}{
returns the \code{format} argument passed at function call.
}

\item{bandwidth}{
a vector containing the bandwidths for each item.
}

\item{DIF}{
a list of \code{ksIRT} objects created for each of the subgroups specified by \code{groups}.
}

\item{groups}{
returns the \code{groups} argument passed at function call.
}
\item{itemcor}{
a vector containing the point polyserial correlation for each item.
}
\item{RCC}{
a list of \code{nsubj} vectors containing the normalized likelihood for each value in \code{evalpoints}. 
}
\item{subjthetaML}{
the maximum likelihood estimate for the expected total score of each subject.
}
}

\examples{
 ## Psych101 data
data(Psych101)
Psych1 <- ksIRT(responses = Psychresponses[1:100,], key = Psychkey, format = 1)
Psych1
    
plot(Psych1,plottype="OCC", item=c(24,25,92,96))
plot(Psych1,plottype="EIS", item=c(24,25,92,96))
plot(Psych1, plottype="tetrahedron", items=c(24,92))
plot(Psych1, plottype="triangle", items=c(24,92))
plot(Psych1, plottype="PCA")
plot(Psych1,plottype="RCC", subjects=c(33,92))
 
PCA(Psych1)
subjEIS(Psych1)
subjETS(Psych1)
subjOCC(Psych1, stype="ObsScore")
subjscore(Psych1)
subjthetaML(Psych1)
subjscoreML(Psych1)
 
plot(Psych1, plottype="expected")
plot(Psych1, plottype="sd")
plot(Psych1, plottype="density")

## HIV data
data(HIV)
HIVsubset <- HIV[c(c(1:50),c(1508:1558),c(2934:2984)),]
gr2 <- as.character(HIVsubset$SITE)
DIF2 <- ksIRT(res=HIVsubset[,-(1:3)], key=HIVkey, format = 2, groups=gr2, miss="omit")

plot(DIF2, plottype="expectedDIF", lwd=2)
plot(DIF2, plottype="densityDIF", lwd=2)
plot(DIF2, plottype="EISDIF",  item=c(6,11))

### Ordinal Survey Data
data(BDI)
BDI1 <- ksIRT(responses=BDIresponses, key=BDIkey, format = 2, miss="omit")

plot(BDI1, plottype="OCC", items=1:4)
plot(BDI1, plottype="sd")
plot(BDI1, plottype="density", ylim=c(0,0.1))

}

\references{
 Mazza A, Punzo A, McGuire B. (2014). 
\pkg{KernSmoothIRT}: An {R} Package for Kernel Smoothing in Item Response Theory.
\emph{Journal of Statistical Software}, 58 6, 1-34. URL: \url{http://www.jstatsoft.org/v58/i06/}.

Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 

   
 Silverman, B.W. (1986). 
 \emph{Density Estimation for Statistics and Data Analysis}. 
 Chapman & Hall, London.

}




