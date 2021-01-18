\name{s2Data}
\alias{s2Data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Data wrapper for \code{s2net}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function preprocess the data to fit a semi-supervised linear joint trained model.
}
\usage{
s2Data(xL, yL, xU = NULL, preprocess = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xL}{
%%     ~~Describe \code{xL} here~~
  The labeled data. Could be a \code{matrix} or \code{data.frame}. 
}
  \item{yL}{
%%     ~~Describe \code{yL} here~~
  The labels associated with \code{xL}. Could be a \code{vector}, \code{matrix} or \code{data.frame}, of \code{factor} or \code{numeric} types.
}
  \item{xU}{
%%     ~~Describe \code{xU} here~~
The unlabeled data (optional). Could be a \code{matrix} or \code{data.frame}. 
}
  \item{preprocess}{
%%     ~~Describe \code{preprocess} here~~
Should the input data be pre-processed? Possible values are: 

\code{TRUE}{ (default) The data is converted to a matrix. Factor variables are automatically coded using \code{\link{model.matrix}}. The data is scaled, and constant columns are removed.}

\code{FALSE}{ Do nothing. Keep in mind that the theoretical framework assumes that \code{xL} is centered. Unless you are absolutely sure, avoid this.}

Another object of class \code{s2Data} that was obtained from similar data (same original variables). This is useful when using train/validation sets, to apply the validation data the same transformation as train data. 
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an object of S3 class \code{s2Data} with fields

\item{xL }{Transformed labeled data}
\item{yL }{Transformed labels. If \code{yL} was a factor, it is converted to numeric, and the base category is kept in \code{base}}
\item{xU }{Tranformed unlabeled data} 
\item{type }{Type of task. This one is inferred from the response labels.}
\item{base }{Base category for classification \code{0 = base}}

In addition the following attributes are stored.

\item{pr:rm_cols }{logical vector of removed columns}
\item{pr:center }{column center}
\item{pr:scale }{column scale}
\item{pr:ycenter }{yL center. Regression}
\item{pr:yscale }{yL scale. Regression}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Juan C. Laria
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{s2Fista}}
}
\examples{
data("auto_mpg")

train = s2Data( xL = auto_mpg$P1$xL,
                  yL = auto_mpg$P1$yL,
                  xU = auto_mpg$P1$xU,
                  preprocess = TRUE )
show(train)

# Notice how ordered factor variable $cylinders is handled 
# .L (linear) .Q (quadratic) .C (cubic) and .^4
head(train$xL) 


#if you want to do validation with the unlabeled data
idx = sample(length(auto_mpg$P1$yU), 200)

train = s2Data(xL = auto_mpg$P1$xL, yL = auto_mpg$P1$yL, xU = auto_mpg$P1$xU[idx, ])

valid = s2Data(xL = auto_mpg$P1$xU[-idx, ], yL = auto_mpg$P1$yU[-idx], preprocess = train)

test = s2Data(xL = auto_mpg$P1$xU[idx, ], yL = auto_mpg$P1$yU[idx], preprocess = train)

train
valid
test
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~data }% __ONLY ONE__ keyword per line
