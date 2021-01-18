\name{plot_table}
\alias{plot_table}

\title{
Plot a Table Using Color Intensity for Counts
}

\description{
A table is visualized as a matrix whose cells are shown with intensity of a given color proportional to the count in each cell. The count in a cell must be real: negative numbers or non-integers are acceptable. It provides a global understanding of the underlying pattern.
}

\usage{
plot_table(table, xlab = "Column", ylab = "Row", col = "green3",
           xaxt = "n", yaxt = "n", main = NULL,
           show.value = TRUE, value.cex = 2,
           highlight=c("row.maxima", "none"),
           highlight.col=col,
           mgp=c(0.5,0,0), mar=c(2,2,3,1.5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{A data frame or a matrix.
}
  \item{xlab}{
  The lable of the horizontal axis.
}
  \item{ylab}{
  The lable of the vertical axis.
}
  \item{col}{
  The color corresponding to the maximum value in the table.
}
  \item{xaxt}{
  The style of the horizontal axis. See \code{\link[graphics]{par}}.
}
  \item{yaxt}{
  The style of the vertical axis. See \code{\link[graphics]{par}}.
}
  \item{main}{
  The title of the plot.
}
  \item{show.value}{
  logical. Show the value of each cell in the table on the plot.
}
  \item{value.cex}{
  Relative magnification factor if values are to be put in the cell.
}
  \item{\dots}{
  Parameters acceptable to \code{\link[graphics]{image}} function in the \pkg{graphics} package.
}
  \item{highlight}{
  Specify to highlight row maxima or no highlight.
  When highlighted, a box is placed around each row maximum.
}
  \item{highlight.col}{
  The color used to highlight a cell in the table.
}
  \item{mgp}{
  The margin (in mex units) for the axis title, labels and line. See \code{\link[graphics]{par}}.
}
  \item{mar}{
  The margins of the four sides of the plot. See \code{\link[graphics]{par}}.
}


}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Joe Song
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
opar <- par(mfrow=c(2,2))
plot_table(matrix(1:6, nrow=2), col="seagreen2")

plot_table(matrix(rnorm(20), nrow=5), col="orange", show.value=FALSE)

plot_table(matrix(rpois(16, 2), nrow=4), col="cornflowerblue", highlight="none")

plot_table(matrix(rbinom(15, 8, 0.5), nrow=3), col="sienna2", highlight="none")
par(opar)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
