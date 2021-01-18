\name{cellMap}
\alias{cellMap}

\title{
Draw a cellmap
}

\description{
This function draws a cellmap, possibly of a subset of rows and columns of the data,
and possibly combining cells into blocks. A cellmap shows which cells are missing and which ones are outlying, marking them in red for unusually large cell values and in blue for unusually low cell values. When cells are combined into blocks, the final color is the average of the colors in the individual cells. 
}

\usage{
cellMap(D, R, indcells = NULL, indrows = NULL,
                   standOD=NULL,showVals=NULL,rowlabels="",
                   columnlabels="",mTitle="", rowtitle="",
                   columntitle="",showrows=NULL, showcolumns=NULL,
                   nrowsinblock=1, ncolumnsinblock=1,autolabel=TRUE,
                   columnangle=90,sizetitles=1.1,adjustrowlabels=1,
                   adjustcolumnlabels=1, colContrast=1,outlyingGrad=TRUE,
                   darkestColor = sqrt(qchisq(0.999,1)),
                   drawCircles = TRUE)
}

\arguments{
  \item{D}{
  The data matrix (required input argument).
  }
  \item{R}{
  Matrix of standardized residuals of the cells (required input argument)
  }
  \item{indcells}{
  Indices of outlying cells. Defaults to \code{NULL}, which indicates the cells for which
\eqn{|\code{R}| > \sqrt(qchisq(0.99,1))}.
  }
  \item{indrows}{
  Indices of outlying rows. By default no rows are indicated.
  }
  \item{standOD}{
  Standardized Orthogonal Distance of each row. Defaults to \code{NULL}, then no rows
  are indicated.
  }
  \item{showVals}{
  Takes the values \code{"D"}, \code{"R"} or \code{NULL} and determines whether or not to show the entries of the data matrix (D) or the residuals (R) in the cellmap. Defaults to \code{NULL}, then no values are shown.
  }
  \item{rowlabels}{
  Labels of the rows.
  }
  \item{columnlabels}{
  Labels of the columns.
  }
  \item{mTitle}{
  Main title of the cellMap. 
  }
  \item{rowtitle}{
  Title for the rows.
  }
  \item{columntitle}{
  Title for the columns.
  }
  \item{showrows}{
  Indices of the rows to be shown. Defaults to \code{NULL} which means all rows are shown.
  }
  \item{showcolumns}{
  Indices of the columns to be shown. Defaults to \code{NULL} which means all columns are shown.
  }
  \item{nrowsinblock}{
  How many rows are combined in a block. Defaults to \eqn{1}.
  }
  \item{ncolumnsinblock}{
  How many columns are combined in a block. Defaults to \eqn{1}.
  }
  \item{autolabel}{
  Automatically combines labels of cells in blocks.  
  If \code{FALSE}, you must provide the final columnlabels and/or rowlabels.
  Defaults to \code{TRUE}.
  }
  \item{columnangle}{
  Angle of the column labels. Defaults to \eqn{90}.
  }
  \item{sizetitles}{
  Size of row title and column title. Defaults to \eqn{1.1}.
  }
  \item{adjustrowlabels}{
  Adjust row labels: 0=left, 0.5=centered, 1=right. Defaults to \eqn{1}.
  }
  \item{adjustcolumnlabels}{
  Adjust column labels: 0=left, 0.5=centered, 1=right. Defaults to \eqn{1}.
  }
  \item{colContrast}{
 Parameter regulating the contrast of colors, should be in \eqn{[1,5]}. Defaults to \eqn{1}.
  }
  \item{outlyingGrad}{
  If \code{TRUE}, the color is gradually adjusted in function of the outlyingness. Defaults to \code{TRUE}.
  }
  \item{darkestColor}{
  Standardized residuals bigger than this will get the darkest color.
  }
  \item{drawCircles}{
  Whether or not to draw black circles indicating the outlying rows.
  }
}


\references{
Rousseeuw, P.J., Van den Bossche W. (2018). Detecting Deviating Data Cells. \emph{Technometrics}, \bold{60}(2), 135-145. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2017.1340909}{(link to open access pdf)}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}

\seealso{
\code{\link{DDC}}
}

\examples{
# For examples of the cellmap, we refer to the vignette:
vignette("DDC_examples")
}

