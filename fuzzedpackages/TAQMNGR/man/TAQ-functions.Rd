%system("R CMD Rdconv --type=html --output=/home/cipollini/local/R/library/TAQMNGR/html/TAQ-functions.html /home/cipollini/Dropbox/francesco/R/MultiPlatform/TAQ/TAQ_R_Project/man/TAQ-functions.Rd")

%%%% Package name
\name{TAQ}

%%%% Aliases
\alias{TAQ.CleanTickByTick}
\alias{TAQ.Aggregate}
\alias{TAQ.Report}
\alias{TAQ.Read}

%\docType{package}

%%%% Title
\title{
  TAQ Manager
}

%%%% Description
\description{
  Manage tick-by-tick transaction data
}


% \author{
% Francesco Calvori <francesco.calvori@gmail.com>\cr
% Fabrizio Cipollini <cipollini.fabrizio@gmail.com>\cr
% Giampiero M. Gallo <giampiero.gallo@gmail.com>\cr
% Maintainer: Francesco Calvori <francesco.calvori@gmail.com>
% }


\usage{
  
  TAQ.CleanTickByTick(dirInput, dirOutput, 
  
    window = 80, deltaTrimmed = 0.10, granularity = 0.04, useCleaned = TRUE)
  
  TAQ.Aggregate(dirInput, symbol, bin, useAggregated = TRUE)
  
  TAQ.Report(dirInput, symbol)

  TAQ.Read(dirInput, symbol, import = NULL, startDate, endDate, bin)
}


\arguments{

  \item{dirInput}{A \emph{character} scalar: the input directory.\cr
  For the function \code{CleanTickByTick()}, \code{dirInput} is the name of the folder containing the raw data files. 
  In this case it is important that \code{dirInput} includes only \code{.gz} files to be cleaned.
  For the remaining functions, \code{dirInput} is the name of the folder including the previously cleaned/aggregated data (appearing as \code{dirOutput} in the function \code{CleanTickByTick()}).}
  \item{dirOutput}{A \emph{character} scalar: the output directory. 
  It must be different from \code{dirInput}.}
  \item{window}{A \emph{numeric} integer scalar: the window size for the trimming procedure of data clean (see Details).}
  
  \item{deltaTrimmed}{A \emph{numeric} scalar into $(0,1)$: the trimming proportion (see Details).}
  
  \item{granularity}{A \emph{numeric} positive scalar: the granularity parameter (see Details).}
  
  \item{useCleaned}{A \emph{logical} scalar: if \code{TRUE}, previously cleaned files (if any) are not cleaned again.}
  
  \item{useAggregated}{A \emph{logical} scalar: if \code{TRUE}, previously aggregated data (if any) are not aggregated again.}
  \item{symbol}{A \emph{character} (vector in \code{TAQ.Aggregate()}; scalar in \code{TAQ.Report()} and \code{TAQ.Read()}): the ticker symbols of interest.}
  
  \item{startDate}{A \emph{numeric} integer scalar: the start date in the \emph{yyyymmdd} format.}
  
  \item{endDate}{A \emph{numeric} integer scalar: the end date in the \emph{yyyymmdd} format.}
  
  \item{bin}{A \emph{numeric} integer scalar: the bin size (in seconds) for aggregating data.}
  
  \item{import}{A \emph{character}: the list of fields to be imported. One or more among:\cr
  \code{"FIRST"}: First price in the bin.\cr
  \code{"MIN"}: Min price in the bin.\cr
  \code{"MAX"}: Max price in the bin.\cr
  \code{"LAST"}: Last price in the bin.\cr
  \code{"SIZE"}: First price in the bin.\cr
  \code{"#TRADES"}: Number of trades in the bin.\cr
  \code{"VWAP"}: Volume Weighted Average Price in the bin.\cr
  If \code{NULL}, all fields are imported.
  }
}

\details{
	The meaning of the arguments \code{window}, \code{deltaTrimmed}, and \code{granularity} is detailed in the reference below.
}

%%\examples{
%% Search where example data are installed
%%x1 <- installed.packages()[, c("Package", "LibPath")]
%%ind <- x1[, "Package"] == "TAQMNGR"
%% Dir
%%dirInput <- paste(x1[ind, "LibPath"], "/data", sep = "")
%%
%%TAQ.CleanTickByTick(dirInput = dirInput, dirOutput = dirInput)
%%}

\examples{
  #### A fake dataset for running the example can be downloaded at 
  #### 'http://local.disia.unifi.it/cipollini/webpage-new/data/data_sample.txt.gz' 
  ## Input
  # dirInput  <- "path of the input folder" 
  # dirOutput <- "path of the output folder" ## Must be different from 'dirInput' 
  ## Clean
  # TAQ.CleanTickByTick(dirInput = dirInput, dirOutput = dirInput)
  ## Make the report (1 at a time)
  # TAQ.Report(dirInput = dirOutput, symbol = c("DOG")) ## A scalar symbol
  # TAQ.Report(dirInput = dirOutput, symbol = c("GNU")) ## A scalar symbol
  ## Aggregate
  # TAQ.Aggregate(dirInput = dirOutput, symbol = c("DOG", "GNU"), bin = 300, 
  #   useAggregated = TRUE)
  ## Import data
  # dog <- TAQ.Read(dirInput = dirOutput, symbol = "DOG", 
  #   startDate = 00010101, endDate = 20141231, bin = 300)
}


\references{
Brownlees, C. T., and Gallo, G. M. (2006). Financial Econometric Analysis at Ultra{--}High Frequency: Data Handling Concerns, \emph{Computational Statistics and Data Analysis} \bold{51}, 2232--2245.
}

% \keyword{ package }
% \seealso{
% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
% }

