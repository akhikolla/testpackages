\name{data_VOC}
\alias{data_VOC}
\docType{data}
\title{
VOC dataset
}
\description{
This dataset contains the data on volatile organic components (VOCs) in urine of children
between 3 and 10 years old. It is composed of pubicly available data from the National
Health and Nutrition Examination Survey (NHANES) and was analyzed in Raymaekers and Rousseeuw (2020). See below for details and references.
}

\usage{data("data_VOC")}

\format{A matrix of dimensions \eqn{512 \times 19}.
The first 16 variables are the VOC, the last 3 are: 
\itemize{
\item \code{SMD460}: number of smokers that live in the same home as the subject
\item \code{SMD470}: number of people that smoke inside the home of the subject
\item \code{RIDAGEYR}: age of the subject
}
Note that the original variable names are kept.
}


\details{
All of the data was collected from the NHANES website, and was part of the 
NHANES 2015-2016 survey. This was the most recent epoch with complete data at the time of extraction. Three datasets were matched in order to assemble this data:
\itemize{
\item UVOC_I:
    contains the information on the Volative organic components in urine
\item DEMO_I:
    contains the demographical information such as age
\item SMQFAM_I:
    contains the data on the smoking habits of family members
}
The dataset was constructed as follows:
\enumerate{
\item Select the relevant VOCs from the UVOC_I data (see column names) and transform by taking the logarithm
\item Match the subjects in the UVOC_I data with their age in the DEMO_I data
\item Select all subjects with age at most 10
\item Match the data on smoking habits with the selected subjects.
}
}

\source{
\url{https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Laboratory&CycleBeginYear=2015}

 \url{https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Demographics&CycleBeginYear=2015}
 
 \url{https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Questionnaire&CycleBeginYear=2015}
}

\references{
J. Raymaekers and P.J. Rousseeuw (2020). Handling cellwise outliers by sparse
regression and robust covariance. \emph{Arxiv: 1912.12446}. \href{https://arxiv.org/abs/1912.12446}{(link to open access pdf)}
}
\examples{
data("data_VOC")
# For an analysis of this data, we refer to the vignette:
vignette("DI_examples")
}
