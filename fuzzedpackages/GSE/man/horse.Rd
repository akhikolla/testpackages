\name{horse}
\alias{horse}
\docType{data}
\title{Horse-colic data}
\description{
This is a modified version of the original data set (taken from UCI repository, see reference), where only quantitative variables are considered.
This data set is about horse diseases where the task is to determine if the lesion of the horse was surgical or not. It contains rows with completely
missing values except for ID and must be removed by the users. They are kept mainly for pedagogical purposes. 
}
\usage{data(horse)}
\format{
  A data frame with 368 observations on the following 7 variables are quantitative and 1 categorical. The first variable is a numeric id. 

  \tabular{ll}{
	\code{Hospital_Number} \tab numeric id, i.e. the case number assigned to the horse (may not be unique if the horse is treated > 1 time) \cr
	\code{Rectal_temperature}\tab rectal temperature in degree celcius \cr
	\code{Pulse}\tab the heart rate in beats per minute; normal rate is 30-40 for adults \cr
	\code{Respiratory_rate}\tab respiratory rate; normal rate is 8 to 10 \cr
	\code{Nasogastric_reflux_PH} \tab scale is from 0 to 14 with 7 being neutral; normal values are in the 3 to 4 range \cr
	\code{Packed_cell_volume} \tab the number of red cells by volume in the blood; normal range is 30 to 50 \cr
	\code{Total_protein} \tab normal values lie in the 6-7.5 (gms/dL) range \cr
	\code{Abdomcentesis_total_protein} \tab Values are in gms/dL \cr
	\code{surgical_leison} \tab was the problem (lesion) surgical?; 1 = yes, 2 = no \cr
  }
}
\references{
Frank, A. & Asuncion, A. (2010). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. 
Irvine, CA: University of California, School of Information and Computer Science.
}
\source{
  The original data have been taken from the Journal of Statistics Education
  Databases at
  \itemize{
    \item \url{http://archive.ics.uci.edu/ml/datasets/Horse+Colic},
  }
}
\examples{
\dontrun{
data(horse)
horse.cts <- horse[,-c(1,9)] ## remove the id and categorical variable
res <- GSE(horse.cts)
plot(res, which="dd", xlog10=TRUE, ylog10=TRUE)
getOutliers(res)
}
}
\keyword{datasets}
