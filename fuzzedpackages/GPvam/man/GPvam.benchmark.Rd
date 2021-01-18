\name{GPvam.benchmark}
\alias{GPvam.benchmark}
\docType{data}
\title{
Benchmarks of the program using simulated data.
}
\description{
The data file GPvam.benchmark that is included with the package contains runtime and peak memory requirements for different persistence settings, using simulated data sets with different values for number of years, number of teachers per year, and number of students per teacher. These have been multiplied to show the total number of teachers in the data set, as well as the total number of students. With student.side="R", the persistence="GP" model is most sensitive to increases in the size of the data set. With student.side="G", the memory requirements increase exponentially with the number of students and teachers, and that model should not be considered scalable to extremely large data sets.

All of these benchmarks were performed with Hessian=TRUE. Calculation of the Hessian accounts for anywhere from 20\% to 75\% of those run times. Unless the standard errors of the variance components are needed, leaving Hessian=FALSE will lead to a faster run time with smaller memory requirements.
}
\usage{data(vam_data)}



\examples{
data(GPvam.benchmark)
print(GPvam.benchmark[1,])
}
\keyword{datasets}
